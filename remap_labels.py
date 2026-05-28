"""
=============================================================================
ABA_v3 Label Image Coarse Remapping Pipeline
=============================================================================

OVERVIEW
--------
This script is stage 3 of a 3-stage pipeline that collapses a fine-grained
Allen Brain Atlas v3 (ABA_v3) parcellation NIfTI into a coarser 12-region
parcellation.  Each non-zero voxel value (an Allen structure ID) is replaced
by the Allen ID of its lowest qualifying ancestor in the ABA_v3 hierarchy.

PIPELINE STAGES
---------------
1. fetch_hierarchy.py
   Downloads the full ABA_v3 region hierarchy from the INCF Scalable Brain
   Atlas API (primary source).  If that is unreachable it falls back to the
   Allen Brain Map REST API, then to a public GitHub CSV mirror.  Result is
   cached to hierarchy_cache.json.

2. build_mapping.py
   Reads hierarchy_cache.json.  For every region ID, walks up the parent
   chain until one of the 12 target regions is reached; records the *lowest*
   (most specific) matching ancestor.  Exports label_mapping.json and
   label_mapping.csv.

3. remap_labels.py  ← this file
   Loads the input NIfTI, applies the label LUT from label_mapping.json via
   numpy, validates that every non-zero output voxel belongs to one of the
   12 targets, prints a per-region voxel-count summary, and writes the
   remapped image alongside the input.

TARGET REGIONS (12)
-------------------
  HB        Hindbrain
  IB        Interbrain
  MB        Midbrain
  CBN       Cerebellar nuclei
  CBX       Cerebellar cortex
  CNU       Cerebral nuclei
  HIP       Hippocampal region
  RHP       Retrohippocampal region
  Isocortex Isocortex
  OLF       Olfactory areas
  CTXsp     Cortical subplate
  VS        Ventricular systems

Output label values are the Allen structure IDs of the 12 target regions so
that the image retains meaningful, look-up-able identifiers.

USAGE
-----
  # Full pipeline (fetch → build → remap):
  pixi run remap <path/to/annotation.nii.gz>

  # This script directly (after fetch + build have run):
  pixi run python remap_labels.py <path/to/annotation.nii.gz>

OUTPUT
------
  <input_stem>_coarse12.nii.gz — written to the same directory as the input.
  Affine transform and NIfTI header metadata are copied verbatim from the
  input image.
=============================================================================
"""

import argparse
import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

MAPPING_JSON = "label_mapping.json"
CACHE_JSON   = "hierarchy_cache.json"

TARGET_ACRONYMS = [
    "HB", "IB", "MB", "CBN", "CBX", "CNU",
    "HIP", "RHP", "Isocortex", "OLF", "CTXsp", "VS",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_mapping() -> dict[int, int]:
    """Load label_mapping.json → {original_id: target_id}."""
    if not Path(MAPPING_JSON).exists():
        sys.exit(
            f"ERROR: {MAPPING_JSON} not found.\n"
            "Run  pixi run build-mapping  (or  pixi run remap <input>  which\n"
            "runs the full pipeline automatically)."
        )
    with open(MAPPING_JSON) as f:
        raw = json.load(f)
    return {int(k): int(v) for k, v in raw.items()}


def _load_hierarchy_index() -> tuple[dict[int, str], dict[int, int | None], set[int]]:
    """Read hierarchy_cache.json.

    Returns:
        target_id_to:  {target_id: acronym}
        parent_of:     {id: parent_id}  for every structure in the cache
        fiber_tract_ids: set of IDs that sit in the 'fiber tracts' branch
    """
    if not Path(CACHE_JSON).exists():
        sys.exit(
            f"ERROR: {CACHE_JSON} not found.\n"
            "Run  pixi run fetch-hierarchy  first."
        )
    with open(CACHE_JSON) as f:
        structures = json.load(f)

    by_acronym  = {s["acronym"]: s for s in structures}
    parent_of   = {s["id"]: s["parent_id"] for s in structures}
    by_id       = {s["id"]: s for s in structures}

    target_id_to = {by_acronym[a]["id"]: a for a in TARGET_ACRONYMS if a in by_acronym}

    # Build the set of IDs that live in the fiber-tracts / white-matter branch.
    # The Allen ontology places fiber tracts under a top-level node called
    # 'fiber tracts' (id=1009) alongside 'grey' (id=8) and 'VS' (id=73).
    fiber_root = by_acronym.get("fiber tracts", by_acronym.get("ft", {})).get("id")
    if fiber_root is None:
        # fallback: find via name
        for s in structures:
            if "fiber tract" in s.get("name", "").lower():
                fiber_root = s["id"]
                break

    fiber_tract_ids: set[int] = set()
    if fiber_root:
        # BFS/DFS over children
        stack = [fiber_root]
        while stack:
            nid = stack.pop()
            fiber_tract_ids.add(nid)
            for s in structures:
                if s["parent_id"] == nid:
                    stack.append(s["id"])

    return target_id_to, parent_of, fiber_tract_ids


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def remap(input_path: str) -> None:
    input_path = Path(input_path)
    if not input_path.exists():
        sys.exit(f"ERROR: Input file not found: {input_path}")

    mapping                             = _load_mapping()
    target_id_to, parent_of, fiber_ids = _load_hierarchy_index()
    target_ids                          = set(target_id_to.keys())
    tree_ids                            = set(parent_of.keys())

    if len(target_ids) != len(TARGET_ACRONYMS):
        missing = set(TARGET_ACRONYMS) - set(target_id_to.values())
        print(
            f"WARNING: {len(missing)} target acronym(s) not found in "
            f"hierarchy cache: {missing}",
            file=sys.stderr,
        )

    # ---- Load NIfTI -------------------------------------------------------
    print(f"Loading: {input_path}")
    img  = nib.load(str(input_path))
    data = np.asarray(img.dataobj, dtype=np.int32)
    print(f"  Shape: {data.shape}  dtype: {img.header.get_data_dtype()}")

    nonzero_mask   = data != 0
    nonzero_labels = np.unique(data[nonzero_mask]).tolist()
    print(f"  Unique non-zero labels in input: {len(nonzero_labels)}")

    # ---- Coverage check ---------------------------------------------------
    unmapped = [l for l in nonzero_labels if l not in mapping]
    if unmapped:
        fiber_unmapped    = [l for l in unmapped if l in fiber_ids]
        grey_intermediate = [l for l in unmapped if l in tree_ids and l not in fiber_ids]
        unknown_ids       = [l for l in unmapped if l not in tree_ids]

        print(
            f"\n{'!'*62}\n"
            f"WARNING: {len(unmapped)} input label ID(s) are NOT in the mapping\n"
            f"and will be set to background (0) in the output.\n"
            f"\n"
            f"  White matter / fiber tracts ({len(fiber_unmapped)} IDs):\n"
            f"    These are legitimately outside all 12 grey-matter targets.\n"
            f"    IDs: {fiber_unmapped[:20]}"
            + ("  …" if len(fiber_unmapped) > 20 else "") + "\n"
            f"\n"
            f"  Intermediate grey-matter nodes ({len(grey_intermediate)} IDs):\n"
            f"    These sit above the 12 targets in the hierarchy (e.g.\n"
            f"    'grey', 'Cerebellum', 'Cerebrum') and cannot be unambiguously\n"
            f"    assigned to a single target region.\n"
            f"    IDs: {grey_intermediate}\n"
            f"\n"
            f"  Unknown IDs — not in the CCF 2017 ontology"
            f" ({len(unknown_ids)} IDs):\n"
            f"    These are likely white-matter structures added in the Allen\n"
            f"    CCF 2022 update, which is not yet publicly mirrored on GitHub.\n"
            f"    IDs: {unknown_ids[:30]}"
            + ("  …" if len(unknown_ids) > 30 else "") + "\n"
            f"{'!'*62}\n",
            file=sys.stderr,
        )

    # ---- Build LUT and apply ----------------------------------------------
    # Use a flat numpy array indexed by label value for speed.
    max_label = max(nonzero_labels) if nonzero_labels else 0
    lut       = np.zeros(max_label + 1, dtype=np.int32)
    for orig, tgt in mapping.items():
        if 0 < orig <= max_label:
            lut[orig] = tgt

    print("Applying label remapping …")
    # Clip ensures we never index outside lut; out-of-range values → 0.
    remapped = np.where(
        (data > 0) & (data <= max_label),
        lut[np.clip(data, 0, max_label)],
        0,
    ).astype(np.int32)

    # ---- Validation -------------------------------------------------------
    out_nonzero = remapped[remapped != 0]
    unexpected  = [int(v) for v in np.unique(out_nonzero) if int(v) not in target_ids]
    assert not unexpected, (
        f"Validation FAILED: output contains unexpected label IDs: {unexpected}"
    )

    # ---- Summary table ----------------------------------------------------
    col_w = max(len(a) for a in TARGET_ACRONYMS) + 2
    hdr   = f"{'Region':<{col_w}} {'Label ID':>10}  {'Voxels':>12}"
    sep   = "-" * len(hdr)
    print(f"\n{hdr}\n{sep}")
    total = 0
    for tid in sorted(target_ids):
        acro  = target_id_to[tid]
        count = int(np.sum(remapped == tid))
        total += count
        print(f"  {acro:<{col_w-2}} {tid:>10}  {count:>12,}")
    print(sep)
    print(f"  {'TOTAL':<{col_w-2}} {'':>10}  {total:>12,}")
    pct = 100.0 * total / remapped.size
    print(f"\n  Brain coverage: {pct:.1f}% of total voxels")

    if unmapped:
        lost = int(np.sum(np.isin(data, unmapped)))
        pct_labeled = 100.0 * lost / int(nonzero_mask.sum())
        print(
            f"  Voxels from unmapped IDs → background: {lost:,}  "
            f"({pct_labeled:.1f}% of labeled voxels; "
            f"primarily white matter / fiber tracts)",
            file=sys.stderr,
        )

    # ---- Save output ------------------------------------------------------
    stem     = input_path.name
    for ext in (".nii.gz", ".nii"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    out_path = input_path.parent / f"{stem}_coarse12.nii.gz"

    out_img = nib.Nifti1Image(remapped, img.affine, img.header)
    out_img.header.set_data_dtype(np.int32)
    nib.save(out_img, str(out_path))
    print(f"\nSaved: {out_path}")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Remap ABA_v3 fine-grained NIfTI labels to 12 coarse regions.\n"
            "Expects hierarchy_cache.json and label_mapping.json in the current\n"
            "directory (produced by fetch_hierarchy.py and build_mapping.py)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        help="Path to the input ABA_v3 NIfTI label image (.nii or .nii.gz)",
    )
    args = parser.parse_args()
    remap(args.input)
