"""
Builds the leaf → coarse-parent remapping for the 12 target ABA_v3 regions.

Reads:  hierarchy_cache.json  (produced by fetch_hierarchy.py)
Writes: label_mapping.json    {str(original_id): target_id, …}
        label_mapping.csv     original_id, original_acronym, original_name,
                              target_id, target_acronym

For each structure in the hierarchy the script walks up the parent chain and
records the *first* (i.e. lowest / most specific) ancestor that is one of the
12 target regions.  Structures that sit at or above all 12 targets (e.g. root,
grey, cerebrum wrappers) produce no mapping entry; their voxels would be
background in the annotation anyway.
"""

import argparse
import csv
import json
import sys

from fetch_hierarchy import fetch_hierarchy

TARGET_ACRONYMS = [
    "HB", "IB", "MB", "CBN", "CBX", "CNU",
    "HIP", "RHP", "Isocortex", "OLF", "CTXsp", "VS",
]

MAPPING_JSON = "label_mapping.json"
MAPPING_CSV  = "label_mapping.csv"


def _resolve_targets(by_acronym: dict) -> dict:
    """Return {target_id: acronym} for the 12 targets."""
    missing = [a for a in TARGET_ACRONYMS if a not in by_acronym]
    if missing:
        raise ValueError(
            f"Target region(s) not found in hierarchy: {missing}\n"
            "Is the cached hierarchy from a compatible ABA version?"
        )
    return {by_acronym[a]["id"]: a for a in TARGET_ACRONYMS}


def _find_target_ancestor(node_id: int, parent_of: dict, target_ids: set):
    """Walk up from node_id; return the first ancestor in target_ids, or None."""
    visited: set = set()
    current = node_id
    while current is not None:
        if current in visited:
            break
        visited.add(current)
        if current in target_ids:
            return current
        current = parent_of.get(current)
    return None


def build_mapping(force: bool = False) -> tuple[dict, dict]:
    """Build and persist the label remapping.

    Returns:
        (mapping, target_ids)
        mapping:    {original_id (int): target_id (int)}
        target_ids: {target_id (int): acronym (str)}
    """
    structures = fetch_hierarchy()

    by_id      = {s["id"]: s for s in structures}
    by_acronym = {s["acronym"]: s for s in structures}

    target_ids = _resolve_targets(by_acronym)
    target_id_set = set(target_ids)

    print(f"\nTarget regions ({len(target_ids)}):")
    for tid in sorted(target_ids, key=lambda x: target_ids[x]):
        print(f"  {target_ids[tid]:12s}  id={tid}")

    parent_of = {s["id"]: s["parent_id"] for s in structures}

    mapping: dict[int, int] = {}
    above_targets: list[str] = []   # structures that are ancestors of targets
    no_path: list[str] = []          # disconnected / orphan nodes

    for s in structures:
        sid = s["id"]
        ancestor = _find_target_ancestor(sid, parent_of, target_id_set)
        if ancestor is not None:
            mapping[sid] = ancestor
        else:
            acro = s["acronym"]
            # Determine why there's no mapping
            if s["parent_id"] is None:
                # Root node — expected to be unmapped
                pass
            else:
                above_targets.append(acro)

    if above_targets:
        print(
            f"\nNOTE: {len(above_targets)} structure(s) lie above all 12 targets "
            f"and are not mapped (expected for high-level nodes):"
        )
        print("  " + ", ".join(above_targets[:30]) +
              ("  …" if len(above_targets) > 30 else ""))

    # Persist JSON
    with open(MAPPING_JSON, "w") as f:
        json.dump({str(k): v for k, v in mapping.items()}, f, indent=2)

    # Persist CSV (sorted for easy diffing)
    with open(MAPPING_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "original_id", "original_acronym", "original_name",
            "target_id", "target_acronym",
        ])
        for sid in sorted(mapping):
            s = by_id.get(sid, {})
            writer.writerow([
                sid,
                s.get("acronym", ""),
                s.get("name", ""),
                mapping[sid],
                target_ids[mapping[sid]],
            ])

    print(f"\nMapping saved:")
    print(f"  {MAPPING_JSON}  ({len(mapping)} entries)")
    print(f"  {MAPPING_CSV}")

    # Per-target counts
    from collections import Counter
    counts = Counter(mapping.values())
    print(f"\nLeaf structures per target region:")
    for tid in sorted(target_ids, key=lambda x: target_ids[x]):
        print(f"  {target_ids[tid]:12s}  {counts.get(tid, 0):4d} structures mapped")

    return mapping, target_ids


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument(
        "--force", action="store_true",
        help="Re-fetch hierarchy even if cached"
    )
    args = parser.parse_args()
    build_mapping(force=args.force)
