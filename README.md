# ABA_v3 Coarse Parcellation Pipeline

Collapses a fine-grained Allen Brain Atlas v3 (ABA_v3) NIfTI label image into
12 top-level grey-matter regions by walking each leaf region up the ABA
hierarchy until one of the target ancestors is reached.

## Target regions

| Acronym | Allen ID | Name |
|---------|----------|------|
| VS | 73 | Ventricular systems |
| MB | 313 | Midbrain |
| Isocortex | 315 | Isocortex |
| CBN | 519 | Cerebellar nuclei |
| CBX | 528 | Cerebellar cortex |
| CNU | 623 | Cerebral nuclei |
| OLF | 698 | Olfactory areas |
| CTXsp | 703 | Cortical subplate |
| RHP | 822 | Retrohippocampal region |
| HB | 1065 | Hindbrain |
| HIP | 1080 | Hippocampal region |
| IB | 1129 | Interbrain |

## Requirements

[Pixi](https://pixi.sh) — dependency and task manager.

```
curl -fsSL https://pixi.sh/install.sh | bash
```

All Python dependencies (`numpy`, `nibabel`, `requests`) are declared in
`pixi.toml` and installed automatically into an isolated environment.

## Usage

```bash
# Full pipeline: fetch hierarchy → build mapping → remap labels
pixi run remap path/to/annotation.nii.gz
```

This writes `<input_stem>_coarse12.nii.gz` alongside the input file and prints
a per-region voxel count table.

Individual stages can also be run separately:

```bash
pixi run fetch-hierarchy          # download + cache ABA hierarchy
pixi run build-mapping            # build leaf→target remapping
pixi run python remap_labels.py path/to/annotation.nii.gz
```

Force a fresh hierarchy download (ignores the cache):

```bash
pixi run python fetch_hierarchy.py --force
```

## Pipeline stages

### 1. `fetch_hierarchy.py`

Downloads the complete ABA_v3 region hierarchy and caches it to
`hierarchy_cache.json`. Tries sources in order:

1. INCF Scalable Brain Atlas API (`scalablebrainatlas.incf.org`)
2. Allen Brain Map REST API (`api.brain-map.org`)
3. IBL Atlas GitHub CSV (`int-brain-lab/iblatlas`)
4. CortexLab GitHub CSV (`cortex-lab/allenCCF`)

### 2. `build_mapping.py`

Reads `hierarchy_cache.json` and, for every structure ID, walks up the parent
chain until the lowest qualifying ancestor among the 12 targets is found.
Exports:

- `label_mapping.json` — `{"original_id": target_id, …}`
- `label_mapping.csv` — human-readable table with acronyms and names

### 3. `remap_labels.py`

Loads the input NIfTI, applies the mapping via a numpy LUT, validates that
every non-zero output voxel belongs to one of the 12 targets, and saves the
remapped image. Affine transform and header metadata are copied verbatim from
the input.

## Included outputs (P56 annotation)

Pre-computed outputs for `atlas/P56_Annotation.nii.gz` are committed to the
repository:

| File | Description |
|------|-------------|
| `atlas/P56_Annotation_coarse12.nii.gz` | 12-region label image |
| `label_mapping.json` | Full leaf→target ID mapping (1105 entries) |
| `label_mapping.csv` | Same mapping with acronyms and names |
| `hierarchy_cache.json` | Cached ABA_v3 hierarchy (1328 structures) |

### Voxel counts

| Region | Allen ID | Voxels |
|--------|----------|--------|
| VS | 73 | 278,223 |
| MB | 313 | 1,222,635 |
| Isocortex | 315 | 8,933,121 |
| CBN | 519 | 29,826 |
| CBX | 528 | 351,145 |
| CNU | 623 | 1,188,003 |
| OLF | 698 | 2,609,133 |
| CTXsp | 703 | 241,555 |
| RHP | 822 | 1,248,457 |
| HB | 1065 | 2,972,674 |
| HIP | 1080 | 887,607 |
| IB | 1129 | 4,287,123 |
| **Total** | | **24,249,502** |

## Notes on unmapped voxels

~24.6% of labeled voxels in the P56 annotation do not map to any of the 12
targets. These fall into three categories, all warned about explicitly at
runtime:

1. **White matter / fiber tracts (106 IDs)** — e.g. corpus callosum,
   spinocerebellar tract, cranial nerves. These are legitimately outside all 12
   grey-matter regions and correctly become background (0).

2. **Intermediate grey-matter wrapper nodes (10 IDs)** — e.g. `grey`,
   `Cerebellum` (CB), `Cerebrum` (CH), `Cortical plate` (CTXpl). These sit
   *above* the 12 targets in the hierarchy and cannot be unambiguously assigned
   to a single target region.

3. **CCF 2022 additions (73 IDs, range 1130–1303)** — structure IDs present in
   the annotation volume but absent from the publicly available CCF 2017
   ontology files. These are likely white-matter delineations added in the
   Allen CCF 2022 update. They will be resolved automatically once a public
   mirror of the 2022 structure tree becomes available; to add support sooner,
   replace `hierarchy_cache.json` with a structure tree that covers these IDs
   and re-run `pixi run build-mapping`.
