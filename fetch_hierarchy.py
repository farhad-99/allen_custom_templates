"""
Fetches the ABA_v3 brain atlas region hierarchy and caches it locally.

Tries source endpoints in order:
  1. INCF Scalable Brain Atlas API  (ABA_v3)
  2. Allen Brain Map REST API
  3. IBL Atlas GitHub CSV  (int-brain-lab/iblatlas — 1327 structures, CCF 2017)
  4. CortexLab GitHub CSV  (cortex-lab/allenCCF — same 1327 structures, CCF 2017)

Note on coverage:
  Both GitHub CSV fallbacks cover the Allen CCF v3 *2017* ontology (1327 unique
  grey-matter and white-matter structures, IDs ≤ 1145 in the uint16 range).
  Annotation volumes produced from the Allen CCF *2022* update may contain
  additional structure IDs (typically 1130–1303) that are absent from the 2017
  ontology files.  These extra IDs correspond primarily to new white-matter /
  fiber-tract delineations added in 2022 and will be listed as unmapped in
  the remapping step.

Output: hierarchy_cache.json — list of dicts with keys:
  id (int), acronym (str), name (str), parent_id (int | None)
"""

import argparse
import csv
import io
import json
import os
import sys

import requests

CACHE_FILE = "hierarchy_cache.json"

# INCF SBA endpoints to try
_SBA_URLS = [
    "https://scalablebrainatlas.incf.org/services/sba/ABA_v3/structures",
    "https://scalablebrainatlas.incf.org/services/sba/ABA_v3/ontology",
]

# Allen Brain Map REST API (flat structure list for ontology_id=1 = Mouse CCF)
_ALLEN_API = (
    "https://api.brain-map.org/api/v2/data/Structure/query.json"
    "?criteria=[ontology_id$eq1]&num_rows=all"
    "&only=id,name,acronym,parent_structure_id"
)

# IBL Atlas GitHub CSV — same content as the Allen CCF 2017 ontology but from
# the International Brain Laboratory, a well-maintained authoritative mirror.
_IBLATLAS_CSV = (
    "https://raw.githubusercontent.com/int-brain-lab/iblatlas/main"
    "/iblatlas/allen_structure_tree.csv"
)

# CortexLab GitHub CSV — Allen CCF 2017 structure tree maintained by
# cortex-lab; same 1327 structures as iblatlas.
_CORTEXLAB_CSV = (
    "https://raw.githubusercontent.com/cortex-lab/allenCCF/master"
    "/structure_tree_safe_2017.csv"
)


# ---------------------------------------------------------------------------
# Normalise different source formats to a flat list of dicts
# ---------------------------------------------------------------------------

def _norm_sba(payload):
    items = payload if isinstance(payload, list) else payload.get("structures", [])
    result = []
    for item in items:
        pid = item.get("parent_id") or item.get("parent_structure_id")
        try:
            result.append({
                "id": int(item["id"]),
                "acronym": str(item.get("acronym", "")),
                "name": str(item.get("name", "")),
                "parent_id": int(pid) if pid else None,
            })
        except (KeyError, ValueError):
            continue
    return result or None


def _norm_allen(payload):
    if not payload.get("success"):
        return None
    result = []
    for item in payload.get("msg", []):
        pid = item.get("parent_structure_id")
        try:
            result.append({
                "id": int(item["id"]),
                "acronym": str(item.get("acronym", "")),
                "name": str(item.get("name", "")),
                "parent_id": int(pid) if pid else None,
            })
        except (KeyError, ValueError):
            continue
    return result or None


def _norm_csv(text):
    reader = csv.DictReader(io.StringIO(text))
    result = []
    for row in reader:
        pid = (row.get("parent_structure_id") or "").strip()
        try:
            sid = int(row["id"])
        except (KeyError, ValueError):
            continue
        # Prefer the cleaned 'safe_name'; fall back to 'name'
        name = (row.get("safe_name") or row.get("name") or "").strip()
        result.append({
            "id": sid,
            "acronym": row.get("acronym", "").strip(),
            "name": name,
            "parent_id": int(pid) if pid else None,
        })
    return result or None


# ---------------------------------------------------------------------------
# Fetch helpers
# ---------------------------------------------------------------------------

def _try_sba(session):
    for url in _SBA_URLS:
        try:
            r = session.get(url, timeout=15)
            r.raise_for_status()
            data = _norm_sba(r.json())
            if data:
                print(f"  OK — SBA ({url})")
                return data
        except Exception:
            continue
    return None


def _try_allen(session):
    try:
        r = session.get(_ALLEN_API, timeout=20)
        r.raise_for_status()
        data = _norm_allen(r.json())
        if data:
            print(f"  OK — Allen Brain Map API")
            return data
    except Exception:
        pass
    return None


def _try_github_csv(session, url: str, label: str):
    try:
        r = session.get(url, timeout=30)
        r.raise_for_status()
        data = _norm_csv(r.text)
        if data:
            print(f"  OK — {label} ({url})")
            return data
    except Exception as exc:
        print(f"  FAILED — {label}: {exc}", file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fetch_hierarchy(force: bool = False) -> list:
    """Return the full ABA region hierarchy, fetching and caching if needed.

    Args:
        force: Re-download even if a local cache exists.

    Returns:
        List of dicts: {id, acronym, name, parent_id}
    """
    if not force and os.path.exists(CACHE_FILE):
        print(f"Using cached hierarchy: {CACHE_FILE}")
        with open(CACHE_FILE) as f:
            return json.load(f)

    session = requests.Session()
    session.headers["User-Agent"] = "allen-coarse-parcellation/1.0"

    print("Fetching ABA_v3 region hierarchy …")
    print("  Trying INCF Scalable Brain Atlas API …")
    data = _try_sba(session)

    if data is None:
        print("  Trying Allen Brain Map REST API …")
        data = _try_allen(session)

    if data is None:
        print("  Trying IBL Atlas GitHub CSV …")
        data = _try_github_csv(session, _IBLATLAS_CSV, "IBL Atlas")

    if data is None:
        print("  Trying CortexLab GitHub CSV …")
        data = _try_github_csv(session, _CORTEXLAB_CSV, "CortexLab")

    if data is None:
        raise RuntimeError(
            "All hierarchy sources failed. Check network connectivity."
        )

    print(f"Fetched {len(data)} structures.")
    with open(CACHE_FILE, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Cached → {CACHE_FILE}")
    return data


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument(
        "--force", action="store_true",
        help="Re-fetch even if a local cache already exists"
    )
    args = parser.parse_args()
    fetch_hierarchy(force=args.force)
