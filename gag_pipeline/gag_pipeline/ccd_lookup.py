"""
Chemical Component Dictionary (CCD) lookup
--------------------------------------------
Automatically classifies HETATM residue codes that are NOT already in
gag_detector's curated lists, by querying RCSB's own Chemical Component
Dictionary REST API (https://data.rcsb.org/rest/v1/core/chemcomp/{code}).

This exists because the curated code list in gag_detector.py can only ever
cover GAG variants someone has manually checked. Every time a new synthetic
heparin-analog drug (e.g. a fondaparinux-like pentasaccharide, or a longer
custom oligosaccharide bridge) gets solved and deposited, it introduces
brand-new three-letter component codes that have never been seen before.
Rather than requiring a manual code-list patch for every new structure,
this module asks RCSB's own database what a given code actually IS, and
classifies it the same way a person would by hand:

  1. Is it classified as a saccharide at all? (chem_comp.type contains
     "SACCHARIDE") -- if not, it's not a GAG candidate, full stop
     (rules out buffer salts, metal ions, drugs, crystallization additives,
     etc. without needing a hardcoded exclusion list).
  2. Does its name/formula show heavy sulfation ("sulfo"/"sulfonato" in the
     name, or a bare sulfur atom in the formula) or uronic-acid/amino-sugar
     chemistry ("uronic", "glucosamine", "idopyranuronic", etc. in the
     name)? This is the same signature used to manually confirm Z9L, Z9K,
     GU4, GU6, U9G, YYB, etc. as real GAG units. If so: UNAMBIGUOUS
     evidence on its own (mirrors gag_detector.UNAMBIGUOUS_GAG_CODES).
  3. Otherwise it's a saccharide but not diagnostic on its own (e.g. a
     plain methylated sugar with no sulfation, like GU8/U9D/U9J/U9M in the
     1SR5 16-mer) -- AMBIGUOUS, same treatment as NAG/GAL/MAN/etc.: only
     counted as GAG if PDB LINK records connect it to a confirmed GAG
     residue.

Results are cached to disk (~/.gag_pipeline_ccd_cache.json by default) so
repeated runs, and repeated structures sharing common codes, don't refetch.
If the network is unavailable, lookups fail closed (treated as "unknown",
never silently assumed to be a GAG) and the pipeline falls back to
curated-list-only detection, exactly like before this module existed.
"""
from __future__ import annotations

import json
import os
import re
import urllib.error
import urllib.request
from typing import Callable, Dict, Optional, Set, Tuple

CCD_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/chemcomp/{code}"
DEFAULT_CACHE_PATH = os.path.join(os.path.expanduser("~"), ".gag_pipeline_ccd_cache.json")

_SULFATION_NAME_KEYWORDS = ("SULFO", "SULFAT", "SULPHO", "SULPHAT")
_URONIC_AMINO_SUGAR_KEYWORDS = (
    "URONIC", "GLUCOSAMINE", "GALACTOSAMINE", "IDOPYRANURONIC",
    "GLUCOPYRANURONIC", "IDURONIC", "SULFAMIDO",
)
_SACCHARIDE_TYPE_MARKER = "SACCHARIDE"


_SULFUR_TOKEN_RE = re.compile(r"^S\d*$")


class ChemCompInfo:
    __slots__ = ("code", "type", "name", "formula")

    def __init__(self, code: str, type_: str, name: str, formula: str):
        self.code = code
        self.type = type_ or ""
        self.name = name or ""
        self.formula = formula or ""


def _default_fetch(code: str, timeout: int = 10) -> Optional[ChemCompInfo]:
    """Real network lookup against RCSB's Chemical Component Dictionary."""
    url = CCD_URL_TEMPLATE.format(code=code.upper())
    req = urllib.request.Request(url, headers={"User-Agent": "gag-pipeline/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError):
        return None

    chem_comp = payload.get("chem_comp", {})
    if not chem_comp:
        return None
    return ChemCompInfo(
        code=code.upper(),
        type_=chem_comp.get("type", ""),
        name=chem_comp.get("name", ""),
        formula=chem_comp.get("formula", ""),
    )


def _has_sulfation_evidence(info: ChemCompInfo) -> bool:
    name_upper = info.name.upper()
    if any(kw in name_upper for kw in _SULFATION_NAME_KEYWORDS):
        return True
    for token in info.formula.split():
        if _SULFUR_TOKEN_RE.match(token):
            return True
    return False


def _has_uronic_or_amino_sugar_evidence(info: ChemCompInfo) -> bool:
    name_upper = info.name.upper()
    return any(kw in name_upper for kw in _URONIC_AMINO_SUGAR_KEYWORDS)


def classify_chem_comp(info: ChemCompInfo) -> str:
    """
    Returns one of: "unambiguous", "ambiguous", "not_gag".
    See module docstring for the reasoning behind each tier.
    """
    if _SACCHARIDE_TYPE_MARKER not in info.type.upper():
        return "not_gag"
    if _has_sulfation_evidence(info) or _has_uronic_or_amino_sugar_evidence(info):
        return "unambiguous"
    return "ambiguous"


def _load_cache(cache_path: str) -> Dict[str, dict]:
    if not cache_path or not os.path.exists(cache_path):
        return {}
    try:
        with open(cache_path, "r") as fh:
            return json.load(fh)
    except (ValueError, OSError):
        return {}


def _save_cache(cache_path: str, cache: Dict[str, dict]) -> None:
    if not cache_path:
        return
    try:
        with open(cache_path, "w") as fh:
            json.dump(cache, fh, indent=2, sort_keys=True)
    except OSError:
        pass


def classify_codes(
    codes: Set[str],
    fetch_fn: Optional[Callable[[str], Optional[ChemCompInfo]]] = None,
    cache_path: Optional[str] = DEFAULT_CACHE_PATH,
) -> Tuple[Set[str], Set[str], Dict[str, dict]]:
    """
    Classify a set of unrecognized residue codes against RCSB's Chemical
    Component Dictionary.

    Parameters
    ----------
    codes : set of residue codes not already covered by gag_detector's
        curated UNAMBIGUOUS_GAG_CODES / AMBIGUOUS_SUGAR_CODES sets.
    fetch_fn : injectable lookup function, code -> ChemCompInfo or None.
        Defaults to a real network call; tests supply a stub here so the
        classification logic can be verified without network access.
    cache_path : where to persist lookup results across runs. Pass None to
        disable disk caching (e.g. in tests).

    Returns
    -------
    (newly_unambiguous, newly_ambiguous, classification_log)
        classification_log maps code -> {"tier": ..., "type": ..., "name": ...}
        for every code that was looked up (including "not_gag" and
        "lookup_failed" entries), so callers can report exactly which codes
        were auto-classified and how, for scientific transparency.
    """
    fetch = fetch_fn or _default_fetch
    cache = _load_cache(cache_path) if cache_path else {}

    newly_unambiguous: Set[str] = set()
    newly_ambiguous: Set[str] = set()
    classification_log: Dict[str, dict] = {}
    cache_dirty = False

    for code in sorted(codes):
        code_upper = code.upper()

        if code_upper in cache:
            entry = cache[code_upper]
        else:
            info = fetch(code_upper)
            if info is None:
                entry = {"tier": "lookup_failed", "type": None, "name": None}
            else:
                tier = classify_chem_comp(info)
                entry = {"tier": tier, "type": info.type, "name": info.name}
            cache[code_upper] = entry
            cache_dirty = True

        classification_log[code_upper] = entry
        if entry["tier"] == "unambiguous":
            newly_unambiguous.add(code_upper)
        elif entry["tier"] == "ambiguous":
            newly_ambiguous.add(code_upper)


    if cache_dirty and cache_path:
        _save_cache(cache_path, cache)

    return newly_unambiguous, newly_ambiguous, classification_log
