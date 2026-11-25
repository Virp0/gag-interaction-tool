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
