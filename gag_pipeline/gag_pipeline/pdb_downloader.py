"""
Structure downloader
---------------------
Task 3: Automatically download a structure for a given PDB ID.

Supports the current RCSB "current format" .pdb endpoint and falls back to
mmCIF and the legacy-format endpoint when a plain PDB file is unavailable
(e.g. for large structures RCSB no longer serves in legacy PDB format).

Network note: this module needs outbound access to files.rcsb.org. If run in
a network-restricted environment, download_structure() will raise
ConnectionError with a clear message; supply a local file path to
pdb_parser.parse_pdb() directly in that case.
"""
from __future__ import annotations

import os
import urllib.request
import urllib.error

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"

RCSB_LEGACY_URL = "https://files.rcsb.org/download/{pdb_id}-pdb-bundle.zip"


class DownloadError(RuntimeError):
    pass


def _fetch(url: str, dest_path: str, timeout: int = 30) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "gag-pipeline/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    with open(dest_path, "wb") as fh:
        fh.write(data)


def download_structure(pdb_id: str, out_dir: str, prefer_format: str = "pdb") -> str:
    """
    Download the original, unmodified structure for `pdb_id` into `out_dir`.

    Parameters
    ----------
    pdb_id : str
        Four-character PDB identifier, e.g. "1FQ9".
    out_dir : str
        Directory to store the downloaded file. Created if missing.
    prefer_format : str
        "pdb" (default) tries the legacy/current .pdb format first, then
        falls back to mmCIF (.cif) if the entry is CIF-only (common for
        very large structures).

    Returns
    -------
    str
        Path to the downloaded, unmodified structure file.

    Raises
    ------
    DownloadError
        If neither format could be retrieved.
    """
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4:
        raise ValueError(f"'{pdb_id}' does not look like a 4-character PDB ID")

    os.makedirs(out_dir, exist_ok=True)

    attempts = []
    if prefer_format == "pdb":
        attempts = [
            (RCSB_PDB_URL.format(pdb_id=pdb_id), os.path.join(out_dir, f"{pdb_id}.pdb")),
            (RCSB_CIF_URL.format(pdb_id=pdb_id), os.path.join(out_dir, f"{pdb_id}.cif")),
        ]
    else:
        attempts = [
            (RCSB_CIF_URL.format(pdb_id=pdb_id), os.path.join(out_dir, f"{pdb_id}.cif")),
            (RCSB_PDB_URL.format(pdb_id=pdb_id), os.path.join(out_dir, f"{pdb_id}.pdb")),
        ]

    last_err = None
    for url, dest in attempts:
        try:
            _fetch(url, dest)
            return dest
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_err = e
            continue

    raise DownloadError(
        f"Could not download structure '{pdb_id}' in any supported format "
        f"(tried PDB and mmCIF at files.rcsb.org). Last error: {last_err}"
    )
