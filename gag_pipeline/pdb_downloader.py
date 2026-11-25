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
