# GAG-Protein Interaction Pipeline

Automated pipeline that takes a PDB ID, identifies glycosaminoglycan (GAG)
molecules, finds the protein residues that contact them within a 3.2-3.5 A
cutoff, and produces an interaction-region PDB, a CSV interaction table,
sequence and summary reports, and method documentation — with no manual
intervention required beyond supplying the PDB ID.

## Install

```bash
pip install -r requirements.txt
```

(`numpy` is required; `scipy` is optional but strongly recommended — the
neighbor search falls back to a slower pure-Python/NumPy brute-force search
if scipy isn't available.)

## Usage

```bash
# Normal use — downloads directly from RCSB
python -m gag_pipeline.main 1FQ9

# Custom cutoff (must be between 3.2 and 3.5 A) and output directory
python -m gag_pipeline.main 1FQ9 --cutoff 3.2 --out-dir results/1FQ9

# Network-restricted environments: analyze an already-downloaded PDB file
python -m gag_pipeline.main 1FQ9 --local-file /path/to/1FQ9.pdb
```

Or from Python:

```python
from gag_pipeline.main import run_pipeline
result = run_pipeline("1FQ9", out_dir="results/1FQ9", cutoff=3.5)
```

## Output files

For each analyzed structure, written to `--out-dir`:

| File | Contents |
|---|---|
| `<ID>.pdb` (or `.cif`) | Original, unmodified downloaded structure |
| `<ID>_interaction_region.pdb` | GAG residue(s) + interacting protein residues only |
| `<ID>_interactions.csv` | Ordered atom-level interaction table |
| `<ID>_sequences.txt` | Interacting-residue sequences by chain, followed by the Phase 2 conserved-residue (composition) summary: basic / acidic / polar / nonpolar breakdown and the basic-residue (Lys/Arg/His) fraction |
| `<ID>_structure_sequence_table.csv` | Phase 2: one row per interacting residue, linking chain/residue number/3-letter code directly to its position in the extracted sequence and its closest GAG contact distance |
| `<ID>_summary.txt` | Interaction summary report |
| `<ID>_summary.json` | Machine-readable summary |

Open `<ID>_interaction_region.pdb` in Chimera, ChimeraX, or PyMOL to
visually confirm the selected residues correspond to the GAG binding site.

## Cross-structure consensus (Phase 2)

To compare interacting residues across multiple already-analyzed
structures -- finding residue types and residue positions that recur across
them. Structures are matched by UniProt accession (read from the PDB DBREF
lines) and positions are compared in UniProt numbering, so chain letters,
protein-name spelling and residue-number offsets between entries don't matter
(e.g. 1AZX vs 1SR5, both antithrombin). If a file has no UniProt reference it
falls back to exact protein name + chain + residue number:

```bash
python -m gag_pipeline.consensus_main --dir results/
# or list specific summary files directly:
python -m gag_pipeline.consensus_main results/1AZX_summary.json results/1SR5_summary.json
```

Writes `consensus_report.txt` (override with `--out`).

## Module layout

Each module is independently testable and importable on its own:

- `pdb_downloader.py` — fetches structures from RCSB (current PDB format,
  falls back to mmCIF for CIF-only entries)
- `pdb_parser.py` — fixed-column PDB parsing; separates protein / HETATM /
  water / ligand atom populations
- `gag_detector.py` — identifies GAG residues via a curated set of wwPDB
  chemical-component codes (heparin/HS, chondroitin/dermatan, keratan,
  hyaluronan)
- `bio_component_selector.py` — restricts analysis to protein chains
  actually near a GAG (excludes unrelated chains, e.g. crystallization
  companions)
- `neighbor_search.py` — atom-atom Euclidean distance search (scipy
  cKDTree, cutoff 3.2-3.5 A only — 4.5 A is deliberately unsupported)
- `residue_contacts.py` — atom contacts -> unique residue-level contacts;
  builds the local interaction region
- `sequence_extractor.py` — builds interacting-residue sequences per chain
  and cross-checks them against the original structure
- `pdb_exporter.py` — writes the interaction-region PDB subset
- `report_generator.py` — CSV / text / JSON / methods report writers
- `motif_analysis.py` / `motif_extra.py` — Cardin-Weintraub, clusters, CPC clip; turn, vWF, Margalit, clamp
- `sequence_input.py` — sequence-only motif scan (pasted FASTA / sequence)
- `msa_consensus.py` — alignment, consensus, conservation CSV, SVG logo, K/R analysis
- `web/app.py` — Flask + JSmol web interface
- `main.py` — CLI orchestrator tying every module together

## Validation

`tests/make_test_pdb.py` generates a small synthetic structure (protein
chain A with two residues placed in contact range of a two-residue GAG
chain B, an unrelated distant chain C, a water, and a non-GAG ligand) used
to validate the full pipeline end-to-end without requiring network access.
Run it and feed the result through `--local-file` to exercise every stage:

```bash
python tests/make_test_pdb.py
python -m gag_pipeline.main TEST1 --local-file tests/TEST1.pdb --out-dir /tmp/out
```

For real structures, validate against known GAG-bound PDB entries (e.g.
1FQ9, 1HPN, 3IRJ) and visually confirm the exported interaction region in
Chimera/ChimeraX per Task 17 of the spec.

## How each output file is made

11 files per structure (the flowchart is shared by all structures in a folder).

Everything below is produced in one run of `python -m gag_pipeline.main <ID>`.

| File | Made by | How |
|---|---|---|
| `<ID>.pdb` | `pdb_downloader.py` | Downloaded from RCSB, saved unmodified. |
| `<ID>_interactions.csv` | `neighbor_search.py` | Every GAG atom is compared to every protein atom (3D distance). Pairs within the cutoff (3.5 A default, 3.2 A allowed) become rows, sorted by GAG residue then distance. |
| `<ID>_interaction_region.pdb` | `pdb_exporter.py` | Copies the full atoms (same coordinates, chain IDs, residue numbers) of the GAG residues plus every protein residue that had at least one contact. Nothing is recalculated or moved. |
| `<ID>_sequences.txt` | `sequence_extractor.py`, `consensus_analysis.py` | Takes the interacting residues, sorts by chain and number, and joins consecutive numbers into blocks (e.g. Chain A, residues 12-18, one-letter sequence); then counts the interacting residues by type and class (the conserved-residue summary). |
| `<ID>_structure_sequence_table.csv` | `consensus_analysis.py` | One row per interacting residue: chain, number, residue name, one-letter code, its block, position in block, closest GAG distance. |
| `<ID>_motif_report.txt` | `motif_analysis.py` | Phase 3. Slides a window along each chain's full sequence looking for the two Cardin-Weintraub motifs and basic-residue clusters; tests every cation-polar-cation triple among the interacting residues against the CPC clip geometry; then checks which interacting residues fall inside a hit. |
| `<ID>_documentation.txt` | `documentation.py`, `report_generator.py` | Phase 4 and the methods description in one file: results summary with this run's protein name and counts, figure/table captions, methods (how the interactions are calculated), what the three counts mean, and why the cutoff is 3.2-3.5 A. |
| `gag_pipeline_workflow.svg` | `documentation.py` | Phase 4. Workflow flowchart. The same for every structure, so it is written once per output folder (open in a browser or insert into a document). |
| `<ID>_summary.txt` / `.json` | `report_generator.py` | Counts and lists from the steps above, GAG-code decisions, the ChimeraX commands, and the full list of output files. The JSON is what the cross-structure consensus reads. |
| `<ID>_chimerax_commands.cxc` | `chimerax_script.py` | ChimeraX commands that select exactly the GAG and contact residues by chain + residue number, then check coordinates against the exported file. Only written when a GAG was found. |
| `consensus_report.txt` | `consensus_main.py` | Separate command. Compares several `_summary.json` files for shared residues. |

### Phase 3 motif rules (from Cardin & Weintraub 1989)

B = Lys, Arg or His. X = any residue. Only consecutive residue numbers are matched.

- **CW-6** `[X-B-B-X-B-X]`: positions -2, -1, +2 must all be basic.
- **CW-8** `[X-B-B-B-X-X-B-X]`: positions -2 and +3 must be basic; at least one of -3 / -1 must be basic (the paper: "either one or both of these positions are used").
- **Basic cluster**: 3 or more basic residues within any 7-residue window (the paper gives no number; this is our stated choice).

### CPC clip motif (Torrent et al. 2012, PLoS ONE 7:e42692)

One polar residue P (Asn, Gln, Ser, Thr or Tyr) and two cationic residues (Arg or Lys). C is the cation nearer P. Reference C-alpha distances: PC 6.0 +- 1.8, PC' 11.6 +- 1.6, CC' 11.4 +- 2.4 A. A triple matches when all three are within 2.0 SD (adjustable; the paper gives no formal cut-off). Only residues that already contact the GAG are used, and they may be on different chains. Side-chain centroid distances (paper: 6.0, 10.6, 10.7 A) are reported alongside.

Not implemented: cell-penetrating-peptide sequence motifs (no detection rule has been supplied).

## Additional motifs (Phase 3, continued)

Each is its own function in `gag_pipeline/motif_extra.py` with tests in `tests/test_motifs_extra.py`.
Results are appended to the existing `*_motif_report.txt` (no extra files).

| Motif | Source | Rule used |
|---|---|---|
| Turn pattern `TXXBXXTBXXXTBB` | Hileman et al. 1998 | B = Lys/Arg/His; T = residue outside HELIX/SHEET records whose C-alpha(i)-C-alpha(i+3) < 7 A |
| von Willebrand `XBBBXXBBBXXBBX` | Hileman et al. 1998 | sequence only |
| Basic residues ~20 A apart | Margalit et al. 1993 | two GAG-contacting Lys/Arg, C-alpha distance 20 +- 3 A |
| Clamp / parallelogram | Cheng et al. 2016 | 3 or 4 Lys/Arg whose C-alpha atoms are ~5 A apart at ~100 deg (paper values), no Asp/Glu C-alpha within 6 A of the patch. Found on the whole chain; the report says how many residues touch the GAG |

The clamp values (5 A, 100 deg) come from Cheng et al. 2016; the tolerances around them (+-1.5 A, +-25 deg, 4th corner within 2.5 A, acidic radius 6 A) are ours. The clamp was checked on the paper's own examples (CXCL12 2NWG: R41-K24-K43; cardiotoxin A5 1KXI: K24-K37-R38; tests/data). The turn and Margalit numbers are also our assumptions and are printed in every motif report so they can be confirmed.

The motif report ends with **sequence prediction vs structure**: residues inside a Cardin-Weintraub motif, basic cluster,
turn pattern, vWF pattern or clamp compared with residues that actually contact the GAG (precision = overlap/predicted,
recall = overlap/contacting). CPC and Margalit are NOT counted here: they are built from the contacting
residues, so they would match by construction.

## Web interface (Phase 5, JSmol)

```
py -3 -m pip install -r requirements.txt
py -3 -m gag_pipeline.web.app
```
Open http://127.0.0.1:5000. Inputs: PDB ID, uploaded PDB file, one protein sequence, or several sequences
(FASTA). Results page: 3-D view (GAG magenta, interacting residues orange sticks, motif residues cyan cartoon,
toggles for each), residue-numbered sequence with contact/motif highlighting and missing-residue gaps,
clickable residues and motifs that zoom the 3-D view, and downloads of every file (individually or as a zip).

* **Compare structures** tab: 2 to 8 PDB IDs -> a consensus page (residues that touch the GAG in every structure of the
  same protein, matched by UniProt position) plus a 3-D page for each structure. Same as `consensus_main.py`.
* Several FASTA sequences -> alignment (`msa_consensus.py`), consensus, conservation CSV, SVG sequence logo,
  and a Lys/Arg report. Already-aligned FASTA (equal length, with `-`) is used as is; unaligned input is aligned
  by a simple center-star method, which is fine for close homologs only.
* A single sequence -> sequence motifs only (CW-6, CW-8, basic clusters, vWF). Turn and geometric motifs need coordinates.
* Environment variables: `GAG_JOBS_DIR` (result folder, default `./gag_web_jobs`) and `GAG_PDB_CACHE`
  (folder of `<ID>.pdb` files used instead of downloading).
* The web downloads leave out the documentation text and workflow flowchart (the command line still writes them;
  `run_pipeline(..., include_documentation=False)` turns them off).
* JSmol (LGPL, in `gag_pipeline/web/static/jsmol`) runs fully offline in the browser. This server is for local
  use; do not expose it to the internet as is.

## Notes / limitations

- GAG detection uses two tiers of evidence (`gag_detector.py`):
  - `UNAMBIGUOUS_GAG_CODES` — residue codes that are effectively GAG-only
    (natural heparin/HS, chondroitin/dermatan building blocks, plus known
    synthetic/chemically-modified heparin-analog codes such as those used
    in fondaparinux-type pentasaccharides). A code match alone is enough.
  - `AMBIGUOUS_SUGAR_CODES` — codes like NAG/GAL/MAN/FUC/SIA that are real
    GAG building blocks in some contexts but are just as commonly ordinary
    N-/O-linked protein glycosylation sugars. These are only counted as
    GAG if PDB `LINK` records show a sugar-to-sugar bond to an already
    confirmed GAG residue — a lone NAG bonded directly to an Asn side
    chain (glycosylation) is correctly excluded.
  - Both sets are extendable: `gag_detector.UNAMBIGUOUS_GAG_CODES.add("XXX")`
    or `gag_detector.AMBIGUOUS_SUGAR_CODES.add("XXX")`.
  - **Automatic classification (`ccd_lookup.py`):** any residue code not
    covered by the two curated lists above is looked up automatically
    against RCSB's own Chemical Component Dictionary
    (`data.rcsb.org/rest/v1/core/chemcomp/{code}`), classified using the
    same saccharide/sulfation criteria a human reviewer would apply. This
    means a brand-new synthetic GAG-analog code (as happened going from
    1AZX to the longer 1SR5 construct) is picked up automatically instead
    of requiring a manual code-list patch. Results are cached to
    `~/.gag_pipeline_ccd_cache.json` by default. Requires network access to
    `data.rcsb.org`; if unavailable, lookups fail closed and detection
    silently falls back to curated-list-only behavior. Disable entirely
    with `--no-auto-classify` for fully offline/deterministic runs.
  - This was validated against real deposited data from PDB entries 1AZX
    (fondaparinux-type pentasaccharide bound to antithrombin, non-standard
    codes Z9L/Z9K/GU6/GU1/Z9H) and 1SR5 (a longer synthetic 16-mer bridging
    antithrombin and thrombin), both of which also contain N-linked
    glycosylation NAG residues that must not be misidentified as GAG. See
    `tests/make_real_1azx_excerpt.py` and the corresponding tests in
    `tests/test_pipeline.py`.
- The biological-component selector uses a 10 A proximity heuristic to
  drop unrelated chains before the strict 3.2-3.5 A contact search; adjust
  `bio_component_selector.select_relevant_chains(..., proximity_cutoff=...)`
  if a complex has unusually large/small interfaces.
- Downloading requires outbound HTTPS access to `files.rcsb.org` (and, for
  auto-classification, `data.rcsb.org`).

## NMR / multi-model files

Files with several MODEL blocks (NMR ensembles) are analyzed using model 1 only, and the summary says so. In ChimeraX such files open as #1.1, #1.2, ...; the generated `.cxc` script addresses `#1.1` automatically.

## Viewer: functional unit, labels, interaction PDB

The results page shows one functional unit by default: the protein chain with the most
residues touching the GAG, plus the GAG chain it touches most (other chains in the file are
crystal copies). Switch to "Interaction only" (GAG + touching residues) or "All chains".
Protein is grey ribbon; each touching residue gets its own carbon colour with element-coloured
N/O/S; GAG carbons are magenta. The Labels checkbox labels touching residues (name + number)
and GAG sugars (name + number : chain). "Interaction PDB" downloads the coordinates of only
the GAG and the touching residues (`<ID>_interaction_region.pdb`).

## PLIP cross-check (optional)

`pip install plip`, then `python -m gag_pipeline.main 1AZX --plip` (the web app does this automatically when PLIP is installed).
Writes `<ID>_plip_check.txt`: protein residues PLIP links to the GAG vs the residues found by the distance cutoff.

## Test cases from the literature

`python run_test_cases.py` (needs internet) runs 1AZX, 1E03, 1XMN, 1FQ9, 1J8I, 1TB6, 1BFB and 2AXM, scores the contacts
against published binding residues (Mosier 2012, Sarkar & Desai 2015, Torrent 2012), counts every motif and PLIP overlap,
and writes `test_case_output/test_case_results.csv`. Offline check: `--local TEST=tests/REAL_1AZX_EXCERPT.pdb`.
