# GAG-protein interaction tool

Takes a PDB ID or a PDB file, finds the glycosaminoglycan (GAG) chains, lists the protein residues within a distance cutoff (3.2-3.5 A, default 3.5 A) of each sugar, scans for published heparin-binding patterns, and shows one functional unit in a JSmol viewer.

## Install

```
pip install -r requirements.txt
```

Optional: `pip install plip` for the PLIP cross-check.

## Run

```
python -m gag_pipeline.main 1AZX
python -m gag_pipeline.main 1AZX --cutoff 3.3 --out-dir results --plip
python -m gag_pipeline.main MYID --local-file my_structure.pdb
python -m gag_pipeline.consensus_main 1AZX 1SR5
python -m gag_pipeline.web.app
python -m pytest -q
python run_test_cases.py
```

The web app is at http://127.0.0.1:5000. `--no-auto-classify` skips the RCSB Chemical Component Dictionary lookup for unknown sugar codes.

## Code files

| File | Purpose |
|---|---|
| `main.py` | Runs the whole analysis and the command line |
| `pdb_downloader.py` | Downloads a structure from RCSB |
| `pdb_parser.py` | Reads ATOM, HETATM, LINK, HELIX, SHEET and DBREF records |
| `gag_detector.py` | Finds GAG sugar residues |
| `ccd_lookup.py` | Classifies unknown sugar codes with the RCSB Chemical Component Dictionary |
| `bio_component_selector.py` | Keeps the protein chains within 10 A of a GAG |
| `neighbor_search.py` | Finds protein atoms within the cutoff of GAG atoms |
| `residue_contacts.py` | Groups atom contacts into residue contacts |
| `sequence_extractor.py` | Builds the sequences of the interacting residues |
| `motif_analysis.py` | Cardin-Weintraub motifs, basic clusters, CPC clip |
| `motif_extra.py` | Turn, von Willebrand, Margalit and clamp patterns, and pattern-versus-contact statistics |
| `plip_check.py` | Compares the contacts with PLIP |
| `pdb_exporter.py` | Writes the interaction-region PDB |
| `report_generator.py` | Builds the report text and the contact CSV |
| `consensus_analysis.py` | Residue table, residue-type counts, and comparison across structures |
| `consensus_main.py` | Command line for the cross-structure comparison |
| `msa_consensus.py` | Alignment, consensus, conservation and sequence logo |
| `sequence_input.py` | Motif scan of a pasted sequence |
| `documentation.py` | Writes the methods text and workflow flowchart |
| `web/app.py` | Flask web interface |
| `web/templates/` | Web pages, including the JSmol viewer |
| `run_test_cases.py` | Runs the literature test cases and writes a results table |
| `tests/` | Unit tests and small test structures |

## Output files

Per structure, in `--out-dir`:

| File | Contents |
|---|---|
| `<ID>.pdb` | The structure as downloaded |
| `<ID>_interaction_region.pdb` | The functional unit's GAG chain and the protein residues touching it, nothing else |
| `<ID>_interactions.csv` | One row per atom pair within the cutoff, with both atoms' coordinates and the distance |
| `<ID>_structure_sequence_table.csv` | One row per touching residue: chain, number, residue name, position in its sequence block, closest distance to the GAG |
| `<ID>_report.txt` | Summary, sequences of the touching stretches, and every motif found |
| `<ID>_plip_check.txt` | With `--plip`: the PLIP comparison |
| `<ID>_documentation.txt`, `gag_pipeline_workflow.svg` | Command line only: methods text and flowchart |

Files with several MODEL blocks (NMR) are analyzed with model 1 only.

## GAG detection

A residue is a GAG sugar if its code is on the curated list, or, for ambiguous sugar codes (NAG, GAL, MAN and similar), if LINK records join it to a confirmed GAG residue. Codes on neither list are looked up in the RCSB Chemical Component Dictionary (needs network access to data.rcsb.org). The rules are heuristic and not benchmarked.

## Limitations

- Only model 1 of NMR files is analyzed.
- All alternate locations of an atom are counted, so the atom-pair count can be higher than the number of atoms.
- Modified residues in a fixed list (MSE, SEP, TPO, PTR, HYP and others) count as protein. Nucleic acids do not.
- The viewer groups residues by chain and number, so two residues that differ only by insertion code are selected together. The tables and CSVs keep them apart.
- Pattern statistics (precision and recall) have no published benchmark. Thresholds we chose ourselves are marked in the report.
- Contacts across crystal symmetry mates are not included.

## Functional unit

The viewer and the page counts use one functional unit: the protein chain with the most contacting residues and the GAG chain it touches most. Other chains in the file are crystal copies. The contact CSV, the report and the sequence table count every chain in the file. `<ID>_interaction_region.pdb` holds only the functional unit. Views: functional unit, interaction only, all chains.

## Motif rules (Cardin & Weintraub 1989)

B = Lys, Arg or His. X = any residue. Only consecutive residue numbers are matched.

- **CW-6** `[X-B-B-X-B-X]`: positions -2, -1, +2 must all be basic.
- **CW-8** `[X-B-B-B-X-X-B-X]`: positions -2 and +3 must be basic; at least one of -3 / -1 must be basic (the paper: "either one or both of these positions are used").
- **Basic cluster**: 3 or more basic residues within any 7-residue window (the paper gives no number; this is our stated choice).

### CPC clip motif (Torrent et al. 2012, PLoS ONE 7:e42692)

One polar residue P (Asn, Gln, Ser, Thr or Tyr) and two cationic residues (Arg or Lys). C is the cation nearer P. Reference C-alpha distances: PC 6.0 +- 1.8, PC' 11.6 +- 1.6, CC' 11.4 +- 2.4 A. A triple matches when all three are within 2.0 SD (adjustable; the paper gives no formal cut-off). Only residues that already contact the GAG are used, and they may be on different chains. Side-chain centroid distances (paper: 6.0, 10.6, 10.7 A) are reported alongside.

Not implemented: cell-penetrating-peptide sequence motifs (no detection rule has been supplied).

### Additional motifs

Each is its own function in `gag_pipeline/motif_extra.py` with tests in `tests/test_motifs_extra.py`.
Results are in the motif section of `<ID>_report.txt`.

| Motif | Source | Rule used |
|---|---|---|
| Turn pattern `TXXBXXTBXXXTBB` | Hileman et al. 1998 | B = Lys/Arg/His/Gln; T = Gly/Pro/Asn/Ser/Asp/Thr, or a residue outside HELIX/SHEET records whose C-alpha(i)-C-alpha(i+3) < 7 A |
| von Willebrand `XBBBXXBBBXXBBX` | Hileman et al. 1998 | sequence only |
| Basic residues ~20 A apart | Margalit et al. 1993 | two GAG-contacting Lys/Arg, C-alpha distance 20 +- 3 A |
| Clamp / parallelogram | Cheng et al. 2016 | 3 or 4 Lys/Arg whose C-alpha atoms are ~5 A apart at ~100 deg (paper values), no Asp/Glu C-alpha within 6 A of the patch. Found on the whole chain; the report says how many residues touch the GAG |

The clamp values (5 A, 100 deg) come from Cheng et al. 2016; the tolerances around them (+-1.5 A, +-25 deg, 4th corner within 2.5 A, acidic radius 6 A) are ours. The clamp was checked on the paper's own examples (CXCL12 2NWG: R41-K24-K43; cardiotoxin A5 1KXI: K24-K37-R38; tests/data). The turn and Margalit numbers are also our assumptions and are printed in every motif report so they can be confirmed.

The motif report ends with **sequence prediction vs structure**: residues inside a Cardin-Weintraub motif, basic cluster,
turn pattern, vWF pattern or clamp compared with residues that actually contact the GAG (precision = overlap/predicted,
recall = overlap/contacting). CPC and Margalit are NOT counted here: they are built from the contacting
residues, so they would match by construction.

## Web interface

Inputs: PDB ID, PDB file upload, several PDB IDs to compare, or a pasted sequence (several FASTA sequences give an alignment, consensus and logo). The structure page shows the JSmol viewer (grey ribbon, each touching residue in its own colour, magenta GAG, labels), the motifs found, the pattern-versus-contact statistics for the functional unit, the sequence, the contact table and the downloads.

## PLIP cross-check

`--plip` writes `<ID>_plip_check.txt`: the protein residues PLIP links to the GAG compared with the residues found by the distance cutoff.

## Test cases

`python run_test_cases.py` (needs internet) runs 1AZX, 1E03, 1XMN, 1FQ9, 1J8I, 2NWG, 1TB6, 1BFB and 2AXM, scores the contacts against published binding residues, counts every motif and the PLIP overlap, and writes `test_case_output/test_case_results.csv`.
