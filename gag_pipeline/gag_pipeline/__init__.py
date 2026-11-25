"""
gag_pipeline
============
Automated GAG-protein interaction extraction and sequence analysis.

Modules:
    pdb_downloader          - fetch structures from RCSB
    pdb_parser               - parse ATOM/HETATM records
    gag_detector              - identify GAG/carbohydrate residues
    bio_component_selector   - restrict analysis to the relevant biological unit
    neighbor_search           - atom-atom distance / contact search
    residue_contacts          - collapse atom contacts to residue-level contacts
    sequence_extractor         - build interacting-residue sequences + consistency checks
    pdb_exporter               - write the GAG + interacting-residue PDB subset
    report_generator           - CSV / text / sequence / JSON reports
"""

__version__ = "1.0.0"
