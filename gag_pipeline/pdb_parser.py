from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Dict, Optional

WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD"}
MODIFIED_PROTEIN_RESIDUES = {"MSE", "SEP", "TPO", "PTR", "HYP", "CSO", "CME", "CSD", "OCS", "MLY", "M3L", "ALY", "KCX", "PCA", "SEC", "PYL"}
NUCLEIC_RESNAMES = {"A", "C", "G", "U", "T", "I", "DA", "DC", "DG", "DT", "DI", "DU"}


@dataclass
class Atom:
    serial: int
    name: str
    alt_loc: str
    res_name: str
    chain_id: str
    res_seq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    temp_factor: float
    element: str
    record_type: str
    raw_line: str = field(repr=False, default="")

    @property
    def coord(self):
        return (self.x, self.y, self.z)

    @property
    def res_key(self):
        return (self.chain_id, self.res_seq, self.icode, self.res_name)


@dataclass
class Link:
    name1: str
    res_name1: str
    chain1: str
    res_seq1: int
    icode1: str
    name2: str
    res_name2: str
    chain2: str
    res_seq2: int
    icode2: str

    @property
    def res_key1(self):
        return (self.chain1, self.res_seq1, self.icode1, self.res_name1)

    @property
    def res_key2(self):
        return (self.chain2, self.res_seq2, self.icode2, self.res_name2)


@dataclass
class DbRef:
    chain_id: str
    seq_begin: int
    seq_end: int
    database: str
    accession: str
    db_begin: int


@dataclass
class ParsedStructure:
    pdb_id: str
    header_title: str
    protein_name: str
    atoms: List[Atom]
    protein_atoms: List[Atom]
    hetatm_atoms: List[Atom]
    water_atoms: List[Atom]
    ligand_atoms: List[Atom]
    gag_atoms: List[Atom]
    chains_present: List[str]
    links: List[Link]
    dbrefs: List[DbRef] = field(default_factory=list)
    n_models: int = 1

    ss_ranges: List[tuple] = field(default_factory=list)


def _parse_atom_line(line: str) -> Optional[Atom]:
    record_type = line[0:6].strip()
    if record_type not in ("ATOM", "HETATM"):
        return None
    try:
        serial = int(line[6:11])
        name = line[12:16].strip()
        alt_loc = line[16:17].strip()
        res_name = line[17:20].strip()
        chain_id = line[21:22].strip() or " "
        res_seq = int(line[22:26])
        icode = line[26:27].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occ_str = line[54:60].strip()
        occupancy = float(occ_str) if occ_str else 1.0
        tf_str = line[60:66].strip()
        temp_factor = float(tf_str) if tf_str else 0.0
        element = line[76:78].strip()
        if not element:

            element = "".join(c for c in name if c.isalpha())[:2].strip()
    except (ValueError, IndexError):
        return None

    return Atom(
        serial=serial,
        name=name,
        alt_loc=alt_loc,
        res_name=res_name,
        chain_id=chain_id,
        res_seq=res_seq,
        icode=icode,
        x=x, y=y, z=z,
        occupancy=occupancy,
        temp_factor=temp_factor,
        element=element,
        record_type=record_type,
        raw_line=line.rstrip("\n"),
    )


def _parse_link_line(line: str) -> Optional[Link]:
    try:
        name1 = line[12:16].strip()
        res_name1 = line[17:20].strip()
        chain1 = line[21:22].strip() or " "
        res_seq1 = int(line[22:26])
        icode1 = line[26:27].strip()

        name2 = line[42:46].strip()
        res_name2 = line[47:50].strip()
        chain2 = line[51:52].strip() or " "
        res_seq2 = int(line[52:56])
        icode2 = line[56:57].strip()
    except (ValueError, IndexError):
        return None

    return Link(
        name1=name1, res_name1=res_name1, chain1=chain1, res_seq1=res_seq1, icode1=icode1,
        name2=name2, res_name2=res_name2, chain2=chain2, res_seq2=res_seq2, icode2=icode2,
    )


def _parse_dbref_line(line: str) -> Optional[DbRef]:
    try:
        return DbRef(
            chain_id=line[12:13].strip() or " ",
            seq_begin=int(line[14:18]),
            seq_end=int(line[20:24]),
            database=line[26:32].strip(),
            accession=line[33:41].strip(),
            db_begin=int(line[55:60]),
        )
    except (ValueError, IndexError):
        return None


def map_to_uniprot(dbrefs: List[DbRef], chain_id: str, res_seq: int):
    for d in dbrefs:
        if d.database == "UNP" and d.chain_id == chain_id and d.seq_begin <= res_seq <= d.seq_end:
            return d.accession, res_seq - d.seq_begin + d.db_begin
    return None


def parse_pdb(filepath: str) -> ParsedStructure:
    pdb_id = ""
    header_title = ""
    compnd_lines = []

    atoms: List[Atom] = []
    links: List[Link] = []
    dbrefs: List[DbRef] = []
    ss_ranges: List[tuple] = []
    n_models = 0
    in_extra_model = False

    with open(filepath, "r", errors="replace") as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec == "MODEL":
                n_models += 1
                continue
            if rec == "ENDMDL":
                in_extra_model = True
                continue
            if rec == "HEADER" and len(line) > 62:
                pdb_id = line[62:66].strip()
            elif rec == "TITLE":
                header_title += line[10:80].strip() + " "
            elif rec == "COMPND":
                compnd_lines.append(line[10:80].strip())
            elif rec in ("ATOM", "HETATM"):

                if in_extra_model:
                    continue
                atom = _parse_atom_line(line)
                if atom is not None:
                    atoms.append(atom)
            elif rec == "DBREF":
                dbref = _parse_dbref_line(line)
                if dbref is not None:
                    dbrefs.append(dbref)
            elif rec in ("HELIX", "SHEET"):
                try:
                    if rec == "HELIX":
                        ss_ranges.append((line[19].strip(), int(line[21:25]), int(line[33:37])))
                    else:
                        ss_ranges.append((line[21].strip(), int(line[22:26]), int(line[33:37])))
                except (ValueError, IndexError):
                    pass
            elif rec == "LINK":
                link = _parse_link_line(line)
                if link is not None:
                    links.append(link)

    protein_name = ""
    for cl in compnd_lines:
        if "MOLECULE:" in cl:
            protein_name = cl.split("MOLECULE:")[1].split(";")[0].strip()
            break

    protein_atoms = [a for a in atoms
                     if a.res_name not in NUCLEIC_RESNAMES
                     and (a.record_type == "ATOM" or a.res_name in MODIFIED_PROTEIN_RESIDUES)]
    hetatm_atoms = [a for a in atoms if a.record_type == "HETATM" and a.res_name not in MODIFIED_PROTEIN_RESIDUES]
    water_atoms = [a for a in hetatm_atoms if a.res_name in WATER_RESNAMES]

    ligand_atoms = [a for a in hetatm_atoms if a.res_name not in WATER_RESNAMES]

    chains_present = sorted({a.chain_id for a in atoms})

    return ParsedStructure(
        pdb_id=pdb_id,
        header_title=header_title.strip(),
        protein_name=protein_name,
        atoms=atoms,
        protein_atoms=protein_atoms,
        hetatm_atoms=hetatm_atoms,
        water_atoms=water_atoms,
        ligand_atoms=ligand_atoms,
        gag_atoms=[],
        chains_present=chains_present,
        links=links,
        dbrefs=dbrefs,
        n_models=max(n_models, 1),
        ss_ranges=ss_ranges,
    )


def group_atoms_by_residue(atoms: List[Atom]) -> Dict[tuple, List[Atom]]:
    groups: Dict[tuple, List[Atom]] = {}
    for a in atoms:
        groups.setdefault(a.res_key, []).append(a)
    return groups
