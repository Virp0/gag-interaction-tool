import os

OUT_PATH = os.path.join(os.path.dirname(__file__), "REAL_1AZX_EXCERPT.pdb")

ATOM_FMT = (
    "{rec:<6}{serial:>5} {name:<4}{alt:<1}{resn:>3} {chain:<1}{resseq:>4}{icode:<1}   "
    "{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{temp:>6.2f}          {elem:>2}"
)

header_lines = [
    "HEADER    SERINE PROTEASE INHIBITOR               23-NOV-97   TEST              ",
    "TITLE     TRIMMED REAL EXCERPT OF 1AZX FOR REGRESSION TESTING",
]

het_block = """\
HET    Z9L  A   1      25
HET    Z9K  A   2      17
HET    GU6  A   3      23
HET    GU1  A   4      14
HET    Z9H  A   5      21
HET    NAG  I 501      14
HETNAM     Z9L METHYL 2,3,6-TRI-O-SULFO-ALPHA-D-GLUCOPYRANOSIDE
HETNAM     Z9K 3-O-METHYL-2-O-SULFO-ALPHA-L-IDOPYRANURONIC ACID
HETNAM     GU6 2,3,6-TRI-O-SULFO-ALPHA-D-GLUCOPYRANOSE
HETNAM     GU1 2,3-DI-O-METHYL-BETA-D-GLUCOPYRANURONIC ACID
HETNAM     Z9H 3,4-DI-O-METHYL-2,6-DI-O-SULFO-ALPHA-D-GLUCOPYRANOSE
HETNAM     NAG 2-ACETAMIDO-2-DEOXY-BETA-D-GLUCOPYRANOSE
"""

link_block = """\
LINK         ND2 ASN I 155                 C1  NAG I 501     1555   1555  1.43
LINK         O4  Z9L A   1                 C1  Z9K A   2     1555   1555  1.42
LINK         O4  Z9K A   2                 C1  GU6 A   3     1555   1555  1.44
LINK         O4  GU6 A   3                 C1  GU1 A   4     1555   1555  1.42
LINK         O4  GU1 A   4                 C1  Z9H A   5     1555   1555  1.43
"""

body_lines = []
serial = 1


protein_residues = [
    (11, "LYS", (0.0, 0.0, 0.0)),
    (13, "ARG", (3.0, 0.0, 0.0)),
    (46, "ARG", (6.0, 0.0, 0.0)),
    (47, "ARG", (9.0, 0.0, 0.0)),
    (155, "ASN", (100.0, 100.0, 100.0)),
]
for resseq, resn, (bx, by, bz) in protein_residues:
    for atom_name, elem, off in [("N", "N", (0.0, 0.0, 0.0)),
                                  ("CA", "C", (1.0, 0.0, 0.0)),
                                  ("C", "C", (2.0, 0.5, 0.0)),
                                  ("O", "O", (2.0, 1.7, 0.0)),
                                  ("CB", "C", (1.0, -1.2, 0.0))]:
        body_lines.append(ATOM_FMT.format(
            rec="ATOM", serial=serial, name=atom_name, alt="", resn=resn,
            chain="I", resseq=resseq, icode="",
            x=bx + off[0], y=by + off[1], z=bz + off[2],
            occ=1.00, temp=20.00, elem=elem,
        ))
        serial += 1
    if resn == "ASN":
        for atom_name, elem, off in [("CG", "C", (2.0, -1.5, 0.3)),
                                      ("OD1", "O", (2.3, -1.2, 1.5)),
                                      ("ND2", "N", (2.2, -2.8, 0.0))]:
            body_lines.append(ATOM_FMT.format(
                rec="ATOM", serial=serial, name=atom_name, alt="", resn=resn,
                chain="I", resseq=resseq, icode="",
                x=bx + off[0], y=by + off[1], z=bz + off[2],
                occ=1.00, temp=20.00, elem=elem,
            ))
            serial += 1
body_lines.append(f"TER   {serial:>5}      ASN I 155")
serial += 1


het_residues = [
    (1, "Z9L", [("C1", "C", (0.5, 0.3, 0.0)), ("O1S", "O", (0.8, 1.0, 0.5))]),
    (2, "Z9K", [("C1", "C", (3.5, 0.3, 0.0)), ("O1S", "O", (3.8, 1.0, 0.5))]),
    (3, "GU6", [("C1", "C", (6.5, 0.3, 0.0)), ("O1S", "O", (6.8, 1.0, 0.5))]),
    (4, "GU1", [("C1", "C", (9.5, 0.3, 0.0)), ("O1", "O", (9.8, 1.0, 0.5))]),
    (5, "Z9H", [("C1", "C", (12.5, 0.3, 0.0)), ("O1S", "O", (12.8, 1.0, 0.5))]),
]
for resseq, resn, atoms in het_residues:
    for name, elem, off in atoms:
        body_lines.append(ATOM_FMT.format(
            rec="HETATM", serial=serial, name=name, alt="", resn=resn,
            chain="A", resseq=resseq, icode="",
            x=off[0], y=off[1], z=off[2],
            occ=1.00, temp=25.00, elem=elem,
        ))
        serial += 1


body_lines.append(ATOM_FMT.format(
    rec="HETATM", serial=serial, name="C1", alt="", resn="NAG",
    chain="I", resseq=501, icode="",
    x=102.2, y=97.2, z=100.0, occ=1.00, temp=30.00, elem="C",
))
serial += 1
body_lines.append(ATOM_FMT.format(
    rec="HETATM", serial=serial, name="C2", alt="", resn="NAG",
    chain="I", resseq=501, icode="",
    x=103.0, y=97.5, z=100.3, occ=1.00, temp=30.00, elem="C",
))
serial += 1
body_lines.append("END")

full_text = "\n".join(header_lines) + "\n" + het_block + link_block + "\n".join(body_lines) + "\n"

with open(OUT_PATH, "w") as fh:
    fh.write(full_text)

print(f"Wrote {OUT_PATH}")
