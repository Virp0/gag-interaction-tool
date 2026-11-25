import os

OUT_PATH = os.path.join(os.path.dirname(__file__), "TEST1.pdb")

ATOM_FMT = (
    "{rec:<6}{serial:>5} {name:<4}{alt:<1}{resn:>3} {chain:<1}{resseq:>4}{icode:<1}   "
    "{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{temp:>6.2f}          {elem:>2}"
)

lines = []
lines.append("HEADER    TEST STRUCTURE FOR GAG PIPELINE          01-JAN-26   TEST")
lines.append("TITLE     SYNTHETIC TEST STRUCTURE")
lines.append("COMPND    MOL_ID: 1; MOLECULE: TEST PROTEIN; CHAIN: A, C;")

serial = 1


seq_a = ["ALA", "ARG", "LYS", "GLY", "SER", "ARG", "LYS", "THR", "TYR", "GLN"]
for i, resn in enumerate(seq_a, start=1):
    if i in (5, 6):
        base_x = 0.0 + (i - 5) * 3.5
        base_y = 0.0
        base_z = 0.0
    else:
        base_x = 50.0 + i * 3.8
        base_y = 50.0
        base_z = 50.0
    for atom_name, elem, off in [("N", "N", (0.0, 0.0, 0.0)),
                                  ("CA", "C", (1.0, 0.0, 0.0)),
                                  ("C", "C", (2.0, 0.5, 0.0)),
                                  ("O", "O", (2.0, 1.7, 0.0))]:
        lines.append(ATOM_FMT.format(
            rec="ATOM", serial=serial, name=atom_name, alt="", resn=resn,
            chain="A", resseq=i, icode="",
            x=base_x + off[0], y=base_y + off[1], z=base_z + off[2],
            occ=1.00, temp=20.00, elem=elem,
        ))
        serial += 1
lines.append(f"TER   {serial:>5}      {seq_a[-1]:>3} A{len(seq_a):>4}")
serial += 1


for i, resn in enumerate(["GLY", "ALA", "VAL"], start=1):
    base_x, base_y, base_z = 200.0 + i * 3.8, 200.0, 200.0
    for atom_name, elem, off in [("N", "N", (0.0, 0.0, 0.0)),
                                  ("CA", "C", (1.0, 0.0, 0.0)),
                                  ("C", "C", (2.0, 0.5, 0.0)),
                                  ("O", "O", (2.0, 1.7, 0.0))]:
        lines.append(ATOM_FMT.format(
            rec="ATOM", serial=serial, name=atom_name, alt="", resn=resn,
            chain="C", resseq=i, icode="",
            x=base_x + off[0], y=base_y + off[1], z=base_z + off[2],
            occ=1.00, temp=20.00, elem=elem,
        ))
        serial += 1
lines.append(f"TER   {serial:>5}      VAL C   3")
serial += 1


ids_atoms = [
    ("C1", "C", (0.5, 0.3, 0.0)),
    ("C2", "C", (1.5, 0.3, 0.0)),
    ("O2", "O", (2.0, 1.0, 0.5)),
    ("C3", "C", (0.5, 1.5, 0.0)),
]
for name, elem, off in ids_atoms:
    lines.append(ATOM_FMT.format(
        rec="HETATM", serial=serial, name=name, alt="", resn="IDS",
        chain="B", resseq=1, icode="",
        x=off[0], y=off[1], z=off[2],
        occ=1.00, temp=25.00, elem=elem,
    ))
    serial += 1

sgn_atoms = [
    ("C1", "C", (4.2, 0.3, 0.0)),
    ("C2", "C", (5.0, 0.4, 0.0)),
    ("N2", "N", (4.0, -0.2, 0.3)),
    ("S1", "S", (5.5, 1.0, 0.0)),
]
for name, elem, off in sgn_atoms:
    lines.append(ATOM_FMT.format(
        rec="HETATM", serial=serial, name=name, alt="", resn="SGN",
        chain="B", resseq=2, icode="",
        x=off[0], y=off[1], z=off[2],
        occ=1.00, temp=25.00, elem=elem,
    ))
    serial += 1


lines.append(ATOM_FMT.format(
    rec="HETATM", serial=serial, name="O", alt="", resn="HOH",
    chain="B", resseq=101, icode="",
    x=10.0, y=10.0, z=10.0, occ=1.00, temp=30.00, elem="O",
))
serial += 1


lines.append(ATOM_FMT.format(
    rec="HETATM", serial=serial, name="ZN", alt="", resn="ZN ",
    chain="B", resseq=102, icode="",
    x=90.0, y=90.0, z=90.0, occ=1.00, temp=30.00, elem="ZN",
))
serial += 1

lines.append("END")

with open(OUT_PATH, "w") as fh:
    fh.write("\n".join(lines) + "\n")

print(f"Wrote {OUT_PATH}")
