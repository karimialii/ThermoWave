"""Generate the four "LM6000-class" utility-scale component maps by scaling
the packaged T100 map shapes to a twin-spool aeroderivative's design points.

Real utility-machine component maps (GE LM6000, 9HA, Siemens SGT5, ...) are
OEM-proprietary and not published. The standard industry workaround — the
same one GasTurb/NPSS practitioners use — is *map scaling*: take a map whose
speed-line topology (surge/choke behavior, efficiency islands) is known, and
rescale its corrected speed, corrected mass flow, pressure ratio, and peak
efficiency to the target machine's published design point:

    A' = A * sA            (corrected speed)
    B' = B * sB            (corrected mass flow)
    PR' = 1 + (PR - 1)*sPR (pressure ratio, scaled about PR = 1)
    eff: the file's E_fact conversion factor is replaced so the map's peak
         matches a utility-class polytropic quality

The shapes come from the packaged Turbec T100 maps; the design points are
LM6000-class, from GE's published data (two-shaft aeroderivative, overall
pressure ratio 29:1, LP spool at the CF6-80C2's 3600 rev/min synchronous
speed for direct drive, ~130 kg/s airflow, ~40+ MWe — see
https://www.geaerospace.com/sites/default/files/datasheet-lm6000.pdf and
GE's aeroderivative design & operating features paper). The result is NOT
GE data — it's a self-consistent, realistically-sized set of maps for a
generic utility-scale twin-spool machine, which is exactly what an open
demonstration needs.

Cycle design points used (ISO inlet, gamma = 1.4 ideal-gas air):
    LPC: PR 2.4  at 3600 rev/min,  127 kg/s          (A 3.54,  B 2128)
    HPC: PR 12.2 at 10000 rev/min, T_in ~383 K       (A 8.52,  B 1021)
    HPT: PR ~4.6 at 10000 rev/min, T_in ~1450 K      (A 4.38,  B  173)
    LPT: PR ~6.2 at 3600 rev/min,  T_in ~1000 K      (A 1.90,  B  655)

Run: .venv/bin/python generate_lm6000_class_maps.py
(Regenerates the four "LM6000-class *.cop/.tur" files in place.)
"""

from pathlib import Path

SENTINEL = "-999999999"

# (source file, output file, anchor A0/B0/PR0 read off the source map near
#  its own design point, target A_d/B_d/PR_d, new E_fact)
SPECS = [
    ("T100 Comp.cop", "LM6000-class LPC.cop",
     63.42, 11.6, 3.70, 3.536, 2128.0, 2.4, 1.103),
    ("T100 Comp.cop", "LM6000-class HPC.cop",
     63.42, 11.6, 3.70, 8.517, 1021.0, 12.2, 1.090),
    ("T100 Turb.tur", "LM6000-class HPT.tur",
     32.42, 6.69, 3.38, 4.377, 173.0, 4.56, 1.023),
    ("T100 Turb.tur", "LM6000-class LPT.tur",
     32.42, 6.69, 3.38, 1.897, 655.0, 6.15, 1.035),
]


def scale_map(source, A0, B0, PR0, A_d, B_d, PR_d, e_fact):
    sA, sB, sPR = A_d / A0, B_d / B0, (PR_d - 1.0) / (PR0 - 1.0)
    lines = Path(source).read_text().splitlines()
    out = []
    in_pr_section = False
    in_eff_section = False
    row_in_group = 0  # 0 = speed line, 1 = mass-flow row, 2 = value row

    def scale_row(line, factor, affine=False):
        tokens = line.split()
        scaled = []
        for token in tokens:
            value = float(token)
            if value <= float(SENTINEL):
                scaled.append(token)
                break
            scaled.append(
                f"{1.0 + (value - 1.0) * factor:.6g}" if affine else f"{value * factor:.6g}"
            )
        return "  " + "  ".join(scaled)

    for line in lines:
        stripped = line.strip()
        if "Pressure Ratio vs Non-Dimensional Mass Flow" in line:
            in_pr_section, in_eff_section, row_in_group = True, False, 0
            out.append(line)
            continue
        if "Efficiency Or Corrected Work vs Non-Dimensional Mass Flow" in line:
            in_pr_section, in_eff_section, row_in_group = False, True, 0
            out.append(line)
            continue
        if "(E_fact)" in line:
            out.append(f"       {e_fact}   (E_fact) To convert to  fraction")
            continue
        if not (in_pr_section or in_eff_section):
            out.append(line)
            continue
        # Inside a data section.
        if "Angle" in line or stripped.startswith(SENTINEL) or not stripped:
            row_in_group = 0
            out.append(line)
            continue
        if row_in_group == 0:  # iso-speed value
            out.append(f"       {float(stripped) * sA:.6g}")
            row_in_group = 1
        elif row_in_group == 1:  # corrected mass flow row
            out.append(scale_row(line, sB))
            row_in_group = 2
        else:  # PR row (PR section) or efficiency row (eff section)
            out.append(scale_row(line, sPR, affine=True) if in_pr_section else line)
            row_in_group = 0
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    for source, target, A0, B0, PR0, A_d, B_d, PR_d, e_fact in SPECS:
        Path(target).write_text(
            scale_map(source, A0, B0, PR0, A_d, B_d, PR_d, e_fact)
        )
        print(f"wrote {target}")
