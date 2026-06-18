import re, sys
from pathlib import Path
import pymupdf
from astropy.table import Table

PDF_PATH = Path("/home/haitong/work/galactic_dynamics/data/Cappellari2006_SAURON_MGE.pdf")
OUT_PATH = Path("/home/haitong/work/galactic_dynamics/data/processed/sauron_mge.ecsv")
EM_DASH = "\u2014"

def _parse_table_text(text):
    result = {}
    lines = [l.strip() for l in text.splitlines()]
    i = 0
    while i < len(lines):
        if re.match(r"NGC\s+\d+", lines[i]):
            group = []
            j = i
            while j < len(lines) and re.match(r"NGC\s+\d+", lines[j]) and len(group) < 4:
                clean = "NGC" + re.search(r"(\d+)", lines[j]).group(1)
                group.append(clean)
                j += 1
            if len(group) >= 2:
                components = {g: [] for g in group}
                k = j
                while k < len(lines):
                    line = lines[k]
                    if re.match(r"NGC\s+\d+", line):
                        break
                    if line.startswith("(x") or line.startswith("(L") or "RAS" in line:
                        break
                    try:
                        idx = int(line)
                        values = []
                        for v in range(12):
                            if k + 1 + v < len(lines):
                                values.append(lines[k + 1 + v])
                        if len(values) == 12:
                            for gidx, gname in enumerate(group):
                                base = gidx * 3
                                logI = values[base]
                                log_s = values[base + 1]
                                q = values[base + 2]
                                if EM_DASH in logI or EM_DASH in log_s or EM_DASH in q:
                                    continue
                                try:
                                    components[gname].append((
                                        10.0 ** float(logI),
                                        10.0 ** float(log_s),
                                        float(q),
                                    ))
                                except ValueError:
                                    pass
                        k += 13
                    except ValueError:
                        k += 1
                for gname in group:
                    if components[gname]:
                        result[gname] = components[gname]
            i = j
        else:
            i += 1
    return result

def main():
    doc = pymupdf.open(PDF_PATH)
    all_mge = {}
    for pg in [23, 24]:
        text = doc[pg].get_text()
        page_data = _parse_table_text(text)
        all_mge.update(page_data)
    doc.close()
    print(f"Extracted MGE data for {len(all_mge)} galaxies")
    rows = []
    for galaxy in sorted(all_mge.keys()):
        for j, (I_val, sigma_val, q_val) in enumerate(all_mge[galaxy], start=1):
            rows.append({"galaxy": galaxy, "j": j, "I": I_val,
                         "sigma": sigma_val, "q": q_val, "pa_twist": 0.0})
    tbl = Table(rows=rows, names=["galaxy","j","I","sigma","q","pa_twist"])
    tbl["I"].unit = "L_sun,I/pc^2"
    tbl["sigma"].unit = "arcsec"
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    tbl.write(OUT_PATH, format="ascii.ecsv", overwrite=True)
    print(f"Written {len(rows)} components to {OUT_PATH}")
    for galaxy in sorted(all_mge.keys()):
        print(f"  {galaxy}: {len(all_mge[galaxy])} components")

if __name__ == "__main__":
    main()
