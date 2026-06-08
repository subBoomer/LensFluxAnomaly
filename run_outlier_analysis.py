import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.resolve()))

from src.catalog_utils import build_unified_catalog, mag_to_flux
from src.compute_rmin import compute_rmin
from src.statistic import compute_rfold
from data.radio_quads import RADIO_QUAD_CATALOG
from data.curated_quads import QUAD_CATALOG

TAIL_THRESHOLD = 0.2


def get_optical_system(name):
    for s in QUAD_CATALOG:
        if s["name"].replace("J", "") in [name.replace("J", ""), name]:
            return s
    return None


def compute_rmin_from_optical(s):
    positions = np.array([[img["x"], img["y"]] for img in s["images"]], dtype=float)
    mags = np.array([img["mag"] for img in s["images"]], dtype=float)
    fluxes = mag_to_flux(mags)
    return compute_rmin(positions[:, 0], positions[:, 1], fluxes)


def main():
    systems = build_unified_catalog(include_radio=True, include_optical=True, deduplicate=True)

    print("=" * 95)
    print(f"{'System':20s} {'Type':8s} {'Band':20s} {'R_min':8s} {'Tail?':6s} {'R_fold':8s} {'z_s':6s} {'z_l':6s}")
    print("-" * 95)

    table_rows = []
    all_rmin = []
    radio_rmin = []
    optical_rmin = []

    for s in systems:
        name = s["name"]
        stype = "radio" if s["source"] == "radio_quads" else "optical"
        band = s["band"]
        z_s = s.get("z_s", None)
        z_l = s.get("z_l", None)

        rmin = compute_rmin(s["x_arcsec"], s["y_arcsec"], s["fluxes"])
        rmin_str = f"{rmin:.4f}" if rmin is not None else "None"
        tail = rmin is not None and rmin > TAIL_THRESHOLD

        rfold = -1.0
        if stype == "radio":
            for rs in RADIO_QUAD_CATALOG:
                if rs["name"] == name and not rs.get("placeholder", False):
                    mu_arr = np.array(rs["parity"], dtype=float)
                    f_arr = np.array(rs["fluxes_mjy"], dtype=float)
                    x_arr = np.array(rs["x_arcsec"], dtype=float)
                    y_arr = np.array(rs["y_arcsec"], dtype=float)
                    rfold = compute_rfold(x_arr, y_arr, mu_arr, f_arr)
                    break
        rfold_str = f"{rfold:.4f}" if rfold >= 0 else "   -  "

        z_s_str = f"{z_s:.3f}" if z_s else "  -  "
        z_l_str = f"{z_l:.3f}" if z_l else "  -  "

        print(f"{name:20s} {stype:8s} {band:20s} {rmin_str:>8s} {'*' if tail else ' ':6s} {rfold_str:>8s} {z_s_str:>6s} {z_l_str:>6s}")

        if rmin is not None:
            all_rmin.append(rmin)
            table_rows.append({"name": name, "type": stype, "band": band, "rmin": rmin, "tail": tail})
            if stype == "radio":
                radio_rmin.append(rmin)
            else:
                optical_rmin.append(rmin)

    print("-" * 95)
    print(f"\nSummary:")
    print(f"  Total systems:   {len(systems)}")
    print(f"  Valid R_min:     {len(all_rmin)}")
    print(f"  Tail (R>{TAIL_THRESHOLD}):     {sum(1 for r in table_rows if r['tail'])}/{len(all_rmin)} = {sum(1 for r in table_rows if r['tail'])/len(all_rmin)*100:.0f}%")
    print(f"  Radio R_min:     mean={np.mean(radio_rmin):.4f} std={np.std(radio_rmin):.4f} (N={len(radio_rmin)})")
    print(f"  Optical R_min:   mean={np.mean(optical_rmin):.4f} std={np.std(optical_rmin):.4f} (N={len(optical_rmin)})")

    tail_systems = [r for r in table_rows if r['tail']]
    if tail_systems:
        print(f"\n  Systems with R_min > {TAIL_THRESHOLD}:")
        for r in sorted(tail_systems, key=lambda x: -x['rmin']):
            print(f"    {r['name']:20s}  R_min={r['rmin']:.4f}  ({r['type']}, {r['band']})")

    print("\n" + "=" * 95)
    print("Cross-overlap: radio vs optical R_min for MG0414+0534 and B1422+231")
    print("-" * 95)
    for name in ["MG0414+0534", "B1422+231"]:
        radio_r, optical_r = None, None
        for s in systems:
            if s["name"] == name and s["source"] == "radio_quads":
                radio_r = compute_rmin(s["x_arcsec"], s["y_arcsec"], s["fluxes"])
        opt = get_optical_system(name)
        if opt:
            optical_r = compute_rmin_from_optical(opt)
        ratio = f"{radio_r / optical_r:.2f}" if (radio_r and optical_r and optical_r > 0) else "N/A"
        ref = opt.get("ref", "") if opt else ""
        print(f"  {name:20s}  radio R_min={radio_r:.4f}  optical R_min={optical_r:.4f}  ratio={ratio}  ({ref})")

    print("\n" + "=" * 95)
    print("Tail systems ranked by R_min (flag = potential microlensing candidates)")
    print("-" * 95)
    for r in sorted(tail_systems, key=lambda x: -x['rmin']):
        note = ""
        rfold_str = ""
        if r['type'] == 'radio':
            for rs in RADIO_QUAD_CATALOG:
                if rs["name"] == r['name']:
                    mu_arr = np.array(rs["parity"], dtype=float)
                    f_arr = np.array(rs["fluxes_mjy"], dtype=float)
                    x_arr = np.array(rs["x_arcsec"], dtype=float)
                    y_arr = np.array(rs["y_arcsec"], dtype=float)
                    rf = compute_rfold(x_arr, y_arr, mu_arr, f_arr)
                    rfold_str = f"  R_fold={rf:.4f}"
                    if rs.get("band", "").startswith("optical"):
                        note = " [microlensing possible]"
                    break
        elif r['band'] == 'F814W':
            note = " [microlensing possible]"
        print(f"    {r['name']:20s}  R_min={r['rmin']:.4f}{rfold_str}{note}")

    print(f"\n  Cross-overlap insight:")
    print(f"    MG0414: radio R_min vs optical R_min - tests microlensing vs intrinsic substructure")
    print(f"    B1422:  radio R_min vs optical R_min - same test")
    print(f"    If radio < optical for both -> microlensing is contaminating optical sample")
    print(f"    If radio ~ optical -> signal is intrinsic (substructure)")


if __name__ == '__main__':
    main()
