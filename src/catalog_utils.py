import numpy as np
from data.radio_quads import RADIO_QUAD_CATALOG
from data.curated_quads import QUAD_CATALOG


def mag_to_flux(mag):
    return 10.0 ** (-0.4 * mag)


def build_unified_catalog(include_radio=True, include_optical=True, deduplicate=True):
    systems = []
    seen_names = set()

    if include_radio:
        for s in RADIO_QUAD_CATALOG:
            if s.get("placeholder", False):
                continue
            systems.append({
                "name": s["name"],
                "survey": s["survey"],
                "band": s["band"],
                "z_l": s.get("z_l"),
                "z_s": s.get("z_s"),
                "x_arcsec": np.array(s["x_arcsec"], dtype=float),
                "y_arcsec": np.array(s["y_arcsec"], dtype=float),
                "fluxes": np.array(s["fluxes_mjy"], dtype=float),
                "flux_unit": "mJy",
                "flux_err": None,
                "has_parity": True,
                "parity": s.get("parity"),
                "source": "radio_quads",
                "n_images": 4,
            })
            seen_names.add(s["name"])

    if include_optical:
        for s in QUAD_CATALOG:
            name = s["name"]
            if deduplicate and name in seen_names:
                print(f'  [skip] {name} already in radio catalog')
                continue
            positions = np.array([[img["x"], img["y"]] for img in s["images"]], dtype=float)
            mags = np.array([img["mag"] for img in s["images"]], dtype=float)
            fluxes = mag_to_flux(mags)
            systems.append({
                "name": name,
                "survey": s["survey"],
                "band": s["band"],
                "z_l": None,
                "z_s": None,
                "x_arcsec": positions[:, 0],
                "y_arcsec": positions[:, 1],
                "fluxes": fluxes,
                "flux_unit": "mag_flux",
                "flux_err": fluxes * 0.05,
                "has_parity": False,
                "parity": None,
                "source": "curated_quads",
                "n_images": s["n_images"],
            })
            seen_names.add(name)

    return systems


def describe_catalog(systems):
    radio = sum(1 for s in systems if s["source"] == "radio_quads")
    optical = sum(1 for s in systems if s["source"] == "curated_quads")
    print(f"  Total: {len(systems)} systems ({radio} radio, {optical} optical)")
    bands = set(s["band"] for s in systems)
    print(f"  Bands: {', '.join(sorted(bands))}")
    sys_per_source = {}
    for s in systems:
        sys_per_source.setdefault(s["source"], []).append(s["name"])
    for source, names in sorted(sys_per_source.items()):
        print(f"  {source}: {', '.join(names)}")
