import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")

COLUMNS = [
    "system_id",
    "image_id",
    "ra_image", "dec_image",
    "x_image", "y_image",
    "x_lens", "y_lens",
    "flux", "flux_err",
    "band", "survey",
    "n_images",
]


def mag_to_flux(mag: float) -> float:
    return 10.0 ** (-0.4 * mag)


def build_literature_catalog() -> pd.DataFrame:
    from data.curated_quads import QUAD_CATALOG

    rows = []
    for quad in QUAD_CATALOG:
        system_id = quad["name"]
        lens = quad["lens"]
        band = quad["band"]
        survey = quad["survey"]
        n_images = quad["n_images"]

        for img in quad["images"]:
            flux = mag_to_flux(img["mag"])
            rows.append({
                "system_id": system_id,
                "image_id": img["id"],
                "ra_image": np.nan,
                "dec_image": np.nan,
                "x_image": img["x"],
                "y_image": img["y"],
                "x_lens": lens["x"],
                "y_lens": lens["y"],
                "flux": flux,
                "flux_err": flux * 0.05,
                "band": band,
                "survey": survey,
                "n_images": n_images,
            })

    return pd.DataFrame(rows, columns=COLUMNS)


def load_literature_catalog() -> pd.DataFrame:
    path = os.path.join(PROCESSED_DIR, "literature_quads.fits")
    csv_path = path.replace(".fits", ".csv")
    if os.path.exists(path):
        from astropy.table import Table
        t = Table.read(path)
        return t.to_pandas()
    elif os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return build_literature_catalog()


def save_literature_catalog(df: pd.DataFrame):
    path = os.path.join(PROCESSED_DIR, "literature_quads.fits")
    try:
        from astropy.table import Table
        t = Table.from_pandas(df)
        t.write(path, overwrite=True)
        print(f"[save] FITS: {path}")
    except Exception as e:
        print(f"[save] FITS failed ({e}), saving CSV instead")
        csv_path = path.replace(".fits", ".csv")
        df.to_csv(csv_path, index=False)
        print(f"[save] CSV: {csv_path}")


if __name__ == "__main__":
    df = build_literature_catalog()
    print(f"Built literature catalog: {len(df)} image rows, {df['system_id'].nunique()} systems")
    print(f"Columns: {list(df.columns)}")
    print(df.to_string(index=False))
    for system_id, group in df.groupby("system_id"):
        total_flux = group["flux"].sum()
        r_min = group["x_image"].min()
        print(f"  {system_id}: {len(group)} images, total_flux={total_flux:.4e}")
    save_literature_catalog(df)
