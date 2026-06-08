import os
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
from typing import Dict, List, Optional


DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
CATALOGS_DIR = os.path.join(DATA_DIR, "catalogs")

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


def _ensure_dirs():
    for d in [RAW_DIR, PROCESSED_DIR, CATALOGS_DIR]:
        os.makedirs(d, exist_ok=True)


def _coord_to_arcsec_offset(
    ra_img: float, dec_img: float,
    ra_lens: float, dec_lens: float,
) -> tuple:
    coord_img = SkyCoord(ra_img, dec_img, unit=u.deg)
    coord_lens = SkyCoord(ra_lens, dec_lens, unit=u.deg)
    sep = coord_img.separation(coord_lens).arcsec
    pa = coord_img.position_angle(coord_lens).rad
    dx = sep * np.sin(pa)
    dy = sep * np.cos(pa)
    return dx, dy


def _make_system_id(name: str, survey: str) -> str:
    return f"{survey.upper()}_{name.replace(' ', '_')}"


class IngestPipeline:
    def __init__(self):
        _ensure_dirs()
        self._catalog = []

    def _add_system(self, row: dict):
        self._catalog.append(row)

    def ingest_castles(self) -> int:
        url = "https://lweb.cfa.harvard.edu/castles/"
        try:
            resp = requests.get(url, timeout=60, headers={"User-Agent": "Mozilla/5.0"})
            resp.raise_for_status()
        except requests.RequestException as e:
            print(f"  CASTLES unavailable: {e}")
            return 0
        soup = BeautifulSoup(resp.text, "html.parser")

        systems_found = 0
        tables = soup.find_all("table")
        for table in tables:
            rows = table.find_all("tr")
            for tr in rows:
                cells = tr.find_all("td")
                if len(cells) < 8:
                    continue
                try:
                    name_cell = cells[0].get_text(strip=True)
                    name = name_cell.split(":")[-1].strip() if ":" in name_cell else name_cell
                    if not name:
                        continue
                    n_images_str = cells[6].get_text(strip=True) if len(cells) > 6 else ""
                    n_images = 4 if "4" in n_images_str else 2
                    sep_str = cells[7].get_text(strip=True) if len(cells) > 7 else ""
                    sep = float(sep_str) if sep_str.replace(".", "").isdigit() else None
                except (ValueError, IndexError):
                    continue

                system_id = _make_system_id(name, "CASTLES")
                for img_id in range(1, n_images + 1):
                    self._add_system({
                        "system_id": system_id,
                        "image_id": f"{name}_img{img_id}",
                        "ra_image": np.nan,
                        "dec_image": np.nan,
                        "x_image": np.nan,
                        "y_image": np.nan,
                        "x_lens": 0.0,
                        "y_lens": 0.0,
                        "flux": np.nan,
                        "flux_err": np.nan,
                        "band": "H",
                        "survey": "CASTLES",
                        "n_images": n_images,
                    })
                systems_found += 1

        return systems_found

    def ingest_slacs(self) -> int:
        Vizier.ROW_LIMIT = -1
        try:
            catalogs = Vizier.get_catalogs("J/ApJ/682/964")
        except Exception:
            return 0

        systems_found = 0
        for cat in catalogs:
            table = cat.to_pandas()
            for _, row in table.iterrows():
                try:
                    name = str(row.get("Name", ""))
                    if not name:
                        continue
                    ra_lens = float(row.get("RAJ2000", row.get("_RA", np.nan)))
                    dec_lens = float(row.get("DEJ2000", row.get("_DE", np.nan)))
                    if np.isnan(ra_lens):
                        continue
                except (ValueError, TypeError):
                    continue

                system_id = _make_system_id(name, "SLACS")
                self._add_system({
                    "system_id": system_id,
                    "image_id": f"{name}_img1",
                    "ra_image": ra_lens,
                    "dec_image": dec_lens,
                    "x_image": 0.0,
                    "y_image": 0.0,
                    "x_lens": 0.0,
                    "y_lens": 0.0,
                    "flux": np.nan,
                    "flux_err": np.nan,
                    "band": "F814W",
                    "survey": "SLACS",
                    "n_images": 2,
                })
                systems_found += 1

        return systems_found

    def ingest_vizier_catalog(self, vizier_id: str, survey: str, band: str = "I") -> int:
        Vizier.ROW_LIMIT = -1
        try:
            catalogs = Vizier.get_catalogs(vizier_id)
        except Exception:
            return 0

        systems_found = 0
        for cat in catalogs:
            table = cat.to_pandas()
            for _, row in table.iterrows():
                try:
                    name = str(row.get("Name", row.get("_tab1_1", "")))
                    if not name:
                        continue
                    ra_lens = float(row.get("RAJ2000", row.get("_RA", np.nan)))
                    dec_lens = float(row.get("DEJ2000", row.get("_DE", np.nan)))
                    if np.isnan(ra_lens):
                        continue
                except (ValueError, TypeError):
                    continue

                system_id = _make_system_id(name, survey)
                self._add_system({
                    "system_id": system_id,
                    "image_id": f"{name}_img1",
                    "ra_image": ra_lens,
                    "dec_image": dec_lens,
                    "x_image": 0.0,
                    "y_image": 0.0,
                    "x_lens": 0.0,
                    "y_lens": 0.0,
                    "flux": np.nan,
                    "flux_err": np.nan,
                    "band": band,
                    "survey": survey,
                    "n_images": 2,
                })
                systems_found += 1

        return systems_found

    def ingest_lenscat(self) -> int:
        lenscat_path = None
        possible = [
            os.path.join(os.path.dirname(__file__), "..", ".venv", "Lib", "site-packages", "lenscat", "data", "catalog.csv"),
            os.path.join(os.path.dirname(__file__), "..", ".venv", "Lib", "site-packages", "lenscat", "data", "catalog.fits"),
        ]
        for p in possible:
            if os.path.exists(p):
                lenscat_path = p
                break

        if lenscat_path is None:
            return 0

        try:
            if lenscat_path.endswith(".csv"):
                df = pd.read_csv(lenscat_path)
            else:
                from astropy.table import Table
                t = Table.read(lenscat_path)
                df = t.to_pandas()
        except Exception:
            return 0

        systems_found = 0
        for _, row in df.iterrows():
            try:
                name = str(row.get("name", ""))
                if not name:
                    continue
                ra = float(row.get("RA", np.nan))
                dec = float(row.get("DEC", np.nan))
                if np.isnan(ra):
                    continue
            except (ValueError, TypeError):
                continue

            system_id = _make_system_id(name, "LENSCAT")
            self._add_system({
                "system_id": system_id,
                "image_id": f"{name}_img1",
                "ra_image": ra,
                "dec_image": dec,
                "x_image": 0.0,
                "y_image": 0.0,
                "x_lens": 0.0,
                "y_lens": 0.0,
                "flux": np.nan,
                "flux_err": np.nan,
                "band": "unknown",
                "survey": "LENSCAT",
                "n_images": 2,
            })
            systems_found += 1

        return systems_found

    def ingest_des_lenses(self) -> int:
        Vizier.ROW_LIMIT = -1
        try:
            catalogs = Vizier.get_catalogs("J/ApJ/900/1")
        except Exception:
            return 0

        systems_found = 0
        for cat in catalogs:
            table = cat.to_pandas()
            for _, row in table.iterrows():
                try:
                    name = str(row.get("Name", row.get("_tab1_1", "")))
                    if not name:
                        continue
                    ra_lens = float(row.get("RAJ2000", row.get("_RA", np.nan)))
                    dec_lens = float(row.get("DEJ2000", row.get("_DE", np.nan)))
                    if np.isnan(ra_lens):
                        continue
                except (ValueError, TypeError):
                    continue

                system_id = _make_system_id(name, "DES")
                self._add_system({
                    "system_id": system_id,
                    "image_id": f"{name}_img1",
                    "ra_image": ra_lens,
                    "dec_image": dec_lens,
                    "x_image": 0.0,
                    "y_image": 0.0,
                    "x_lens": 0.0,
                    "y_lens": 0.0,
                    "flux": np.nan,
                    "flux_err": np.nan,
                    "band": "r",
                    "survey": "DES",
                    "n_images": 2,
                })
                systems_found += 1

        return systems_found

    def run_all(self) -> pd.DataFrame:
        results = {}

        print("[ingest] CASTLES...")
        results["CASTLES"] = self.ingest_castles()

        print("[ingest] SLACS (Vizier)...")
        results["SLACS"] = self.ingest_slacs()

        print("[ingest] BELLS (Vizier)...")
        results["BELLS"] = self.ingest_vizier_catalog("J/ApJ/744/41", "BELLS", "F814W")

        print("[ingest] DES lenses (Vizier)...")
        results["DES"] = self.ingest_des_lenses()

        print("[ingest] KiDS strong lenses (Vizier)...")
        results["KiDS"] = self.ingest_vizier_catalog("J/A+A/688/A34", "KiDS", "r")

        print("[ingest] lenscat (cross-reference)...")
        results["LENSCAT"] = self.ingest_lenscat()

        for survey, count in results.items():
            print(f"  {survey}: {count} systems")

        df = pd.DataFrame(self._catalog, columns=COLUMNS)
        return df


def save_catalog(df: pd.DataFrame):
    path = os.path.join(PROCESSED_DIR, "quad_catalog.fits")
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


def load_catalog() -> pd.DataFrame:
    fits_path = os.path.join(PROCESSED_DIR, "quad_catalog.fits")
    csv_path = fits_path.replace(".fits", ".csv")
    if os.path.exists(fits_path):
        from astropy.table import Table
        t = Table.read(fits_path)
        return t.to_pandas()
    elif os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    else:
        return pd.DataFrame()
