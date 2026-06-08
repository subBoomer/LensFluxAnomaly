import numpy as np
import pandas as pd
from typing import Dict, List, Optional


def is_quad(system: pd.Series) -> bool:
    return int(system.get("n_images", 0)) == 4


def filter_quad_systems(catalog: pd.DataFrame) -> pd.DataFrame:
    mask = catalog["n_images"] == 4
    return catalog[mask].copy()


def has_valid_fluxes(system: pd.DataFrame, flux_col: str = "flux") -> bool:
    fluxes = system[flux_col].values
    return bool(np.all(np.isfinite(fluxes)) and np.all(fluxes > 0))


def has_unique_image_ids(system: pd.DataFrame, id_col: str = "image_id") -> bool:
    ids = system[id_col].values
    return len(ids) == len(set(ids))


def normalize_fluxes(system: pd.DataFrame, flux_col: str = "flux") -> pd.DataFrame:
    out = system.copy()
    total = out[flux_col].sum()
    if total > 0:
        out[flux_col] = out[flux_col] / total
    return out


def filter_catalog(
    catalog: pd.DataFrame,
    quad_only: bool = True,
    discard_missing_fluxes: bool = True,
    discard_ambiguous_ids: bool = True,
    flux_col: str = "flux",
    id_col: str = "image_id",
) -> pd.DataFrame:
    if quad_only:
        catalog = filter_quad_systems(catalog)

    valid_systems = []
    for system_id, group in catalog.groupby("system_id"):
        valid = True
        if discard_missing_fluxes and not has_valid_fluxes(group, flux_col):
            valid = False
        if discard_ambiguous_ids and not has_unique_image_ids(group, id_col):
            valid = False
        if valid:
            valid_systems.append(system_id)

    return catalog[catalog["system_id"].isin(valid_systems)].copy()


def apply_all_filters(
    raw_catalog: pd.DataFrame,
    config: Optional[Dict] = None,
) -> pd.DataFrame:
    cfg = config or {}
    filtered = filter_catalog(
        raw_catalog,
        quad_only=cfg.get("quad_only", True),
        discard_missing_fluxes=cfg.get("discard_missing_fluxes", True),
        discard_ambiguous_ids=cfg.get("discard_ambiguous_ids", True),
    )

    normalized_list = []
    for system_id, group in filtered.groupby("system_id"):
        norm = normalize_fluxes(group)
        norm["system_id"] = system_id
        normalized_list.append(norm)

    return pd.concat(normalized_list, ignore_index=True)
