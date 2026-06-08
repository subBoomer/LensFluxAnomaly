---
tags:
  - data
  - catalog
  - lensfluxanomaly
---
# Data Sources

## Radio Sample (8 systems)
CLASS/MG/VLA surveys, 5-8.5 GHz, resolved flux ratios.

See `data/radio_quads.py` for full details.

| System | Survey | Frequency | References |
|--------|--------|-----------|------------|
| MG0414+0534 | MG | 8.4 GHz | Katz+1997, Hewitt+1992 |
| B0128+437 | CLASS | 5 GHz | Phillips+2000 |
| B0712+472 | CLASS | 5 GHz | King+1998, Jackson+1998 |
| B1422+231 | CLASS | 8.4 GHz | Patnaik+1992 |
| B1555+375 | CLASS | 5 GHz | Marlow+2001 |
| B1608+656 | CLASS | 8.4 GHz | Myers+1995, King+1998 |
| B1933+503 | CLASS | 8.4 GHz | Sykes+1998 |
| B2045+265 | CLASS | 8.5 GHz | Fassnacht+1999 |

All have published parity from lens models. See `radio_quads.py` `parity_source` field for provenance.

## Optical Sample (13 systems)
CASTLES HST WFPC2 F814W imaging + literature.

See `data/curated_quads.py` for full details.

| System | Survey | Band | References |
|--------|--------|------|------------|
| HE0230-2130 | CASTLES | F814W | Pooley+2007 |
| MG0414+0534 | CASTLES | F814W | Pooley+2007 |
| RXJ0911+0551 | CASTLES | F814W | Pooley+2007 |
| SDSSJ0924+0219 | CASTLES | F814W | Pooley+2007 |
| PG1115+080 | CASTLES | F814W | Pooley+2007 |
| RXJ1131-1231 | CASTLES | F814W | Pooley+2007 |
| H1413+117 | CASTLES | F814W | Pooley+2007 |
| B1422+231 | CASTLES | F814W | Pooley+2007 |
| WFI2033-4723 | CASTLES | F814W | Pooley+2007 |
| Q2237+0305 | CASTLES | F814W | Pooley+2007 |
| HE0435-1223 | CASTLES | F814W | CASTLES |
| SDSSJ1138+0314 | CASTLES | F814W | CASTLES |
| HS0810+2554 | CASTLES | F814W | CASTLES |

Optical entries have **no parity** — only positions and Vega magnitudes.

## Duplicates
MG0414+0534 and B1422+231 appear in both samples. The radio version is primary (has parity). The optical version is used only for microlensing cross-check.

## Sample Completeness
The CLASS radio quad sample is likely **complete** at 5-8 GHz. A thorough search of CLASS/JVAS (Browne+2003, Myers+2003, King+1998) identified 7 radio quads. B1555+375 was activated as the 8th. No other simple radio quads with resolved flux ratios exist.

Previous studies: Jackson+2015 (4 optical quads at VLA), Jackson+2024 (70 lens survey, but excludes quads). The Jackson+2015 data covers HS0810+2554, RXJ0911+0551, HE0435-1223, SDSSJ0924+0219 (total flux only, not resolved per-image).

See [[03_Catalog_Table]] for the consolidated catalog.
