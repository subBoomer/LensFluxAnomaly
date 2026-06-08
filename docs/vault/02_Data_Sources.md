# Data Sources

## Compiled Systems (15 total)

### Radio Systems (7)
| Name | Survey | Band | z_l | z_s | Reference |
|------|--------|------|-----|-----|-----------|
| MG0414+0534 | MG | 8.4 GHz | 0.96 | 2.64 | Katz+1997, Hewitt+1992 |
| B0128+437 | CLASS | 5 GHz | 0.74 | 3.13 | Phillips+2000 |
| B0712+472 | CLASS | 5 GHz | 0.41 | 1.34 | King+1998 |
| B1422+231 | CLASS | 8.4 GHz | 0.34 | 3.62 | Patnaik+1992 |
| B1608+656 | CLASS | 8.4 GHz | 0.63 | 1.39 | Myers+1995 |
| B1933+503 | CLASS | 8.4 GHz | 0.76 | 2.62 | Sykes+1998 |
| B2045+265 | CLASS | 8.5 GHz | 0.87 | 1.28 | Fassnacht+1999 |

### Optical Systems (8, CASTLES F814W)
| Name | Band | Reference |
|------|------|-----------|
| HE0230-2130 | F814W | Pooley+2007 |
| MG0414+0534 | F814W | Pooley+2007 |
| RXJ0911+0551 | F814W | Pooley+2007 |
| SDSSJ0924+0219 | F814W | Pooley+2007 |
| PG1115+080 | F814W | Pooley+2007 |
| RXJ1131-1231 | F814W | Pooley+2007 |
| H1413+117 | F814W | Pooley+2007 |
| B1422+231 | F814W | Pooley+2007 |
| WFI2033-4723 | F814W | Pooley+2007 |
| Q2237+0305 | F814W | Pooley+2007 |

## Data Format
Each system record contains: image positions (x, y in arcsec), fluxes (mJy for radio, 10^(-0.4*mag) for optical), band, survey, and where available lens redshifts and parities.

## Data Files
- `data/radio_quads.py` — curated radio quad catalog
- `data/curated_quads.py` — optical quad catalog from Pooley+2007 / CASTLES
