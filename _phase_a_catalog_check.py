import sys; sys.path.insert(0, '.')
from src.catalog_utils import build_unified_catalog
systems = build_unified_catalog(include_radio=True, include_optical=True, deduplicate=True)
print(f'Current catalog: {len(systems)} systems')
print(f'Radio: {sum(1 for s in systems if s["source"]=="radio_quads")}')
print(f'Optical: {sum(1 for s in systems if s["source"]=="curated_quads")}')
print()
print('Names:')
for s in systems:
    print(f'  {s["name"]:20s}  band={s["band"]:20s}  source={s["source"]}')
