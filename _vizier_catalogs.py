from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = -1
try:
    cats = Vizier.get_catalogs('J/AJ/143/119')
    print(f'SQLS catalog: {len(cats)} tables')
    for i, c in enumerate(cats):
        name = c.meta.get('name', c.meta.get('title', '?'))
        print(f'  Table {i}: {name} ({len(c)} rows)')
        for col in c.colnames[:10]:
            print(f'    {col}: {c[col].dtype}')
except Exception as e:
    print(f'SQLS error: {e}')

try:
    from astroquery.vizier import VizierClass
    v = VizierClass(keywords=["SQLS"], row_limit=-1)
    results = v()
    print(f'\nSQLS keyword search: {len(results)} catalogs')
except Exception as e:
    print(f'SQLS keyword search error: {e}')

try:
    cats2 = Vizier.get_catalogs('J/A+A/688/A34')
    print(f'\nKiDS/HSC catalog: {len(cats2)} tables')
    for i, c in enumerate(cats2):
        name = c.meta.get('name', c.meta.get('title', '?'))
        print(f'  Table {i}: {name} ({len(c)} rows)')
        for col in c.colnames[:10]:
            print(f'    {col}: {c[col].dtype}')
except Exception as e:
    print(f'HSC error: {e}')
