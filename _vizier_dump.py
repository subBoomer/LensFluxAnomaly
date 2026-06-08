from astroquery.vizier import Vizier
import numpy as np

Vizier.ROW_LIMIT = -1

cats = Vizier.get_catalogs('J/AJ/143/119')

t3 = cats[1]  # table3 - lens systems
print('=== Table 3: Lens Systems (54 rows) ===')
print(f'Columns: {t3.colnames}')
t3.pprint(max_lines=60, max_width=200)

print('\n\n=== Table 5: Confirmed Lenses (26 rows) ===')
t5 = cats[3]
print(f'Columns: {t5.colnames}')
t5.pprint(max_lines=30, max_width=200)

print('\n\n=== Table 6: Extended candidates (36 rows) ===')
t6 = cats[4]
print(f'Columns: {t6.colnames}')
t6.pprint(max_lines=40, max_width=200)
