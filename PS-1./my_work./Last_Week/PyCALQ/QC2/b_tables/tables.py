#!/usr/bin/env python3
"""
Combined tables module for box matrix elements
"""

from b1 import TABLE_B1
from b2 import TABLE_B2
from b3 import TABLE_B3
from b4 import TABLE_B4
from b5 import TABLE_B5
from b6 import TABLE_B6
from b7 import TABLE_B7
from b8 import TABLE_B8
ALL_TABLES = {}
ALL_TABLES.update(TABLE_B1)
ALL_TABLES.update(TABLE_B2)
ALL_TABLES.update(TABLE_B3)
ALL_TABLES.update(TABLE_B4)
ALL_TABLES.update(TABLE_B5)
ALL_TABLES.update(TABLE_B6)
ALL_TABLES.update(TABLE_B7)
ALL_TABLES.update(TABLE_B8)

__all__ = [
    'TABLE_B1', 'TABLE_B2', 'TABLE_B3', 'TABLE_B4',
    'TABLE_B5', 'TABLE_B6', 'TABLE_B7', 'TABLE_B8',
    'ALL_TABLES'
]

def get_table_entry(irrep, jp, lp, np, j, l, n):
    # Get natrix element form combined table
    key = (irrep, jp, lp, np, j, l, n)
    return ALL_TABLES.get(key, None)

def get_irrep_table(irrep):
    # Enteries for specefic irreps
    return {k: v for k, v in ALL_TABLES.items() if k[0] == irrep}

if __name__ == "__main__":
    for table_name in ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8']:
        table = globals()[f'TABLE_{table_name}']
        print(f"  Table {table_name}: {len(table)} entries")
    
    print(f"\nTotal entries: {len(ALL_TABLES)}")
