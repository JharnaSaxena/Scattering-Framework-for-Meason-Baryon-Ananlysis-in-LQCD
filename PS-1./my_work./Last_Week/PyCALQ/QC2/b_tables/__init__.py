#!/usr/bin/env python3
"""
Box matrix tables Morningstar Appendix B
All tables imported from individual modules
"""
from .b1 import TABLE_B1
from .b2 import TABLE_B2
from .b3 import TABLE_B3
from .b4 import TABLE_B4
from .b5 import TABLE_B5
from .b6 import TABLE_B6
from .b7 import TABLE_B7
from .b8 import TABLE_B8
#Combine all tables
ALL_TABLES = {
    'B1': TABLE_B1,
    'B2': TABLE_B2,
    'B3': TABLE_B3,
    'B4': TABLE_B4,
    'B5': TABLE_B5,
    'B6': TABLE_B6,
    'B7': TABLE_B7,
    'B8': TABLE_B8,
}

#Table selector based on irrep and spin
IRREP_TO_TABLE = {
    # P=0, S=0 (Tables B.1, B.2)
    'A1g': 'B1',
    'A2g': 'B1',
    'Eg': 'B1',
    'T1g': 'B1',
    'T2g': 'B1',
    'A1u': 'B1',
    'A2u': 'B1',
    'Eu': 'B1',
    'T1u': 'B2',
    'T2u': 'B1',
    
    # P=0, S=1/2 (Tables B.3)
    'G1g': 'B3',
    'G1u': 'B3',
    'G2g': 'B3',
    'G2u': 'B3',
    'Hg': 'B3',
    'Hu': 'B3',
    
    # Moving frames, S=0 (Tables B.4, B.5)
    'G1': 'B4/B5',
    'G2': 'B4/B5',
    'H': 'B4/B5',
    
    # Moving frames, S=1/2 (Tables B.6, B.7, B.8)
    'G1': 'B6', 
    'G2': 'B6', 
    'G': 'B7',  
    'H': 'B7', 
}

def get_table(irrep: str, S: int = 0, d2: int = 0):
    """    
    Argumentss:
    1. irrep: Irreducible representation
    2. S: Total spin (0, 1/2, 1, etc.)
    3. d2: Momentum squared (for moving frames)
    Returns: Table dictionary
    """
    # Determine table
    if S == 0:
        if d2 == 0:
            return ALL_TABLES.get('B1' if irrep not in ['T1u', 'T2u'] else 'B2', {})
        else:
            # Moving frames S=0
            return ALL_TABLES.get('B4' if irrep in ['G1', 'G2'] else 'B5', {})
    elif S == 1:
        # P=0, S=1
        return ALL_TABLES.get('B1', {})  # Table B.1 also has S=1
    elif S == 0.5:
        if d2 == 0:
            return ALL_TABLES.get('B3', {})
        elif d2 == 1:
            return ALL_TABLES.get('B6', {})
        elif d2 == 2:
            return ALL_TABLES.get('B7', {})
        elif d2 == 3:
            return ALL_TABLES.get('B8', {})
    
    return {}

__all__ = [
    'TABLE_B1',
    'TABLE_B2',
    'TABLE_B3',
    'TABLE_B4',
    'TABLE_B5',
    'TABLE_B6',
    'TABLE_B7',
    'TABLE_B8',
    'ALL_TABLES',
    'IRREP_TO_TABLE',
    'get_table'
]
