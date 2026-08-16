#!/usr/bin/env python3
"""
Table B.7: Lambda_B = A2 partial sector, P=(0,n,n), S=0
From arXiv:1707.05817
"""

import numpy as np
#############################################################################
TABLE_B7 = {
    # Row Block 1: (2,2,1)
    ('A2', 2, 2, 1, 2, 2, 1): [
        (1.0, 0, 0, False),
        (-2.0*np.sqrt(5)/7, 2, 0, False),
        (-4.0/7.0, 4, 0, False),
        (2.0*np.sqrt(10)/7, 4, 2, False),
        (np.sqrt(30)/7, 2, 1, True),
        (-8.0*np.sqrt(5)/7, 4, 1, True),
    ],
    ('A2', 2, 2, 1, 3, 3, 1): [
        (-np.sqrt(42)/7, 1, 0, True),
        (-np.sqrt(2)/6, 3, 0, True),
        (-np.sqrt(15)/3, 3, 2, True),
        (25.0*np.sqrt(154)/462, 5, 0, True),
        (-2.0*np.sqrt(165)/33, 5, 2, True),
        (-np.sqrt(55)/11, 5, 4, True),
    ],
    ('A2', 2, 2, 1, 4, 4, 1): [
        (-np.sqrt(105)/7, 2, 0, False),
        (-2.0*np.sqrt(210)/77, 4, 2, False),
        (-32.0*np.sqrt(65)/143, 6, 2, False),
        (-2.0*np.sqrt(105)/11, 4, 1, True),
        (15.0*np.sqrt(26)/143, 6, 1, True),
        (-18.0*np.sqrt(65)/143, 6, 3, True),
    ],
    ('A2', 2, 2, 1, 4, 4, 2): [
        (-np.sqrt(15)/7, 2, 0, True),
        (20.0*np.sqrt(3)/77, 4, 0, True),
        (4.0*np.sqrt(30)/77, 4, 2, True),
        (-40.0*np.sqrt(39)/143, 6, 0, True),
        (-16.0*np.sqrt(455)/143, 6, 2, True),
        (-2.0*np.sqrt(10)/7, 2, 1, True),
        (26.0*np.sqrt(15)/77, 4, 1, True),
        (5.0*np.sqrt(182)/143, 6, 1, True),
        (-6.0*np.sqrt(455)/143, 6, 3, True),
    ],
    ('A2', 2, 2, 1, 5, 5, 1): [
        (5.0*np.sqrt(22)/33, 3, 0, True),
        (-2.0*np.sqrt(165)/33, 3, 2, True),
        (5.0*np.sqrt(14)/156, 5, 0, True),
        (np.sqrt(15)/39, 5, 2, True),
        (-11.0*np.sqrt(5)/26, 5, 4, True),
        (-7.0*np.sqrt(2310)/572, 7, 0, True),
        (51.0*np.sqrt(110)/572, 7, 2, True),
        (-3.0*np.sqrt(5)/26, 7, 4, True),
        (-3.0*np.sqrt(130)/52, 7, 6, True),
    ],
    ('A2', 2, 2, 1, 5, 5, 2): [
        (-5.0*np.sqrt(66)/66, 3, 0, False),
        (-np.sqrt(55)/11, 3, 2, False),
        (5.0*np.sqrt(42)/39, 5, 0, False),
        (2.0*np.sqrt(5)/13, 5, 2, False),
        (-2.0*np.sqrt(15)/13, 5, 4, False),
        (-3.0*np.sqrt(770)/286, 7, 0, False),
        (-np.sqrt(330)/286, 7, 2, False),
        (3.0*np.sqrt(15)/13, 7, 4, False),
        (np.sqrt(390)/26, 7, 6, False),
    ],
    ('A2', 2, 2, 1, 6, 6, 1): [
        (15.0*np.sqrt(10010)/2002, 4, 0, False),
        (-30.0*np.sqrt(1001)/1001, 4, 2, False),
        (np.sqrt(770)/55, 6, 0, False),
        (16.0*np.sqrt(66)/165, 6, 2, False),
        (-3.0*np.sqrt(170170)/2431, 8, 0, False),
        (32.0*np.sqrt(4862)/2431, 8, 2, False),
        (-19.0*np.sqrt(1105)/1105, 8, 4, False),
        (2.0*np.sqrt(165)/55, 6, 1, True),
        (-np.sqrt(66)/11, 6, 3, True),
        (7.0*np.sqrt(85085)/1105, 8, 1, True),
        (np.sqrt(663)/13, 8, 3, True),
    ],
    ('A2', 2, 2, 1, 6, 6, 2): [
        (5.0*np.sqrt(546)/143, 4, 2, False),
        (-8.0/55.0, 6, 2, False),
        (-4.0*np.sqrt(23205)/1105, 8, 0, False),
        (6.0*np.sqrt(72930)/1105, 8, 4, False),
        (-20.0*np.sqrt(273)/143, 4, 1, True),
        (9.0*np.sqrt(10)/55, 6, 1, True),
        (-2.0/11.0, 6, 3, True),
        (-22.0*np.sqrt(46410)/1105, 8, 1, True),
        (-10.0*np.sqrt(4862)/221, 8, 3, True),
    ],
    ('A2', 2, 2, 1, 6, 6, 3): [
        (5.0*np.sqrt(182)/154, 4, 0, True),
        (2.0*np.sqrt(455)/1001, 4, 2, True),
        (-5.0*np.sqrt(14)/11, 6, 0, True),
        (-16.0*np.sqrt(30)/55, 6, 2, True),
        (np.sqrt(3094)/221, 8, 0, True),
        (-3.0*np.sqrt(2431)/221, 8, 4, True),
        (8.0*np.sqrt(910)/143, 4, 1, True),
        (2.0*np.sqrt(3)/11, 6, 1, True),
        (np.sqrt(30)/11, 6, 3, True),
        (-3.0*np.sqrt(1547)/221, 8, 1, True),
        (np.sqrt(36465)/221, 8, 3, True),
    ],
    
    # Row Block 2: (3,3,1)
    ('A2', 3, 3, 1, 3, 3, 1): [
        (1.0, 0, 0, False),
        (-2.0/11.0, 4, 0, False),
        (5.0*np.sqrt(10)/11, 4, 2, False),
        (-60.0*np.sqrt(13)/143, 6, 0, False),
        (-8.0*np.sqrt(1365)/429, 6, 2, False),
    ],
    ('A2', 3, 3, 1, 4, 4, 1): [
        (-np.sqrt(2)/2, 1, 0, True),
        (-np.sqrt(42)/11, 3, 0, True),
        (5.0*np.sqrt(66)/572, 5, 0, True),
        (-6.0*np.sqrt(385)/143, 5, 2, True),
        (5.0*np.sqrt(1155)/286, 5, 4, True),
        (119.0*np.sqrt(10)/572, 7, 0, True),
        (-3.0*np.sqrt(210)/52, 7, 2, True),
        (-np.sqrt(1155)/286, 7, 4, True),
        (np.sqrt(30030)/572, 7, 6, True),
    ],
    ('A2', 3, 3, 1, 4, 4, 2): [
        (-3.0*np.sqrt(14)/14, 1, 0, True),
        (2.0*np.sqrt(6)/11, 3, 0, True),
        (-2.0*np.sqrt(5)/11, 3, 2, True),
        (-5.0*np.sqrt(462)/364, 5, 0, True),
        (-20.0*np.sqrt(55)/143, 5, 2, True),
        (np.sqrt(165)/286, 5, 4, True),
        (21.0*np.sqrt(70)/572, 7, 0, True),
        (7.0*np.sqrt(30)/572, 7, 2, True),
        (-21.0*np.sqrt(165)/286, 7, 4, True),
        (-7.0*np.sqrt(4290)/572, 7, 6, True),
    ],
}

__all__ = ['TABLE_B7']
