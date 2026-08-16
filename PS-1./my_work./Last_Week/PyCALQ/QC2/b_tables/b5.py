#!/usr/bin/env python3
"""
Table B.5: Lambda_B = B1 partial sector, P=(0,0,n), S=0
From arXiv:1707.05817
"""

import numpy as np

TABLE_B5 = {
    ('B1', 2, 2, 1, 2, 2, 1): [
        (1.0, 0, 0, False),
        (-2.0*np.sqrt(5)/7, 2, 0, False),
        (1.0/7.0, 4, 0, False),
        (np.sqrt(70)/7, 4, 4, False),
    ],
    ('B1', 2, 2, 1, 3, 3, 1): [
        (np.sqrt(21)/7, 1, 0, False),
        (-2.0/3.0, 3, 0, False),
        (5.0*np.sqrt(77)/231, 5, 0, False),
        (np.sqrt(110)/11, 5, 4, False),
    ],
    ('B1', 2, 2, 1, 4, 4, 1): [
        (np.sqrt(15)/7, 2, 0, False),
        (-30.0*np.sqrt(3)/77, 4, 0, False),
        (-2.0*np.sqrt(210)/77, 4, 4, False),
        (5.0*np.sqrt(39)/143, 6, 0, False),
        (5.0*np.sqrt(546)/143, 6, 4, False),
    ],
    ('B1', 2, 2, 1, 5, 5, 1): [
        (5.0*np.sqrt(11)/33, 3, 0, False),
        (-10.0*np.sqrt(7)/39, 5, 0, False),
        (-2.0*np.sqrt(10)/13, 5, 4, False),
        (np.sqrt(1155)/143, 7, 0, False),
        (3.0*np.sqrt(10)/13, 7, 4, False),
    ],
    ('B1', 2, 2, 1, 6, 6, 1): [
        (15.0*np.sqrt(143)/143, 4, 4, False),
        (-2.0*np.sqrt(55)/55, 6, 4, False),
        (np.sqrt(1105)/1105, 8, 4, False),
        (2.0*np.sqrt(119)/17, 8, 8, False),
    ],
    ('B1', 2, 2, 1, 6, 6, 2): [
        (5.0*np.sqrt(182)/143, 4, 0, False),
        (np.sqrt(65)/143, 4, 4, False),
        (-2.0*np.sqrt(14)/11, 6, 0, False),
        (-6.0/11.0, 6, 4, False),
        (np.sqrt(3094)/221, 8, 0, False),
        (3.0*np.sqrt(2431)/221, 8, 4, False),
    ],
    ('B1', 3, 3, 1, 3, 3, 1): [
        (1.0, 0, 0, False),
        (-7.0/11.0, 4, 0, False),
        (np.sqrt(70)/11, 4, 4, False),
        (10.0*np.sqrt(13)/143, 6, 0, False),
        (10.0*np.sqrt(182)/143, 6, 4, False),
    ],
    ('B1', 3, 3, 1, 4, 4, 1): [
        (2.0*np.sqrt(7)/7, 1, 0, False),
        (-np.sqrt(3)/11, 3, 0, False),
        (-40.0*np.sqrt(231)/1001, 5, 0, False),
        (4.0*np.sqrt(330)/143, 5, 4, False),
        (7.0*np.sqrt(35)/143, 7, 0, False),
        (7.0*np.sqrt(330)/143, 7, 4, False),
    ],
    ('B1', 3, 3, 1, 5, 5, 1): [
        (np.sqrt(55)/11, 2, 0, False),
        (-10.0*np.sqrt(11)/143, 4, 0, False),
        (-2.0*np.sqrt(770)/143, 4, 4, False),
        (-7.0*np.sqrt(143)/143, 6, 0, False),
        (np.sqrt(2002)/143, 6, 4, False),
        (56.0*np.sqrt(187)/2431, 8, 0, False),
        (12.0*np.sqrt(238)/221, 8, 4, False),
    ],
    ('B1', 3, 3, 1, 6, 6, 1): [
        (np.sqrt(91)/13, 5, 4, False),
        (-6.0*np.sqrt(91)/221, 7, 4, False),
        (np.sqrt(133)/323, 9, 4, False),
        (14.0*np.sqrt(323)/323, 9, 8, False),
    ],
    ('B1', 3, 3, 1, 6, 6, 2): [
        (20.0*np.sqrt(182)/429, 3, 0, False),
        (-7.0*np.sqrt(286)/429, 5, 0, False),
        (-np.sqrt(5005)/143, 5, 4, False),
        (-70.0*np.sqrt(390)/2431, 7, 0, False),
        (6.0*np.sqrt(5005)/2431, 7, 4, False),
        (63.0*np.sqrt(494)/4199, 9, 0, False),
        (3.0*np.sqrt(7315)/323, 9, 4, False),
    ],
    ('B1', 4, 4, 1, 4, 4, 1): [
        (1.0, 0, 0, False),
        (8.0*np.sqrt(5)/77, 2, 0, False),
        (-27.0/91.0, 4, 0, False),
        (81.0*np.sqrt(70)/1001, 4, 4, False),
        (-2.0*np.sqrt(13)/13, 6, 0, False),
        (6.0*np.sqrt(182)/143, 6, 4, False),
        (196.0*np.sqrt(17)/2431, 8, 0, False),
        (42.0*np.sqrt(2618)/2431, 8, 4, False),
    ],
    ('B1', 4, 4, 1, 5, 5, 1): [
        (np.sqrt(77)/11, 1, 0, False),
        (2.0*np.sqrt(33)/143, 3, 0, False),
        (-np.sqrt(21)/13, 5, 0, False),
        (np.sqrt(30)/13, 5, 4, False),
        (-64.0*np.sqrt(385)/2431, 7, 0, False),
        (20.0*np.sqrt(30)/221, 7, 4, False),
        (252.0*np.sqrt(4389)/46189, 9, 0, False),
        (42.0*np.sqrt(7410)/4199, 9, 4, False),
    ],
    ('B1', 4, 4, 1, 6, 6, 1): [
        (-4.0*np.sqrt(429)/143, 4, 4, False),
        (9.0*np.sqrt(165)/187, 6, 4, False),
        (-18.0*np.sqrt(3315)/4199, 8, 4, False),
        (-12.0*np.sqrt(357)/323, 8, 8, False),
        (7.0*np.sqrt(3)/323, 10, 4, False),
        (42.0*np.sqrt(17)/323, 10, 8, False),
    ],
    ('B1', 4, 4, 1, 6, 6, 2): [
        (2.0*np.sqrt(2730)/143, 2, 0, False),
        (-4.0*np.sqrt(195)/143, 4, 4, False),
        (-np.sqrt(42)/17, 6, 0, False),
        (23.0*np.sqrt(3)/187, 6, 4, False),
        (-18.0*np.sqrt(9282)/3553, 8, 0, False),
        (222.0*np.sqrt(7293)/46189, 8, 4, False),
        (315.0*np.sqrt(26)/4199, 10, 0, False),
        (21.0*np.sqrt(165)/323, 10, 4, False),
    ],
    ('B1', 5, 5, 1, 5, 5, 1): [
        (1.0, 0, 0, False),
        (2.0*np.sqrt(5)/13, 2, 0, False),
        (-1.0/13.0, 4, 0, False),
        (np.sqrt(70)/13, 4, 4, False),
        (-24.0*np.sqrt(13)/221, 6, 0, False),
        (8.0*np.sqrt(182)/221, 6, 4, False),
        (-28.0*np.sqrt(17)/247, 8, 0, False),
        (42.0*np.sqrt(2618)/4199, 8, 4, False),
        (360.0*np.sqrt(21)/4199, 10, 0, False),
        (36.0*np.sqrt(10010)/4199, 10, 4, False),
    ],
}

__all__ = ['TABLE_B5']
