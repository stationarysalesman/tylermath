"""Biochemistry module. Holds calculations useful to working in a 
    biochemical lab."""

import numpy as np

def mspConcFromAbsorption(a280, a450):
    """Estimate [MSP] based on spec measurements."""
    ed280 = 64500 # dye extinction at 280nm
    ed450 = 16600 # dye extinction at 450nm
    ep280 = 21430 # protein extinction at 280
    factor1 = float(ed280 / ed450)
    numerator = a280 - (a450 * factor1)
    return numerator / ep280

