import numpy as np

notnum = -999.

zsun = 0.0134 # Asplund 2009
# zsum = 0.014 #CMB constant ref ?
ohsun = 8.69  # Allende Prieto 2001 and Asplund 2009 (12 + log10(O/H))sun

vol_pm = 542.16**3. # Volume of the P-Millenium simulation

h = 0.6777 # Hubble constant H = 100 h km s**-1 Mpc**-1

unemods = ['kashino20']

photmods = ['gutkin16']

zmet = {
    "gutkin16" : np.array([0.0001,
                         0.0002,
                         0.0005,
                         0.001,
                         0.002,
                         0.004,
                         0.006,
                         0.008,
                         0.010,
                         0.014,
                         0.017,
                         0.020,
                         0.030,
                         0.040])
}

zmet_reduced = {
    "gutkin16" : np.array([0.0001,
                         0.002,
                         0.014,
                         0.030])
}
