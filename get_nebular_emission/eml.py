from get_nebular_emission.eml_io import get_data

def eml(infile, m_sfr_z=[0,1,2], h0=None , verbose=False):
    mstars, ssfr, loh12 = get_data(infile, m_sfr_z, h0=h0, verbose=verbose)
    print(mstars)

