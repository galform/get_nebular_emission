import get_nebular_emission.eml_io as eml

exdir = 'example_data/'
exfile = exdir+'example_data.dat'

def test_check_file():
    assert eml.check_file(exfile) is True

def test_cdreate_dir():
    assert eml.create_dir(exdir) is True

    
def test_get_nheader():
    assert eml.get_nheader(exfile) == 6
