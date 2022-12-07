# Define base classes that represent any kind of radioactive source


class Sample(object):
    """
    A Sample is a radioactive source that can contain any number of isotopes with given activities
    """
    pass


class Source(Sample):
    """
    A Source is a radiactive source that contains only ONE isotope
    """
    pass

