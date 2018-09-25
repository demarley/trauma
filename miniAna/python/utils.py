"""
Created:        --
Last Updated:    2 March 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

File that holds any and all misc. functions 
to be called from other python scripts.
(All information in one file => one location to update!)
"""


def read_config(filename,separation=" "):
    """
    Read configuration file with data stored like:
       'config option'
    And the 'config' and 'option' are separated by a character, e.g., " "
    """
    data = file2list(filename)
    cfg = {}
    for i in data:
        j = i.split(separation)
        cfg[j[0]] = j[1]
    return cfg


def str2bool(param):
    """Convert a string to a boolean"""
    return (param in ['true','True','1'])


## THE END ##
