"""
Extra labels for plots
"""
import os
import sys
from array import array


def hist1d(nbins,bin_low,bin_high):
    """
    Set the binning for a given histogram.
    @param nbins	  Number of bins in histogram
    @param bin_low    Lower bin edge
    @param bin_high   Upper bin edge
    """
    binsize = float(bin_high-bin_low)/nbins
    arr     = array('d',[i*binsize+bin_low for i in xrange(nbins+1)])
    return arr


class Sample(object):
    """Class for organizing plotting information about physics samples"""
    def __init__(self,label='',color=''):
        self.label = label
        self.color = color

class Variable(object):
    """Class for organizing plotting information about variables"""
    def __init__(self,binning=[],label=''):
        self.binning = binning
        self.label   = label


def variable_labels():
    """Dictionaries that contain Variables objects."""
    _phi  = r'$\phi$'
    _eta  = r'$\eta$'
    _T    = r'$_{\text{T}}$ [GeV]'

    variables = {}

    variables["mu_pt"]  = Variable(binning=hist1d(40,0,200), label=r'Muon p'+_T)
    variables["mu_eta"] = Variable(binning=hist1d(10,-2.5,2.5), label=r'Muon '+_eta)
    variables["mu_phi"] = Variable(binning=hist1d(16,-3.2,3.2), label=r'Muon '+_phi)
    variables["mu_charge"]  = Variable(binning=hist1d(3,-1.5,1.5), label=r'Muon Charge')
    variables["tk_pt"]      = Variable(binning=hist1d(40,0,200), label=r'Track p'+_T)
    variables["tk_eta"]     = Variable(binning=hist1d(10,-2.5,2.5), label=r'Track '+_eta)
    variables["tk_phi"]     = Variable(binning=hist1d(16,-3.2,3.2), label=r'Track '+_phi)
    variables["tk_sinheta"] = Variable(binning=hist1d(10,-5,5), label=r'Track sinh('+_eta+')')
    variables["tk_rinv"]    = Variable(binning=10, label=r'Track Rinv')
    variables["tk_chi2"]    = Variable(binning=hist1d(10,0,200), label=r'Track $\chi^\text{2}$')
    variables["tk_z0"]      = Variable(binning=hist1d(30,-15,15), label=r'Track z$_\text{0}$')
    variables["tk_d0"]      = Variable(binning=hist1d(30,-15,15), label=r'Track d$_\text{0}$')
    variables["tk_charge"]  = Variable(binning=hist1d(3,-1.5,1.5), label=r'Track Charge')
    variables["tkmu_deltaR2"] = Variable(binning=hist1d(20,0,1), label=r'$\Delta$R(track,muon)')

    return variables



def sample_labels():
    """Dictionaries that contain Samples objects.
       > The key values match those in config/sampleMetadata.txt.
         (The functions in util.py are setup to read the information this way.)
         If you want something unique, then you may need to specify 
         it manually in your plotting script
    """
    ## Sample information
    samples = {}

    samples['signal'] = Sample(label='Signal',color='b')
    samples['bckg']   = Sample(label='Bckg',color='r')

    ttbar = r't$\bar{\text{t}}$'
    samples['multijet'] = Sample(label=r'Multi-jet',  color='purple')
    samples['BQ']       = Sample(label=ttbar+' (QB)', color='red')
    samples['W']        = Sample(label=ttbar+' (W)',  color='blue')
    samples['ttbckg']   = Sample(label=ttbar+' bckg.',color='green')

    return samples
