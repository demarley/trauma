# Created: 3 March 2018
# 
# Dan Marley
# daniel.edison.marley@cernSPAMNOT.ch
# Texas A&M University
# ---
# Tar the current CMSSW directory for running condor jobs


# set path for where to put the tarball
set eos_path = "/store/user/demarley/correlator/trauma/flatNtuples/"

echo " Make tarball ..."

tar --exclude-caches-all --exclude-vcs -zcf $CMSSW_VERSION.tgz -C $CMSSW_BASE/.. $CMSSW_VERSION \
    --exclude=$CMSSW_VERSION/src/trauma/production \
    --exclude=$CMSSW_VERSION/src/trauma/miniAna/batch \
    --exclude=$CMSSW_VERSION/src/trauam/miniAna/plots \
    --exclude=$CMSSW_VERSION/src/DataFormats \
    --exclude=$CMSSW_VERSION/src/L1Trigger \
    --exclude=$CMSSW_VERSION/include/slc6_amd64_gcc630/DataFormats \
    --exclude=$CMSSW_VERSION/tmp/slc6_amd64_gcc630/src/DataFormats \
    --exclude=$CMSSW_VERSION/include/slc6_amd64_gcc630/L1Trigger \
    --exclude=$CMSSW_VERSION/tmp/slc6_amd64_gcc630/src/L1Trigger \
    --verbose


echo " Copy tarball to new location on EOS "$eos_path
xrdcp $CMSSW_VERSION.tgz root://cmseos.fnal.gov/"$eos_path"
