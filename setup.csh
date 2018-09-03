echo ""
echo " > Setup CMS "
cmsenv

echo " > Setup CRAB "
source /cvmfs/cms.cern.ch/crab3/crab.csh

# Set grid proxy if not provided
echo ""
echo " > Setup Grid "
voms-proxy-info -exists > /dev/null
setenv global_proxy_ok $?
if ( ${global_proxy_ok} != 0 ) then
    echo "   - No valid grid proxy found. Creating new one."
    voms-proxy-init -voms cms
    if ( $? != 0 ) then
        echo "   - Failed to create grid proxy."
    endif
else
    echo "   - Grid proxy already set."
endif
