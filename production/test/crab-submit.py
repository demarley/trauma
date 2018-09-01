"""
Created 1 September 2018

Script to submit multiple CRAB jobs
 https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial#1_CRAB_configuration_file_to_run

To run:
  cd test/
  python crab-submit-multiple.py <datasets>

where <datasets> is a text file that contains the datasets you want to process.
See 'test/crab-datasets-mc.txt' for an example.
If no argument is provided, a default option is selected
"""
import os
import sys
from multiprocessing import Process



def file2list(filename):
    """Load text file and dump contents into a list"""
    listOfFiles = open( filename,'r').readlines()
    listOfFiles = [i.rstrip('\n') for i in listOfFiles if not i.startswith("#")]
    return listOfFiles


def main(input_datasets="crab-datasets.txt"):

    from CRABClient.UserUtilities import config
    config = config()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    input_datasets = file2list(input_datasets)


    def submit(config):
        try:
            crabCommand('submit', config=config, dryrun=False)
            print ' Executed crabCommand() '
        except HTTPException, hte:
            print ' ERROR :: Cannot execute command! '
            print hte.headers


    for id in input_datasets:
        if (not id or id.startswith('#')): continue

        name,dataset = id.split(" ")
        primary_dataset = dataset.split('/')[1]

        # General
        config.General.requestName = 'trauma_'+name
        config.General.workArea    = 'crab_trauma_'+name
        config.General.transferOutputs = True
        config.General.transferLogs    = True

        # JobType
        config.JobType.pluginName  = 'Analysis'
        config.JobType.psetName    = 'trauma.py'
        config.JobType.allowUndistributedCMSSW = True

        # Data
        #config.Data.splitting    = 'Automatic'
        config.Data.splitting     = 'FileBased'
        config.Data.unitsPerJob   = 3
        config.Data.outLFNDirBase = '/store/user/dmarley/correlator/'
        config.Data.publication   = False
        config.Data.inputDataset  = dataset

        # Site
        config.Site.storageSite   = "T3_US_FNALLPC"

        print '\n Configuration :'
        print config
        try :
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
        except :
            print ' ERROR :: Not submitted!'


if __name__=='__main__':
    try:
        main(sys.argv[1])    # pass datasets file and year as command line arguments
    except IndexError:
        print " Not enough arguments! "
        main()

## THE END ##
