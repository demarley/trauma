#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "TROOT.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <sstream>
#include <stdlib.h>

#include "trauma/miniAna/interface/tools.h"


class configuration {
  public:
    // Default - so root can load based on a name;
    configuration( const std::string &configFile );
    //configuration( const configuration& );
    configuration& operator=( const configuration& rhs );

    // Default - so we can clean up;
    virtual ~configuration();

    // Run once at the start of the job;
    virtual void initialize();
    std::string getConfigOption( std::string item );

    // Print configuration
    virtual void print();

    virtual bool isMC() const {return m_isMC;}

    // functions about the TTree
    void setTreename(std::string treeName);
    std::string treename() const{ return m_treename;}

    // functions about the file
    std::vector<std::string> filesToProcess() const{ return m_filesToProcess;}
    void setFilename(std::string fileName);
    std::string filename() const{ return m_filename;}
    std::string primaryDataset() const{return m_primaryDataset;}
    unsigned int NTotalEvents() const{return m_NTotalEvents;}

    // return some values from config file
    bool DNNtraining() const{ return m_DNNtraining;}
    bool DNNinference() const{ return m_DNNinference;}
    std::string dnnFile(){ return m_dnnFile;}
    std::string dnnKey(){ return m_dnnKey;}
    std::string verboseLevel() const{ return m_verboseLevel;}
    std::string selection() const{ return m_selection;}
    std::string cutsfile() const{ return m_cutsfile;}
    std::string outputFilePath() const{ return m_outputFilePath;}
    std::string customFileEnding() const{ return m_customFileEnding;}
    std::string configFileName() const{ return m_configFile;}
    std::string getAbsolutePath() const{ return m_cma_absPath;}
    int nEventsToProcess() const{ return m_nEventsToProcess;}
    unsigned long long firstEvent() const{ return m_firstEvent;}
    bool makeHistograms() const{ return m_makeHistograms;}

  protected:

    std::map<std::string,std::string> m_map_config;
    const std::string m_configFile;

    bool m_isMC;

    // return some values from config file
    bool m_DNNtraining;
    bool m_DNNinference;
    std::string m_dnnFile;
    std::string m_dnnKey;
    std::string m_selection;
    std::string m_cutsfile;
    std::string m_treename;
    std::string m_filename;
    std::string m_primaryDataset;
    unsigned int m_NTotalEvents;
    std::string m_verboseLevel;
    int m_nEventsToProcess;
    unsigned long long m_firstEvent;
    std::string m_outputFilePath;
    std::string m_customFileEnding;
    bool m_makeHistograms;
    std::string m_cma_absPath;
    std::vector<std::string> m_filesToProcess;


    std::map<std::string,std::string> m_defaultConfigs = {
             {"DNNtraining",           "false"},
             {"DNNinference",          "false"},
             {"DNNfile",               ""},
             {"DNNkey",                "dnn"},
             {"makeHistograms",        "true"},
             {"NEvents",               "-1"},
             {"firstEvent",            "0"},
             {"selection",             "none"},
             {"output_path",           "./"},
             {"customFileEnding",      ""},
             {"cutsfile",              "config/cuts_none.txt"},
             {"inputfile",             "config/inputfiles.txt"},
             {"treename",              "stopTreeMaker/AUX"},
             {"metadataFile",          "config/sampleMetaData.txt"},
             {"verboseLevel",          "INFO"} };
};

#endif
