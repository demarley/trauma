/*
Created:         3 September 2018
Last Updated:    3 September 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Configuration class
  -- Read config file and use functions
     to return configurations later
*/
#include "trauma/miniAna/interface/configuration.h"


configuration::configuration(const std::string &configFile) : 
  m_configFile(configFile),
  m_isMC(true),
  m_DNNtraining(false),
  m_DNNinference(false),
  m_dnnFile("SetMe"),
  m_dnnKey("SetMe"),
  m_selection("SetMe"),
  m_cutsfile("SetMe"),
  m_treename("SetMe"),
  m_filename("SetMe"),
  m_primaryDataset("SetMe"),
  m_NTotalEvents(0),
  m_verboseLevel("SetMe"),
  m_nEventsToProcess(0),
  m_firstEvent(0),
  m_outputFilePath("SetMe"),
  m_customFileEnding("SetMe"),
  m_makeHistograms(false),
  m_cma_absPath("SetMe"){
    m_map_config.clear();
  }

configuration::~configuration() {}

configuration &configuration::operator=(const configuration &rhs) { return *this; }

void configuration::initialize() {
    /* Initialize the configurations */
    std::vector<std::string> configurations; 
    cma::read_file( m_configFile, configurations ); // read config file into vector

    // fill map with values from configuration file
    for (const auto& config : configurations){
        // split config items by space
        std::istringstream cfg(config);
        std::istream_iterator<std::string> start(cfg), stop;
        std::vector<std::string> tokens(start, stop);

        m_map_config.insert( std::pair<std::string,std::string>(tokens.at(0),tokens.at(1)) );
    }

    // Protection against default settings missing in custom configuration
    // -- map of defaultConfigs defined in header (can't use 'verbose' tools, not defined yet!)
    for (const auto& defaultConfig : m_defaultConfigs){
        if ( m_map_config.find(defaultConfig.first) == m_map_config.end() ){ // item isn't in config file
            std::cout << " WARNING :: CONFIG : Configuration " << defaultConfig.first << " not defined" << std::endl;
            std::cout << " WARNING :: CONFIG : Setting value to default " << defaultConfig.second << std::endl;
            m_map_config[defaultConfig.first] = defaultConfig.second;
        }
    }


    // Set the verbosity level (the amount of output to the console)
    std::map<std::string,unsigned int> verboseMap = cma::verboseMap(); // load mapping of string to integer
    m_verboseLevel = getConfigOption("verboseLevel");
    if (verboseMap.find(m_verboseLevel)==verboseMap.end()){
        m_verboseLevel = "INFO";
        cma::setVerboseLevel(m_verboseLevel);

        cma::WARNING( "CONFIG : Verbose level selected, "+m_verboseLevel+", is not supported " );
        cma::WARNING( "CONFIG : Please select one of the following: " );
        for (const auto& dm : verboseMap)
            cma::WARNING( "CONFIG :          "+dm.first);
        cma::WARNING( "CONFIG : Continuing; setting verbose level to "+m_verboseLevel);
    }
    else{
        cma::setVerboseLevel(m_verboseLevel);
    }


    // Get the absolute path to leopard for loading
    char* cma_path = getenv("PWD");
    m_cma_absPath = cma_path;
    cma::DEBUG("CONFIG : path set to: "+m_cma_absPath );

    // Assign values
    m_DNNtraining      = cma::str2bool( getConfigOption("DNNtraining") );
    m_DNNinference     = cma::str2bool( getConfigOption("DNNinference") );
    m_dnnFile          = getConfigOption("DNNfile");
    m_dnnKey           = getConfigOption("DNNkey");
    m_nEventsToProcess = std::stoi(getConfigOption("NEvents"));
    m_firstEvent       = std::stoi(getConfigOption("firstEvent"));
    m_selection        = getConfigOption("selection");
    m_outputFilePath   = getConfigOption("output_path");
    m_customFileEnding = getConfigOption("customFileEnding");
    m_cutsfile         = getConfigOption("cutsfile");
    m_makeHistograms   = cma::str2bool( getConfigOption("makeHistograms") );

    cma::read_file( getConfigOption("inputfile"), m_filesToProcess );
    m_treename = getConfigOption("treename");

    return;
}


void configuration::print(){
    // -- Print the configuration
    std::cout << " ** Leopard ** " << std::endl;
    std::cout << " --------------- " << std::endl;
    std::cout << " CONFIGURATION :: Printing configuration " << std::endl;
    std::cout << " " << std::endl;
    for (const auto& config : m_map_config){
        std::cout << " " << config.first << "\t\t\t" << config.second << std::endl;
    }
    std::cout << " --------------- " << std::endl;

    return;
}


std::string configuration::getConfigOption( std::string item ){
    /* Check that the item exists in the map & return it; otherwise throw exception  */
    std::string value("");

    try{
        value = m_map_config.at(item);
    }
    catch(const std::exception&){
        cma::ERROR("CONFIG : Option "+item+" does not exist in configuration.");
        cma::ERROR("CONFIG : This does not exist in the default configuration either.");
        cma::ERROR("CONFIG : Returing an empty string.");
    }

    return value;
}


void configuration::setTreename(std::string treeName){
    m_treename = treeName;
    return;
}

void configuration::setFilename(std::string fileName){
    m_filename = fileName;
    return;
}

// THE END
