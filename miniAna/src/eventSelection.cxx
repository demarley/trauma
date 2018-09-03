/*
Created:        3 September 2018
Last Updated:   3 September 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Event Selection script
*/
#include "trauma/miniAna/interface/eventSelection.h"



eventSelection::eventSelection(configuration &cmaConfig, const std::string &level) :
  m_config(&cmaConfig),
  m_level(level),
  m_selection("SetMe"),
  m_cutsfile("SetMe"),
  m_numberOfCuts(0),
  m_training(false){
    m_cuts.resize(0);
    m_cutflowNames.clear();

    m_selection = m_config->selection();
    m_cutsfile  = m_config->cutsfile();
  }

eventSelection::~eventSelection() {}


void eventSelection::initialize() {
    /* Build the cuts using the cut file from configuration */
    initialize( m_cutsfile,m_selection );
    return;
}

void eventSelection::initialize(const std::string& cutsfile, const std::string& selection) {
    /* Load cut values using specific name for cutsfile */
    m_cutsfile  = cutsfile;
    m_selection = selection;

    std::ifstream file = cma::open_file(cutsfile);

    // Read one line at a time into the vector of Cut structs:
    // this only stores information, but can be expanded
    m_cuts.clear();
    std::string line;
    if (file.is_open()){
        while(std::getline(file, line)){
            std::stringstream  lineStream(line);
            Cut tmp_cut;
            lineStream >> tmp_cut.name >> tmp_cut.comparison >> tmp_cut.value;  // read line into struct
            m_cuts.push_back(tmp_cut);
        }
        file.close();
    } // end reading cuts file

    // Get the number of cuts (for cutflow histogram binning)
    m_numberOfCuts = m_cuts.size();

    // Get the names of cuts (for cutflow histogram bin labeling)
    m_cutflowNames.clear();
    getCutNames();

    // Identify the selection this instance will apply
    identifySelection();

    return;
}


void eventSelection::identifySelection(){
    /* Set the booleans for applying the selection below */
    m_training = m_selection.compare("training")==0;
    return;
}


void eventSelection::setCutflowHistograms(TFile& outputFile){
    /* Set the cutflow histograms to use in the framework -- 
       can modify this function to generate histograms with different names
       e.g., based on the name of the TTree 
       Two cutflows:  
         "cutflow"            event weights
         "cutflow_unweighted" no event weights -> raw event numbers
    */
    outputFile.cd();

    m_cutflow     = new TH1D( (m_selection+"_cutflow").c_str(),(m_selection+"_cutflow").c_str(),m_numberOfCuts+1,0,m_numberOfCuts+1);
    m_cutflow_unw = new TH1D( (m_selection+"_cutflow_unweighted").c_str(),(m_selection+"_cutflow_unweighted").c_str(),m_numberOfCuts+1,0,m_numberOfCuts+1);

    m_cutflow->GetXaxis()->SetBinLabel(1,"INITIAL");
    m_cutflow_unw->GetXaxis()->SetBinLabel(1,"INITIAL");

    for (unsigned int c=1;c<=m_numberOfCuts;++c){
        m_cutflow->GetXaxis()->SetBinLabel(c+1,m_cutflowNames.at(c-1).c_str());
        m_cutflow_unw->GetXaxis()->SetBinLabel(c+1,m_cutflowNames.at(c-1).c_str());
    }

    return;
}


bool eventSelection::applySelection(const Event &event) {
    /* Apply cuts 
       Example Cut:
          if (n_jets==3 && n_ljets<1)  FAIL (return false)
          else :                       PASS & fill cutflows (continue to next cut)
    */
    m_nominal_weight = event.nominal_weight();
    double first_bin(0.5);            // first bin value in cutflow histogram ("INITIAL")
                                      // easily increment by 1 for each cut (don't need to remember bin number)

    if(!event.isValidRecoEntry()) return false;  // check event is valid

    // fill cutflow histograms with initial value (before any cuts)
    fillCutflows(first_bin);

    std::vector<TkMu> tkmus = event.trackmuons();

    if (m_training){
        // Cut 1 :: >=1 Track-Muon pair
        if (tkmus.size()<1)
            return false;
        else
            fillCutflows(first_bin+1);
    }

    return true;
}



void eventSelection::fillCutflows(double cutflow_bin){
    /* Fill cutflow histograms with weight at specific bin */
    m_cutflow->Fill(cutflow_bin,m_nominal_weight);
    m_cutflow_unw->Fill(cutflow_bin);
    return;
}


void eventSelection::getCutNames(){
    /* Get the cut names (for labeling bins in cutflow histograms) and store in vector */
    m_cutflowNames.clear();
    for (const auto& cut : m_cuts)
        m_cutflowNames.push_back( cut.name );

    return;
}

// THE END //
