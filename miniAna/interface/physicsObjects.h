#ifndef PHYSICSOBJECTS_H
#define PHYSICSOBJECTS_H

/* 
   Physics objects to be used in analyses
   This structure allows the Event class
   and other classes to access these objects
*/
#include "TLorentzVector.h"
#include <map>
#include <string>


// base object (consistent reference to TLorentzVector)
struct TBase {
    TLorentzVector p4;
    bool isGood;
    int charge;
    int index;
};


// Extra muon attributes
struct Muon : TBase{
    bool OMTF;
    bool EMTF;
    bool BMTF;
};


// Extra TT track attributes
struct Track : TBase{
    float d0;
    float z0;
    float chi2;
    int nstub;

    int genuine;
    int loose;
    int unknown;
    int combinatoric;
    int fake;           // = {0=,1=,2=}

    // tp match
    int matchtp_pdgid;
    float matchtp_pt;
    float matchtp_eta;
    float matchtp_phi;
    float matchtp_z0;
    float matchtp_dxy;
};


// Track-Muon Object
struct TkMu : TBase{
    // Reference indices in vector<Track> & vector<Muon> 
    int track;
    int muon;
    std::map<std::string,double> features;
    float dnn;
};

#endif
