/// C++ headers
#include <iostream>

/// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
//#include "TVector3.hh"
#include <assert.h>
#include <TTree.h>
#include <vector>
#include <TLorentzVector.h>
#include <TClonesArray.h>

/// PicoDst headers
// #include "../StPicoEvent/StPicoDst.h"
// #include "../StPicoEvent/StPicoEvent.h"
// #include "../StPicoEvent/StPicoTrack.h"
// #include "../StRefMultCorr/StRefMultCorr.h"
// #include "../StRefMultCorr/CentralityMaker.h"
#include "./particle.h"
#include "./particle_all.h"
#include "./particle_dau.h"
#include "./event.h"

/// Analysis headers
// #include "./StV0Type.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/structs.h"
#include "./gv.h"

namespace gv{
	
	//Structure 2 - particles
    particle *p_lambda = new particle[250];
    particle *p_lambda_rot = new particle[250];
    particle_all *p_all = new particle_all[2000];
    particle_dau *p_dau = new particle_dau[500];
    particle_dau *p_dau_rot = new particle_dau[500];
    particle *p_proton = new particle[1000];
    particle *p_pion = new particle[2000];

    //Structure 3 - Event details
    event *details = new event();

// //	TTree *corr_tree = new TTree("corr_tree", "Tree for Correlations");
// //    TTree *mXiTree = new TTree("XiPicoDst", "XiPicoDst from StXiMaker");

//     int check = 1;
    
    // struct histograms_created all_hists;
//     struct histograms_xi all_hists_xi;
    
    // StPicoEvent *picoEvent;
    // StPicoDst *picoDst;
    
    // StPicoTrack *track_forp;

    // TClonesArray *mEpdHits;
    // int check = 1;
}
