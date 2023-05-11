#ifndef gv_def
#define gv_def

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

namespace gv
{
    // extern StPicoEvent *picoEvent;
    // extern StPicoDst *picoDst;

    //Structure 2 - particles
    extern particle *p_lambda;
    extern particle *p_lambda_rot;
    extern particle *p_proton;
    extern particle *p_pion;
    extern particle_all *p_all;
    extern particle_dau *p_dau;
    extern particle_dau *p_dau_rot;

    //Structure 3 - Event details
    extern event *details;

    //Global Variables
//    extern TTree *corr_tree;
    // extern struct histograms_created all_hists;
//     extern struct histograms_xi all_hists_xi;
    
// //    extern TTree *mXiTree;

//     extern int check;
    
    // extern StPicoTrack *track_forp;

    // extern TClonesArray *mEpdHits;
    // extern int check;
}

#endif
