#ifndef StKFParticleAnalysisMaker_h
#define StKFParticleAnalysisMaker_h

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"

#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "StMaker.h"
#include "TString.h"
#include "TObject.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "TF1.h"
#include "StTrackHelix.h"
#include "MyConstant.h" // must include
// #include "Gamma_112_module.C"

#include "particle.h"
#include "particle_all.h"
#include "particle_dau.h"
#include "event.h"
#include "gv.h"
#include "./namespaces/gv_gamma.h"

class StPicoDst;
class StPicoDstMaker;
class TString;
class KFParticle;
class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class TH1F;
class TH2F;
class TH2D;
class TH3F;
class TH3D;
class TF1;
class TProfile;
class TProfile2D;
class TProfile3D;
class CentralityMaker;
class StRefMultCorr;
class TTree;
// class StMyAna

class StKFParticleAnalysisMaker : public StMaker
{
public:
    StKFParticleAnalysisMaker(const char *name, const char *outName);
    virtual ~StKFParticleAnalysisMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt = "");
    virtual Int_t Finish();

    void    setRunEnergyAndListDir(int run, double energy, char ListDir[256]);
    void set_cen_curr(int c);
    //    void set_lambdatype_curr(TString lam_type);
    void set_JobIDName_curr(TString Job_ID);
    void set_opt_weight_curr(int opt_weight);
    void set_sys_err_opt_curr(int sys_err_opt);
    void SetEpdHitsArray(TClonesArray *mEpdHitstmp);

    int nLambda;
    int nLambdaRot;
    int nproton;
    int npion;

    int EPD_Weight_Read = 0;
    // StEpdEpFinder *mEpFinder;
    std::vector<StEpdEpFinder *> mEpFinder;
    TProfile2D *Read_BBC_EP_east = 0, *Read_BBC_EP_west = 0, *Read_EPD_EP_east = 0, *Read_EPD_EP_west = 0, *Read_EPD_EP1_east = 0, *Read_EPD_EP1_west = 0;
    TH2D *wt = new TH2D("Order1etaWeight", "Order1etaWeight", 100, 1.5, 6.5, 9, 0, 9);
    TH2D *wt2 = new TH2D("Order2etaWeight", "Order2etaWeight", 100, 1.5, 6.5, 9, 0, 9);

private:
    // KFParticle
    StKFParticleInterface *KFParticleInterface;
    StKFParticlePerformanceInterface *KFParticlePerformanceInterface;
    void SetupKFParticle();
    void SetDaughterTrackPointers(int iKFParticle);
    bool InterfaceCantProcessEvent;
    int ProtonTrackIndex, PionTrackIndex;
    vector<int> trackMap;
    StPicoDst *PicoDst;
    StPicoTrack *ProtonTrack, *PionTrack;
    StPicoEvent *mEvent;
    void BookVertexPlots();

    StPicoDstMaker *mPicoDstMaker;
    StRefMultCorr *mRefMultCorr;

    int        mRun;
    int        mcen_curr;
    int        mopt_weight_curr;
    int        msys_err_opt_curr;
    double     mEnergy;
    TString    mListDir;
    //    TString    mlambdatype_curr;
    TString    mJobIDName_curr;
    TClonesArray *mEpdHits_maker;

    TString    mOutName;
    double     PI;
    double     twoPI;

    int        mJob;

    const int order = 10;

    bool debug_3, debug_1, debug_4, debug_2;

    //int cen_curr;
    //TString lambdatype_curr;
    //TString JobIDName_curr;
    //int opt_weight_curr;
    //int sys_err_opt_curr;

    double magnet_global;

    ////////////////
    TH1F *hNRefMult;
    TH1F *hNRefMultA;
    TH1F *hNRefMultB;
    TH2F *hVertexXY_kf;
    TH1F *hVertexZ_kf;
    TH2F *hVertex2D;
    TH1F *hDiffVz  ;
    TH1F *hcent;
    TH1F *hcentw;
    TH1D *hlamdist;
    TH1D *hantilamdist;
    TH1D *hdau_p_azim_dist;

    TProfile *hcentRefM ;
    TProfile *hcentRefW ;

    TH1F *V0Mass_cent[9][17];
    TH1F *V0Mass_anti_cent[9][17];

    TH2D *Ref_TOF_before;
    TH2D *Ref_TOF_after;

    TH1F *pTDist_b4_TOF = new TH1F("pTDist_b4_TOF", "Transverse Momentum of Primary Protons before TOF Cut", 100, 0, 10);
    TH1F *pTDist_after_TOF = new TH1F("pTDist_after_TOF", "Transverse Momentum of Primary Protons after TOF Cut", 100, 0, 10);
    TH2D *m2_proton_vs_pT = new TH2D("m2_proton_vs_pT", "M2 vs. pT of Protons", 200, 0, 20, 100, 0, 4);

    TTree *corr_tree;

    TFile *tree_file;

    event *details_tmp = new event();
    std::vector<particle> *p_lambda_tmp = new vector<particle>;
    std::vector<particle> *p_proton_tmp = new vector<particle>;
    std::vector<particle> *p_pion_tmp = new vector<particle>;
    std::vector<particle_all> *p_all_tmp = new vector<particle_all>;
    std::vector<particle_dau> *p_dau_tmp = new vector<particle_dau>;
    // std::vector<particle> *p_lambda_rot_tmp = new vector<particle>;
    // std::vector<particle_dau> *p_dau_rot_tmp = new vector<particle_dau>;
    

    // Event variable to save in tree
    int dt_cent, dt_num_trk, dt_nLambda, dt_nLambdaRot, dt_Run, dt_TOFMult, dt_RefMult, dt_n_Proton, dt_n_Pion, dt_EventID;
    float dt_VPDvz, dt_PVtxz, dt_PVtxx, dt_PVtxy, dt_Eweight, dt_EPD_EP1_east, dt_EPD_EP1_west, dt_EPD_EP_east, dt_EPD_EP_west, dt_Magn;
    int n_lam, n_dau;

    // Lambda Particle to save in tree
    std::vector<float> *lambda_px = new vector<float>;
    std::vector<float> *lambda_py = new vector<float>;
    std::vector<float> *lambda_pz = new vector<float>;
    std::vector<float> *lambda_Charge = new vector<float>;
    std::vector<float> *lambda_dcaglobal = new vector<float>;
    std::vector<float> *lambda_nsigma = new vector<float>;
    std::vector<float> *lambda_mass = new vector<float>;
    std::vector<int> *lambda_trk_id = new vector<int>;
    std::vector<float> *lambda_hits_ratio = new vector<float>;
    std::vector<float> *lambda_nhitsfit = new vector<float>;
    std::vector<float> *lambda_nhitsmax = new vector<float>;

    // // Proton Particle to save in tree
    // std::vector<float> *proton_px = new vector<float>;
    // std::vector<float> *proton_py = new vector<float>;
    // std::vector<float> *proton_pz = new vector<float>;
    // std::vector<float> *proton_Charge = new vector<float>;
    // std::vector<float> *proton_dcaglobal = new vector<float>;
    // std::vector<float> *proton_nsigma = new vector<float>;
    // std::vector<float> *proton_mass = new vector<float>;
    // std::vector<int> *proton_trk_id = new vector<int>;
    // std::vector<float> *proton_hits_ratio = new vector<float>;
    // std::vector<float> *proton_nhitsfit = new vector<float>;
    // std::vector<float> *proton_nhitsmax = new vector<float>;

    // // Pion Particle to save in tree
    // std::vector<float> *pion_px = new vector<float>;
    // std::vector<float> *pion_py = new vector<float>;
    // std::vector<float> *pion_pz = new vector<float>;
    // std::vector<float> *pion_Charge = new vector<float>;
    // std::vector<float> *pion_dcaglobal = new vector<float>;
    // std::vector<float> *pion_nsigma = new vector<float>;
    // std::vector<float> *pion_mass = new vector<float>;
    // std::vector<int> *pion_trk_id = new vector<int>;
    // std::vector<float> *pion_hits_ratio = new vector<float>;
    // std::vector<float> *pion_nhitsfit = new vector<float>;
    // std::vector<float> *pion_nhitsmax = new vector<float>;

    // All Particles to save in tree
    std::vector<float> *all_px = new vector<float>;
    std::vector<float> *all_py = new vector<float>;
    std::vector<float> *all_pz = new vector<float>;
    std::vector<float> *all_Charge = new vector<float>;
    std::vector<float> *all_dcaglobal = new vector<float>;
    std::vector<float> *all_nSigmaProton = new vector<float>;
    std::vector<float> *all_nSigmaPion = new vector<float>;
    std::vector<bool> *all_isTofTrack = new vector<bool>;
    std::vector<int> *all_trk_id = new vector<int>;
    std::vector<bool> *all_is_pion = new vector<bool>;
    std::vector<bool> *all_is_proton = new vector<bool>;
    std::vector<bool> *all_is_all = new vector<bool>;
    std::vector<float> *all_nhitsmax = new vector<float>;
    std::vector<float> *all_nhitsfit = new vector<float>;
    
    // Daughter Particles to save in tree
    std::vector<float> *dau_px = new vector<float>;
    std::vector<float> *dau_py = new vector<float>;
    std::vector<float> *dau_pz = new vector<float>;
    std::vector<float> *dau_dcaglobal = new vector<float>;
    std::vector<float> *dau_nSigma = new vector<float>;
    std::vector<int> *dau_trk_id = new vector<int>;
    std::vector<float> *dau_nHitsFit = new vector<float>;
    std::vector<float> *dau_nHitsMax = new vector<float>;


    /////////////////////////////////////

    int mStps;

    std::vector<Int_t> badList;
    std::vector<Int_t> runList;

    Int_t findCentrality(int mult);
    Int_t CheckrunNumber(int runnumber);
    bool  readRunList();
    bool  readBadList();
    bool  removeBadID(int runnumber);
    Int_t openFile();

    void  DeclareHistograms();
    void  WriteHistograms();
    void cleanup();

    bool isGoodObs(double obs);

    int ReadEPDWeight(char *InFileName);
    void init_tree();

    void set_debug();
    void fill_tree();

    ClassDef(StKFParticleAnalysisMaker, 1)
};

#endif


