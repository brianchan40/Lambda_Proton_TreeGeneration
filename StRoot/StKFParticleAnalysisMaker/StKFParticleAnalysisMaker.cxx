#include "StKFParticleAnalysisMaker.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "StMessMgr.h"
#include <algorithm>
#include <TMath.h>
#include <SystemOfUnits.h>

#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"

#include "StTrackHelix.h"
#include "StLambdaDecayPair.h"
#include "MyToolkit.h"
// #include "Gamma_112_module.C"

#include "particle.h"
#include "particle_all.h"
#include "particle_dau.h"
#include "event.h"
#include "gv.h"
#include "./namespaces/gv_gamma.h"

using namespace std;

/// C++ headers

#include "stdio.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include <string>
#include <cstring>
#include <TLorentzVector.h>
#include <deque>

// #include "Gamma_112_module.C"

#define pi TMath::Pi()
#define LambdaPdgMass 1.11568
#define ProtonPdgMass 0.938272
#define PionPdgMass 0.139570
#define LambdaPdg 3122
#define ProtonPdg 2212
#define PionPdg -211

//-----------------------------------------------------------------------------
ClassImp(StKFParticleAnalysisMaker)

    //-----------------------------------------------------------------------------
    StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name, const char *outName)
    : StMaker(name)
{
    mOutName = outName;

    mRun = 0;
    mEnergy = 0;
    mListDir = "./";

    mcen_curr = 0;
    //    mlambdatype_curr = "";
    mJobIDName_curr = "";
    mopt_weight_curr = 0;
    msys_err_opt_curr = 0;
}

//-----------------------------------------------------------------------------
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker()
{
    SafeDelete(KFParticleInterface);
    SafeDelete(KFParticlePerformanceInterface);
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::openFile()
{
    // ======= StRefMultCorr ======= //
    // set up StRefMultCorr for centrality
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
    // ======= StRefMultCorr end ======= //

    cout << "-----------------------------------" << endl;
    cout << "------- user's files loaded -------" << endl;
    cout << "-----------------------------------" << endl;

    return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::Init()
{
    set_debug();

    cout << mOutName << endl;
    const char *sDir = mOutName.Data();
    int lDir = mOutName.Length();
    int iDir = lDir - 6;
    while (!std::isdigit(sDir[iDir]))
        iDir--;
    while (std::isdigit(sDir[iDir]))
        iDir--;
    mJob = std::atoi(&sDir[iDir + 1]);
    cout << "current job id: " << mJob << endl;

    PI = M_PI;
    twoPI = 2 * M_PI;

    badList.clear();
    runList.clear();

    if (!readRunList())
        return kStFatal;
    if (!readBadList())
        return kStFatal;

    DeclareHistograms();
    Int_t openFileStatus = openFile();
    if (openFileStatus == kStFatal)
        return kStFatal;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    init_tree();

    float lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
    float cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
    float par1[9] = {474.651, 474.651, 474.651, 474.651, 474.651, 3.27243e+02, 1.72351, 1.72351, 1.72351};
    float par2[9] = {3.55515, 3.55515, 3.55515, 3.55515, 3.55515, 3.56938, -7.69075, -7.69075, -7.69075};
    float par3[9] = {1.80162, 1.80162, 1.80162, 1.80162, 1.80162, 1.67113, 4.72770, 4.72770, 4.72770};

    for (int iy = 1; iy <= 9; iy++)
    {
        for (int ix = 1; ix < 101; ix++)
        {
            double eta = wt->GetXaxis()->GetBinCenter(ix);
            // cout << "wt = " << ((fabs(eta) > 3.8) ? lin[iy - 1] * eta + cub[iy - 1] * pow(eta, 3) : 0) << endl;
            // cout << "wt2 = " << ((fabs(eta) < 3.4) ? sqrt(1 - 1 / par1[iy - 1] / par1[iy - 1] / cosh(eta) / cosh(eta)) / (1 + exp((abs(eta) - par2[iy - 1]) / par3[iy - 1])) : 0 ) << endl;
            wt->SetBinContent(ix, iy, (fabs(eta) > 3.8) ? lin[iy - 1] * eta + cub[iy - 1] * pow(eta, 3) : 0);
            wt2->SetBinContent(ix, iy, (fabs(eta) < 3.4) ? sqrt(1 - 1 / par1[iy - 1] / par1[iy - 1] / cosh(eta) / cosh(eta)) / (1 + exp((abs(eta) - par2[iy - 1]) / par3[iy - 1])) : 0);
            //                              wt2->SetBinContent(ix,iy, (fabs(eta)<3.4)? 1:0);
        }
    }

    if(debug_1) cout << "before EPD!" << endl;

    char EPD_fname_old[200], EPD_fname_new[200];

    for (int ct = 0; ct < 9; ct++){
        sprintf(EPD_fname_new, "cen%d.weight_1%d%d_nop_EPD_new.root", ct, 1, 2);
        sprintf(EPD_fname_old, "cen%d.weight_1%d%d_nop_EPD.root", ct, 1, 2);
        if(debug_1) cout << "EPD_fname_old = " << EPD_fname_old << endl;
        EPD_Weight_Read = ReadEPDWeight(EPD_fname_old);

        if(debug_1) cout << "EPD_Weight_Read = " << EPD_Weight_Read << endl;

        // Prepare EPD
        StEpdEpFinder* mEpFinder_tmp = new StEpdEpFinder(9, EPD_fname_new, EPD_fname_old);
        mEpFinder_tmp->SetnMipThreshold(0.3); // recommended by EPD group
        mEpFinder_tmp->SetMaxTileWeight(2.0); // recommended by EPD group
        mEpFinder_tmp->SetEpdHitFormat(2);    // 2=pico
        mEpFinder_tmp->SetEtaWeights(1, wt);  // eta weight for 1st-order EP
        mEpFinder_tmp->SetEtaWeights(2, wt2); // eta weight for 2nd-order EP, select different eta range

        mEpFinder.push_back(mEpFinder_tmp);
        // delete mEpFinder_tmp;
        // mEpFinder_tmp = NULL;
    }

    TFile *f = GetTFile(); // These two lines need to be HERE (though I don't know /why/)- don't throw in another function
    if (f)
    {
        f->cd();
        BookVertexPlots();
    }

    // init_Gamma112(mcen_curr, mJobIDName_curr);

    return kStOK;
}

void StKFParticleAnalysisMaker::init_tree()
{
    int split_level = 1;
    int buffer_size = 5000000;
    corr_tree = new TTree("corr_tree", "Tree for Correlations");
    gROOT->ProcessLine("#include <vector>");
    // Event Variables
    // ignored: dt_nLambda, dt_nLambdaRot
    corr_tree->Branch("dt_cent", &dt_cent, "dt_cent/I");
    corr_tree->Branch("dt_num_trk", &dt_num_trk, "dt_num_trk/I");
    corr_tree->Branch("n_lam", &n_lam, "n_lam/I");
    corr_tree->Branch("n_dau", &n_dau, "n_dau/I");
    corr_tree->Branch("dt_Run", &dt_Run, "dt_Run/I");
    corr_tree->Branch("dt_TOFMult", &dt_TOFMult, "dt_TOFMult/I");
    corr_tree->Branch("dt_RefMult", &dt_RefMult, "dt_RefMult/I");
    corr_tree->Branch("dt_n_Proton", &dt_n_Proton, "dt_n_Proton/I");
    corr_tree->Branch("dt_n_Pion", &dt_n_Pion, "dt_n_Pion/I");
    corr_tree->Branch("dt_EventID", &dt_EventID, "dt_EventID/I");

    corr_tree->Branch("dt_VPDvz", &dt_VPDvz, "dt_VPDvz/F");
    corr_tree->Branch("dt_PVtxz", &dt_PVtxz, "dt_PVtxz/F");
    corr_tree->Branch("dt_PVtxx", &dt_PVtxx, "dt_PVtxx/F");
    corr_tree->Branch("dt_PVtxy", &dt_PVtxy, "dt_PVtxy/F");
    corr_tree->Branch("dt_Eweight", &dt_Eweight, "dt_Eweight/F");
    corr_tree->Branch("dt_EPD_EP1_east", &dt_EPD_EP1_east, "dt_EPD_EP1_east/F");
    corr_tree->Branch("dt_EPD_EP1_west", &dt_EPD_EP1_west, "dt_EPD_EP1_west/F");
    corr_tree->Branch("dt_EPD_EP_east", &dt_EPD_EP_east, "dt_EPD_EP_east/F");
    corr_tree->Branch("dt_EPD_EP_west", &dt_EPD_EP_west, "dt_EPD_EP_west/F");
    corr_tree->Branch("dt_Magn", &dt_Magn, "dt_Magn/F");

    // Lambda Particle
    corr_tree->Branch("p_lambda_px", &lambda_px);
    corr_tree->Branch("p_lambda_py", &lambda_py);
    corr_tree->Branch("p_lambda_pz", &lambda_pz);
    corr_tree->Branch("p_lambda_Charge", &lambda_Charge);
    corr_tree->Branch("p_lambda_dcaglobal", &lambda_dcaglobal);
    corr_tree->Branch("p_lambda_nsigma", &lambda_nsigma);
    corr_tree->Branch("p_lambda_mass", &lambda_mass);
    corr_tree->Branch("p_lambda_trk_id", &lambda_trk_id);
    corr_tree->Branch("p_lambda_hits_ratio", &lambda_hits_ratio);
    corr_tree->Branch("p_lambda_nhitsfit", &lambda_nhitsfit);
    corr_tree->Branch("p_lambda_nhitsmax", &lambda_nhitsmax);

    // // Proton Particle
    // corr_tree->Branch("p_proton_px", &proton_px);
    // corr_tree->Branch("p_proton_py", &proton_py);
    // corr_tree->Branch("p_proton_pz", &proton_pz);
    // corr_tree->Branch("p_proton_Charge", &proton_Charge);
    // corr_tree->Branch("p_proton_dcaglobal", &proton_dcaglobal);
    // corr_tree->Branch("p_proton_nsigma", &proton_nsigma);
    // corr_tree->Branch("p_proton_mass", &proton_mass);
    // corr_tree->Branch("p_proton_trk_id", &proton_trk_id);
    // corr_tree->Branch("p_proton_hits_ratio", &proton_hits_ratio);
    // corr_tree->Branch("p_proton_nhitsfit", &proton_nhitsfit);
    // corr_tree->Branch("p_proton_nhitsmax", &proton_nhitsmax);

    // // Pion Particle
    // corr_tree->Branch("p_pion_px", &pion_px);
    // corr_tree->Branch("p_pion_py", &pion_py);
    // corr_tree->Branch("p_pion_pz", &pion_pz);
    // corr_tree->Branch("p_pion_Charge", &pion_Charge);
    // corr_tree->Branch("p_pion_dcaglobal", &pion_dcaglobal);
    // corr_tree->Branch("p_pion_nsigma", &pion_nsigma);
    // corr_tree->Branch("p_pion_mass", &pion_mass);
    // corr_tree->Branch("p_pion_trk_id", &pion_trk_id);
    // corr_tree->Branch("p_pion_hits_ratio", &pion_hits_ratio);
    // corr_tree->Branch("p_pion_nhitsfit", &pion_nhitsfit);
    // corr_tree->Branch("p_pion_nhitsmax", &pion_nhitsmax);

    // All Particles
    corr_tree->Branch("p_all_px", &all_px);
    corr_tree->Branch("p_all_py", &all_py);
    corr_tree->Branch("p_all_pz", &all_pz);
    corr_tree->Branch("p_all_Charge", &all_Charge);
    corr_tree->Branch("p_all_dcaglobal", &all_dcaglobal);
    corr_tree->Branch("p_all_nSigmaProton", &all_nSigmaProton);
    corr_tree->Branch("p_all_nSigmaPion", &all_nSigmaPion);
    corr_tree->Branch("p_all_isTofTrack", &all_isTofTrack);
    corr_tree->Branch("p_all_trk_id", &all_trk_id);
    corr_tree->Branch("p_all_is_pion", &all_is_pion);
    corr_tree->Branch("p_all_is_proton", &all_is_pion);
    corr_tree->Branch("p_all_is_all", &all_is_pion);
    corr_tree->Branch("p_all_nhitsmax", &all_nhitsmax);
    corr_tree->Branch("p_all_nhitsfit", &all_nhitsfit);

    // Daughter Particles
    corr_tree->Branch("p_dau_px", &dau_px);
    corr_tree->Branch("p_dau_py", &dau_py);
    corr_tree->Branch("p_dau_pz", &dau_pz);
    corr_tree->Branch("p_dau_dcaglobal", &dau_dcaglobal);
    corr_tree->Branch("p_dau_nSigma", &dau_nSigma);
    corr_tree->Branch("p_dau_trk_id", &dau_trk_id);
    corr_tree->Branch("p_dau_nHitsFit", &dau_nHitsFit);
    corr_tree->Branch("p_dau_nHitsMax", &dau_nHitsMax);

    
    // corr_tree->Branch("p_lambda", &p_lambda_tmp);
    // corr_tree->Branch("p_proton", &p_proton_tmp);
    // corr_tree->Branch("p_pion", &p_pion_tmp);
    // corr_tree->Branch("p_all", &p_all_tmp);
    // corr_tree->Branch("p_dau", &p_dau_tmp);
}

int StKFParticleAnalysisMaker::ReadEPDWeight(char *InFileName)
{
    TFile *fWgt_EPD = new TFile(InFileName, "READ");
    if (!fWgt_EPD->IsOpen())
        return 0;
    if (fWgt_EPD->IsOpen())
    {
        if(debug_1) cout << "OPENED EPD Weight FILES" << endl;
        Read_EPD_EP1_east = (TProfile2D *)fWgt_EPD->Get("pEPD_EP1_east");
        Read_EPD_EP1_west = (TProfile2D *)fWgt_EPD->Get("pEPD_EP1_west");
        Read_EPD_EP_east = (TProfile2D *)fWgt_EPD->Get("pEPD_EP_east");
        Read_EPD_EP_west = (TProfile2D *)fWgt_EPD->Get("pEPD_EP_west");

        return 1;
    }
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::Finish()
{
    if(debug_1) cout << "in finish()" << endl;
    if(debug_1) cout << "name = " << mOutName.Data() << endl;
        
    
    if (mOutName != "")
    {
        TFile *fout = new TFile(mOutName.Data(), "RECREATE");
        if(debug_1) cout << "file created" << endl;
        fout->cd();
        if(debug_1) cout << "before writing" << endl;
        WriteHistograms();
        if(debug_1) cout << "after writing" << endl;
        fout->Close();
    }

    if(debug_1) cout << "finish_Gamma112" << endl;

    tree_file = new TFile(TString::Format("sched%s_tree.root", mJobIDName_curr.Data()).Data(), "recreate");
    tree_file->cd();
    cout << "tree written" << endl;
    corr_tree->Write();
    tree_file->Close();

    cout << "file closed" << endl;

    // finish_Gamma112(mcen_curr, mopt_weight_curr, mJobIDName_curr);

    // mEpdHits_maker->Clear();
    // delete mEpFinder;
    // mEpFinder = NULL;

    for (int ct = 0; ct < 9; ct++){
        mEpFinder.at(ct)->Finish();
    }

    return kStOK;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::cleanup()
{
    nLambda = 0;
    nLambdaRot = 0;
    nproton = 0;
    npion = 0;

    gv::details->nLambda = 0;
    gv::details->nLambdaRot = 0;
    gv::details->num_trk = 0;
    gv::details->n_Proton = 0;
    gv::details->n_Pion = 0;

    p_dau_tmp->clear();
    // p_dau_rot_tmp->clear();
    p_lambda_tmp->clear();
    // p_lambda_rot_tmp->clear();
    p_all_tmp->clear();
    p_proton_tmp->clear();
    p_pion_tmp->clear();

    lambda_px->clear();
    lambda_py->clear();
    lambda_pz->clear();
    lambda_Charge->clear();
    lambda_dcaglobal->clear();
    lambda_nsigma->clear();
    lambda_mass->clear();
    lambda_trk_id->clear();
    lambda_hits_ratio->clear();
    lambda_nhitsfit->clear();
    lambda_nhitsmax->clear();
    
    // proton_px->clear();
    // proton_py->clear();
    // proton_pz->clear();
    // proton_Charge->clear();
    // proton_dcaglobal->clear();
    // proton_nsigma->clear();
    // proton_mass->clear();
    // proton_trk_id->clear();
    // proton_hits_ratio->clear();
    // proton_nhitsfit->clear();
    // proton_nhitsmax->clear();

    // pion_px->clear();
    // pion_py->clear();
    // pion_pz->clear();
    // pion_Charge->clear();
    // pion_dcaglobal->clear();
    // pion_nsigma->clear();
    // pion_mass->clear();
    // pion_trk_id->clear();
    // pion_hits_ratio->clear();
    // pion_nhitsfit->clear();
    // pion_nhitsmax->clear();

    all_px->clear();
    all_py->clear();
    all_pz->clear();
    all_Charge->clear();
    all_dcaglobal->clear();
    all_nSigmaProton->clear();
    all_nSigmaPion->clear();
    all_trk_id->clear();
    all_isTofTrack->clear();
    all_is_pion->clear();
    all_is_proton->clear();
    all_is_all->clear();
    all_nhitsmax->clear();
    all_nhitsfit->clear();

    dau_px->clear();
    dau_py->clear();
    dau_pz->clear();
    dau_dcaglobal->clear();
    dau_nSigma->clear();
    dau_trk_id->clear();
    dau_nHitsFit->clear();
    dau_nHitsMax->clear();
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::DeclareHistograms()
{

    const int nRuns = runList.size();
    int mRunL = runList.at(0);
    int mRunH = runList.at(nRuns - 1);
    int mDayL = (int)((mRunL) / 1000 % 1000);
    int mDayH = (int)((mRunH) / 1000 % 1000) + 1;
    const int nDays = mDayH - mDayL;
    int nTmps = nDays; //>40?40:nDays;
    mStps = ceil((nDays * 1.0) / (nTmps * 1.0));
    const int nTims = ceil((nDays * 1.0) / (mStps * 1.0));

    cout << nDays << " " << mDayL << " " << mDayH << " " << nTims << " " << mStps << endl;

    hNRefMult = new TH1F("RefMult", "Reference Multiplicity", 1000, 0.0, 1000.0);
    hNRefMultA = new TH1F("RefMultA", "Reference MultiplicityA", 1000, 0.0, 1000.0);
    hNRefMultB = new TH1F("SelectRefMultM0", "Reference MultiplicityB", 1000, 0.0, 1000.0);
    hVertexXY_kf = new TH2F("VertexXY_kf", "Vertex XY Position (KF)", 200, -10.0, 10.0, 200, -10., 10);
    hVertexZ_kf = new TH1F("VertexZ_kf", "Event Vertex Z Position (KF)", 200, -100.0, 100.0);
    hVertex2D = new TH2F("Vertex2D", "VertexZ vs VPD Vz", 200, -100.0, 100.0, 200, -100., 100);
    hDiffVz = new TH1F("VertexZdiff", "VertexZ-VPDVz diff", 100, -10.0, 10.0);
    hcent = new TH1F("centrality", "centrality", nCent, 0., nCent);
    hcentw = new TH1F("centralityw", "centralityw", nCent, 0., nCent);
    hlamdist = new TH1D("lam_dist", "lam_dist", 200, 1.115684 - 0.07, 1.115684 + 0.07);
    hantilamdist = new TH1D("antilam_dist", "antilam_dist", 200, 1.115684 - 0.07, 1.115684 + 0.07);

    hdau_p_azim_dist = new TH1D("dau_p_azim_dist", "dau_p_azim_dist", 800, -4, 4);

    hcentRefM = new TProfile("hcentRefM", "hcentRefM", nCent, 0., nCent, 0, 1000);
    hcentRefW = new TProfile("hcentRefW", "hcentRefW", nCent, 0., nCent, 0, 1000);

    Ref_TOF_before = new TH2D("Ref_TOF_before", "Ref_TOF_before", 500, 0.5, 500.5, 5000, 0.5, 5000.5);
    Ref_TOF_after = new TH2D("Ref_TOF_after", "Ref_TOF_after", 500, 0.5, 500.5, 5000, 0.5, 5000.5);

    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 17; j++)
        {
            V0Mass_cent[i][j] = new TH1F(TString::Format("V0Mass_%d_%d", i, j), TString::Format("V0 Inv. Mass for Centrality-%d Pt Bin - %d", i, j), 200, 1.115684 - 0.07, 1.115684 + 0.07);

            V0Mass_anti_cent[i][j] = new TH1F(TString::Format("V0Mass_anti_%d_%d", i, j), TString::Format("V0 Inv. Mass (AntiLam) for Centrality-%d Pt Bin - %d", i, j), 200, 1.115684 - 0.07, 1.115684 + 0.07);
        }
    }

    cout << "----------------------------------" << endl;
    cout << "------- histograms claimed -------" << endl;
    cout << "----------------------------------" << endl;

    return;
}

void StKFParticleAnalysisMaker::set_debug()
{
    debug_1 = false;
    debug_2 = false;
    debug_3 = false;
    debug_4 = false;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::WriteHistograms()
{

    ///////////////////
    hNRefMult->Write();
    hNRefMultA->Write();
    hNRefMultB->Write();
    hVertexXY_kf->Write();
    hVertexZ_kf->Write();
    hVertex2D->Write();
    hDiffVz->Write();
    hcent->Write();
    hcentw->Write();
    hlamdist->Write();
    hantilamdist->Write();
    hcentRefM->Write();
    hcentRefW->Write();
    hdau_p_azim_dist->Write();

    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 17; j++)
        {
            V0Mass_cent[i][j]->Write();

            V0Mass_anti_cent[i][j]->Write();
        }
    }

    pTDist_b4_TOF->Write();
    pTDist_after_TOF->Write();
    m2_proton_vs_pT->Write();

    Ref_TOF_before->Write();
    Ref_TOF_after->Write();

    return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::setRunEnergyAndListDir(int run, double energy, char ListDir[256])
{
    mRun = run;
    mEnergy = energy;
    mListDir = ListDir;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::set_cen_curr(int c)
{
    mcen_curr = c;
}

//-----------------------------------------------------------------------------
/*void StKFParticleAnalysisMaker::set_lambdatype_curr(TString lam_type)
{
    mlambdatype_curr = lam_type;
}*/

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::set_JobIDName_curr(TString Job_ID)
{
    mJobIDName_curr = Job_ID;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::set_opt_weight_curr(int opt_weight)
{
    mopt_weight_curr = opt_weight;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::set_sys_err_opt_curr(int sys_err_opt)
{
    msys_err_opt_curr = sys_err_opt;
}

void StKFParticleAnalysisMaker::SetEpdHitsArray(TClonesArray *mEpdHitstmp)
{
    mEpdHits_maker = mEpdHitstmp;
}

//-----------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::readBadList()
{
    if (0)
        return kTRUE;
    TString inf = mListDir + "/badList/";
    inf += Form("badrun%dList%.1f.list", mRun, mEnergy);
    ifstream inrun;
    inrun.open(inf);
    if (inrun.fail())
    {
        cout << "cannot open " << inf.Data() << endl;
        return kFALSE;
    }
    Int_t runid;
    while (inrun >> runid)
    {
        badList.push_back(runid);
    }
    inrun.close();
    sort(badList.begin(), badList.end());

    vector<int>::iterator it;
    it = std::unique(badList.begin(), badList.end());
    badList.resize(std::distance(badList.begin(), it));

    cout << "badrun list :" << inf.Data() << " loaded." << endl;
    cout << "Total       :" << badList.size() << " bad runs. " << endl;
    return kTRUE;
}

//------------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::readRunList()
{
    if (0)
        return kTRUE;
    TString inf = mListDir + "/runList/";
    inf += Form("run%dList%.1f.list", mRun, mEnergy);
    ifstream inrun;
    inrun.open(inf);
    if (inrun.fail())
    {
        cout << "cannot open " << inf.Data() << endl;
        return kFALSE;
    }
    Int_t runid;
    while (inrun >> runid)
    {
        runList.push_back(runid);
    }
    inrun.close();
    sort(runList.begin(), runList.end());

    vector<int>::iterator it;
    it = std::unique(runList.begin(), runList.end());
    runList.resize(std::distance(runList.begin(), it));

    cout << "Run list :" << inf.Data() << " loaded. " << endl;
    cout << "Total    :" << runList.size() << " runs. " << endl;

    if (runList.size() < 1)
    {
        cout << "no run number found!!!" << endl;
        return kFALSE;
    }

    return kTRUE;
}

//-----------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::removeBadID(int runnumber)
{
    for (std::vector<int>::iterator it = badList.begin(); it != badList.end(); ++it)
    {
        if (runnumber == *it)
        {
            return kTRUE;
        }
    }
    return kFALSE;
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::CheckrunNumber(int runnumber)
{
    int pointer = -999;
    int id = 0;
    for (std::vector<int>::iterator it = runList.begin(); it != runList.end(); ++it)
    {
        if (runnumber == *it)
            pointer = id;
        id++;
    }

    if (pointer == -999)
        cout << "Run number are not found! " << runnumber << endl;
    return pointer;
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::findCentrality(int mult)
{
    // NOT in use
    int centrality = 0;
    float centFull[9] = {30, 46, 65, 89, 118, 154, 200, 262, 305};
    if (mult > centFull[8])
        centrality = 9;
    else if (mult > centFull[7])
        centrality = 8;
    else if (mult > centFull[6])
        centrality = 7;
    else if (mult > centFull[5])
        centrality = 6;
    else if (mult > centFull[4])
        centrality = 5;
    else if (mult > centFull[3])
        centrality = 4;
    else if (mult > centFull[2])
        centrality = 3;
    else if (mult > centFull[1])
        centrality = 2;
    else if (mult > centFull[0])
        centrality = 1;

    return centrality;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::Clear(Option_t *opt)
{
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::Make()
{
    if(debug_1) cout << "in make!" << endl;
    cleanup();

    if(debug_1) cout << "cleaned up" << endl;

    PicoDst = StPicoDst::instance();
    StPicoDst *mPicoDst = PicoDst;
    if (!mPicoDst)
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStOK;
    }

    if(debug_1) cout << "got pico" << endl;

    //     pass event
    /////////////////////////////////////////////////////////
    mEvent = (StPicoEvent *)mPicoDst->event();
    if (!mEvent)
        return kStOK;

    const int runID = mEvent->runId();
    const int evtID = mEvent->eventId();
    const int refMult = mEvent->refMult();
    const int grefMult = mEvent->grefMult();
    const int ranking = mEvent->ranking();
    const int tofMatch = mEvent->nBTOFMatch();

    const double magnet = mEvent->bField();

    if ((!mEvent->isTrigger(610001)) && (!mEvent->isTrigger(610011)) && (!mEvent->isTrigger(610021)) && (!mEvent->isTrigger(610031)) && (!mEvent->isTrigger(610041)) && (!mEvent->isTrigger(610051)))
        return kStOK;

    if(debug_1) cout << "got past triggers" << endl;

    const TVector3 Vertex3D = mEvent->primaryVertex();
    const double VertexX = Vertex3D.x();
    const double VertexY = Vertex3D.y();
    const double VertexZ = Vertex3D.z();
    const double VertexR = sqrt(VertexX * VertexX + VertexY * VertexY);
    const double vpdVz = mEvent->vzVpd();

    gv::details->Run = runID;
    gv::details->PVtxz = VertexZ;
    gv::details->PVtxx = VertexX;
    gv::details->PVtxy = VertexY;
    gv::details->VPDvz = vpdVz;
    gv::details->RefMult = refMult;

    // event cut
    // if(refMult <=2 || refMult > 1000) return kStOK;
    if (removeBadID(runID))
        return kStOK;
    if (mRefMultCorr->isBadRun(runID))
        return kStOK; // reject bad run of StRefMultCorr
    if (!mRefMultCorr->passnTofMatchRefmultCut(1. * refMult, 1. * tofMatch))
        return kStOK; // reject pileup of StRefMultCorr

    hVertex2D->Fill(VertexZ, vpdVz);
    hDiffVz->Fill(VertexZ - vpdVz);

    const double DVz = VertexZ - vpdVz;

    // if(fabs(VertexZ) > 80) return kStOK;
    if (fabs(VertexZ) > 70)
        return kStOK;
    if (sqrt(pow(VertexX, 2.) + pow(VertexY, 2.)) > 2.0)
        return kStOK;
    // if(fabs(VertexZ - vpdVz) > 3.) return kStOK;   // no vpd cut in low energy?
    // if(fabs(VertexZ - vpdVz) > 4.) return kStOK;   // no vpd cut in low energy?

    Ref_TOF_before->Fill(refMult, tofMatch);

    if (refMult > (11.070446 + 1.384789 * tofMatch - 0.002402 * pow(tofMatch, 2) + 7.514934 * pow(10, -6) * pow(tofMatch, 3) - 9.148756 * pow(10, -9) * pow(tofMatch, 4)))
        return kStOK;
    if (refMult < (-12.764095 + 0.594566 * tofMatch + 0.001040 * pow(tofMatch, 2) - 2.006525 * pow(10, -6) * pow(tofMatch, 3) + 1.127352 * pow(10, -9) * pow(tofMatch, 4)))
        return kStOK;

    Ref_TOF_after->Fill(refMult, tofMatch);

    // check run number
    int runnumberPointer = -999;
    runnumberPointer = CheckrunNumber(runID);
    if (runnumberPointer == -999)
        return kStOK;

    int dayPointer = (int)((runID) / 1000 % 1000);
    int mRunL = runList.at(0);
    int mDayL = (int)((mRunL) / 1000 % 1000);
    dayPointer -= mDayL;
    int timePointer = dayPointer / mStps;

    // StRefMultCorr
    mRefMultCorr->init(runID);
    mRefMultCorr->initEvent(refMult, VertexZ, mEvent->ZDCx());
    double refmultWght = mRefMultCorr->getWeight();
    double refmultCorr = mRefMultCorr->getRefMultCorr();
    int centrality = mRefMultCorr->getCentralityBin9(); // 0 - 8  be careful !!!!!!!!
    // double mWght    = 1;
    // double mult_corr= refMult;
    // int centrality  = findCentrality(mult_corr)-1;  // 0-8
    // if( centrality < 0 || centrality >= (nCent - 1)) return kStOK;
    int cent = centrality + 1;

    if ((cent < 1) || (cent > 9)){
        cout << "cent = " << cent << endl;
        return kStOK;
    }

    // if (centrality != 5)
        // return kStOK;

    if(debug_1) cout << "Vertex3D = " << Vertex3D.x() << ", " << Vertex3D.y() << ", " << Vertex3D.z() << endl;

    // EPD EP
    StEpdEpInfo result = mEpFinder.at(centrality)->Results(mEpdHits_maker, Vertex3D, centrality);
    // StEpdEpInfo result = mEpFinder->Results(mPicoDst->picoArray(StPicoArrays::EpdHit), Vertex3D, centrality);
    gv::details->EPD_EP1_east = result.EastPhiWeightedPsi(1);
    gv::details->EPD_EP1_west = result.WestPhiWeightedPsi(1);
    gv::details->EPD_EP_east = result.EastPhiWeightedPsi(2);
    gv::details->EPD_EP_west = result.WestPhiWeightedPsi(2);
    if (debug_1)
        std::cout << "EPD_EP1_east = " << gv::details->EPD_EP1_east << endl;
    if (debug_1)
        std::cout << "EPD_EP1_west = " << gv::details->EPD_EP1_west << endl;
    if (debug_1)
        std::cout << "EPD_EP_east = " << gv::details->EPD_EP_east << endl;
    if (debug_1)
        std::cout << "EPD_EP_west = " << gv::details->EPD_EP_west << endl;
    if ((gv::details->EPD_EP_east == gv::details->EPD_EP_west || (gv::details->EPD_EP1_east > 0.0264 && gv::details->EPD_EP1_east < 0.0265) || (gv::details->EPD_EP1_west > 0.0264 && gv::details->EPD_EP1_west < 0.0265)))
        return kStOK;

    //    gv::details->Eweight = refmultWght;
    //    gv::details->cent = centrality;

    double mWght = refmultWght;
    double mult_corr = refmultCorr;

    gv::details->TOFMult = tofMatch;
    gv::details->cent = centrality;
    gv::details->Eweight = refmultWght;

    ///////////////////////////
    hNRefMult->Fill(grefMult);
    hNRefMultB->Fill(mult_corr, mWght);
    hVertexXY_kf->Fill(VertexX, VertexY);
    hVertexZ_kf->Fill(VertexZ);
    hcent->Fill(cent);         // 1-9
    hcentw->Fill(cent, mWght); // 1-9

    ///////////////
    hcentRefM->Fill(cent, mult_corr);
    hcentRefW->Fill(cent, mult_corr, mWght);
    hcentRefM->Fill(0., mult_corr);
    hcentRefW->Fill(0., mult_corr, mWght);

    //    if(debug_1) cout << "centrality = " << centrality << endl;
    //    if(debug_1) cout << "mcen_curr = " << mcen_curr << endl;

    // if (centrality != mcen_curr)
    // return kStOK;
    ///////////////
    
    // ======= KFParticle ======= //
    vector<StLambdaDecayPair> KFParticleLambdaDecayPair;

    SetupKFParticle();
    if (InterfaceCantProcessEvent)
        return;

    hNRefMultA->Fill(mult_corr);

    if (debug_3)
        cout << "KFParticle Loop" << endl;

    for (int iKFParticle = 0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++)
    {
        //        if(debug_1) cout << "iKFParticle = " << iKFParticle << endl;
        const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
        int upQ;
        if (particle.GetPDG() == LambdaPdg)
            upQ = 1;
        else if (particle.GetPDG() == -1 * LambdaPdg)
            upQ = -1;
        else
            continue;
        int eLambda = -(upQ - 1) / 2; // 0 if Lambda, 1 if AntiLambda

        SetDaughterTrackPointers(iKFParticle);
        //        if(debug_1) cout << "ProtonTrackIndex = " << ProtonTrackIndex << endl;
        //        if(debug_1) cout << "PionTrackIndex = " << PionTrackIndex << endl;
        if (ProtonTrackIndex == -99999 || PionTrackIndex == -99999)
            continue;
        if (!ProtonTrack)
            continue;
        if (!PionTrack)
            continue;

        //        if(debug_1) cout << "passed all cuts in KFParticle Loop" << endl;

        double dmass = -999;             // just a placeholder
        TLorentzVector p4Pair, p4Proton; // just a placeholder
        StLambdaDecayPair TmpLambdaDecayPair(p4Pair, p4Proton, ProtonTrackIndex, PionTrackIndex, (eLambda == 0), dmass);
        KFParticleLambdaDecayPair.push_back(TmpLambdaDecayPair);

        p4Pair.Clear();
        p4Proton.Clear();
        particle.Clear();
    } // End loop over KFParticles
    // ======= KFParticle end ======= //

    //    if(debug_1) cout << "KFParticleLambdaDecayPair.size() = " << KFParticleLambdaDecayPair.size() << endl;

    // ======= Lambda loop ======= //
    for (int j = 0; j < KFParticleLambdaDecayPair.size(); j++)
    {
        if (debug_3)
            cout << "Lambda Loop Index = " << j << endl;
        int i = KFParticleLambdaDecayPair[j].get_idxProton();
        int k = KFParticleLambdaDecayPair[j].get_idxPion();
        if (k == i)
            continue;

        StPicoTrack *mTrackI = (StPicoTrack *)mPicoDst->track(i);

        int mchgI = mTrackI->charge();
        int mhitI = mTrackI->nHitsFit();
        double mMomI = mTrackI->gMom().Mag();
        double mp0xI = mTrackI->gMom().X();
        double mp0yI = mTrackI->gMom().Y();
        double mp0zI = mTrackI->gMom().Z();
        double mpt0I = mTrackI->gMom().Perp();
        double mphiI = mTrackI->gMom().Phi();
        double metaI = mTrackI->gMom().PseudoRapidity();
        double mdcaI = mTrackI->gDCA(Vertex3D).Mag();

        if (mphiI < 0)
            mphiI += 2 * M_PI;
        if (mphiI >= 2 * M_PI)
            mphiI -= 2 * M_PI;

        StPicoTrack *mTrackK = (StPicoTrack *)mPicoDst->track(k);

        int mchgK = mTrackK->charge();
        int mhitK = mTrackK->nHitsFit();
        double mMomK = mTrackK->gMom().Mag();
        double mp0xK = mTrackK->gMom().X();
        double mp0yK = mTrackK->gMom().Y();
        double mp0zK = mTrackK->gMom().Z();
        double mpt0K = mTrackK->gMom().Perp();
        double mphiK = mTrackK->gMom().Phi();
        double metaK = mTrackK->gMom().PseudoRapidity();
        double mdcaK = mTrackK->gDCA(Vertex3D).Mag();

        if (mphiK < 0)
            mphiK += 2 * M_PI;
        if (mphiK >= 2 * M_PI)
            mphiK -= 2 * M_PI;

        bool isPP = mchgI > 0 && mchgK > 0;
        bool isNN = mchgI < 0 && mchgK < 0;
        bool isPN = mchgI > 0 && mchgK < 0;
        bool isNP = mchgI < 0 && mchgK > 0;
        bool isSelfLambda = isPN;
        bool isAntiLambda = isNP;

        // particle_dau p_dau_1(mp0xI, mp0yI, mp0zI, mTrackI->id(), mhitI, mTrackI->nSigmaProton(), mdcaI);
        // particle_dau p_dau_2(mp0xK, mp0yK, mp0zK, mTrackK->id(), mhitK, mTrackK->nSigmaPion(), mdcaK);

        // remove the SS cases
        if (isPP || isNN)
            continue;

        // reconstruction of V0, the parent particle
        TVector3 xv0, op1, op2;
        double dca1to2 = closestDistance(mTrackI, mTrackK, magnet, Vertex3D, xv0, op1, op2);
        TVector3 pv0 = op1 + op2;
        TVector3 xv0toPV = xv0 - Vertex3D;
        double rdotp = xv0toPV.Dot(pv0);
        double dcav0toPV = rdotp * rdotp / pv0.Mag2();
        dcav0toPV = sqrt(xv0toPV.Mag2() - dcav0toPV);
        double v0decaylength = xv0toPV.Mag();
        double v0cosrdotp = rdotp / v0decaylength / pv0.Mag();

        TLorentzVector p4ProtonI, p4PionK;
        p4ProtonI.SetPxPyPzE(op1.X(), op1.Y(), op1.Z(), sqrt(op1.Mag2() + pmass * pmass));
        p4PionK.SetPxPyPzE(op2.X(), op2.Y(), op2.Z(), sqrt(op2.Mag2() + pimass * pimass));

        // ProtonI & PionK
        TLorentzVector p4Pair = p4ProtonI + p4PionK;
        double massPair = p4Pair.M();
        double ptPair = p4Pair.Pt();

        if ((massPair > 1.077) && (massPair < 1.078))
        {
            hdau_p_azim_dist->Fill(p4ProtonI.Phi());
            hdau_p_azim_dist->Fill(p4PionK.Phi());
        }

        if ((ptPair < 0.4) || (ptPair > 3.0))
            continue;

        double etaPair = p4Pair.Eta();
        double phiPair = p4Pair.Phi();

        Float_t charge_tmp = 0;
        if (isSelfLambda)
            charge_tmp = 1;
        else if (isAntiLambda)
            charge_tmp = -1;

        if (isSelfLambda)
        {
            hlamdist->Fill(massPair);

            if ((int)floor((ptPair - 0.5) / 0.1) >= 17)
            {
                V0Mass_cent[centrality][16]->Fill(massPair);
            }
            else if ((ptPair >= 0.5))
            {
                V0Mass_cent[centrality][(int)floor((ptPair - 0.5) / 0.1)]->Fill(massPair);
            }
        }
        else if (isAntiLambda)
        {
            hantilamdist->Fill(massPair);

            if ((int)floor((ptPair - 0.5) / 0.1) >= 17)
            {
                V0Mass_anti_cent[centrality][16]->Fill(massPair);
            }
            else if ((ptPair >= 0.5))
            {
                V0Mass_anti_cent[centrality][(int)floor((ptPair - 0.5) / 0.1)]->Fill(massPair);
            }
        }

        if (debug_3)
            cout << "massPair = " << massPair << endl;

        // if(!((massPair < 1.113) || (massPair > 1.119)))
        // {
        if (debug_3)
            cout << "In real!" << endl;
        // particle Lambda_particle(p4Pair.Px(), p4Pair.Py(), p4Pair.Pz(), charge_tmp, dcav0toPV, massPair, nLambda, 0);

        nLambda++;

        gv::p_lambda[gv::details->nLambda / 2].px = p4Pair.Px();
        gv::p_lambda[gv::details->nLambda / 2].py = p4Pair.Py();
        gv::p_lambda[gv::details->nLambda / 2].pz = p4Pair.Pz();
        gv::p_lambda[gv::details->nLambda / 2].mass = massPair;
        gv::p_lambda[gv::details->nLambda / 2].dcaglobal = 0;
        gv::p_lambda[gv::details->nLambda / 2].Charge = charge_tmp;

        gv::p_dau[gv::details->nLambda].px = mp0xI;
        gv::p_dau[gv::details->nLambda].py = mp0yI;
        gv::p_dau[gv::details->nLambda].pz = mp0zI;
        gv::p_dau[gv::details->nLambda].nHitsFit = mhitI;
        gv::p_dau[gv::details->nLambda].nSigma = mTrackI->nSigmaProton();
        gv::p_dau[gv::details->nLambda].dcaglobal = mdcaI;
        gv::p_dau[gv::details->nLambda].nHitsMax = mTrackI->nHitsMax();
        gv::p_dau[gv::details->nLambda].trk_id = mTrackI->id();

        gv::details->nLambda++;

        gv::p_dau[gv::details->nLambda].px = mp0xK;
        gv::p_dau[gv::details->nLambda].py = mp0yK;
        gv::p_dau[gv::details->nLambda].pz = mp0zK;
        gv::p_dau[gv::details->nLambda].nHitsFit = mhitK;
        gv::p_dau[gv::details->nLambda].nSigma = mTrackK->nSigmaProton();
        gv::p_dau[gv::details->nLambda].dcaglobal = mdcaK;
        gv::p_dau[gv::details->nLambda].nHitsMax = mTrackK->nHitsMax();
        gv::p_dau[gv::details->nLambda].trk_id = mTrackK->id();

        gv::details->nLambda++;

        // p_dau_tmp->push_back(p_dau_1);
        // p_dau_tmp->push_back(p_dau_2);

        // p_lambda_tmp->push_back(Lambda_particle);
        // }
        if (!((massPair < 1.09) || (massPair > 1.14)))
        {
            if (!((massPair < 1.125) && (massPair > 1.105)))
            {
                if (debug_3)
                    cout << "In bkg!" << endl;

                // particle Lambda_particle(p4Pair.Px(), p4Pair.Py(), p4Pair.Pz(), charge_tmp, dcav0toPV, massPair, nLambdaRot, 0);

                nLambdaRot++;

                gv::p_lambda_rot[gv::details->nLambdaRot / 2].px = p4Pair.Px();
                gv::p_lambda_rot[gv::details->nLambdaRot / 2].py = p4Pair.Py();
                gv::p_lambda_rot[gv::details->nLambdaRot / 2].pz = p4Pair.Pz();
                gv::p_lambda_rot[gv::details->nLambdaRot / 2].mass = massPair;
                gv::p_lambda_rot[gv::details->nLambdaRot / 2].dcaglobal = 0;
                gv::p_lambda_rot[gv::details->nLambdaRot / 2].Charge = charge_tmp;

                gv::p_dau_rot[gv::details->nLambdaRot].px = mp0xI;
                gv::p_dau_rot[gv::details->nLambdaRot].py = mp0yI;
                gv::p_dau_rot[gv::details->nLambdaRot].pz = mp0zI;
                gv::p_dau_rot[gv::details->nLambdaRot].nHitsFit = mhitI;
                gv::p_dau_rot[gv::details->nLambdaRot].nSigma = mTrackI->nSigmaProton();
                gv::p_dau_rot[gv::details->nLambdaRot].dcaglobal = mdcaI;
                gv::p_dau_rot[gv::details->nLambdaRot].nHitsMax = mTrackI->nHitsMax();
                gv::p_dau_rot[gv::details->nLambdaRot].trk_id = mTrackI->id();

                gv::details->nLambdaRot++;

                gv::p_dau_rot[gv::details->nLambdaRot].px = mp0xK;
                gv::p_dau_rot[gv::details->nLambdaRot].py = mp0yK;
                gv::p_dau_rot[gv::details->nLambdaRot].pz = mp0zK;
                gv::p_dau_rot[gv::details->nLambdaRot].nHitsFit = mhitK;
                gv::p_dau_rot[gv::details->nLambdaRot].nSigma = mTrackK->nSigmaProton();
                gv::p_dau_rot[gv::details->nLambdaRot].dcaglobal = mdcaK;
                gv::p_dau_rot[gv::details->nLambdaRot].nHitsMax = mTrackK->nHitsMax();
                gv::p_dau_rot[gv::details->nLambdaRot].trk_id = mTrackK->id();

                gv::details->nLambdaRot++;

                // p_dau_rot_tmp->push_back(p_dau_1);
                // p_dau_rot_tmp->push_back(p_dau_2);

                // p_lambda_rot_tmp->push_back(Lambda_particle);
            }
        }

        p4ProtonI.Clear();
        p4PionK.Clear();
        xv0.Clear();
        op1.Clear();
        op2.Clear();
    }
    // ======= Lambda loop ends ======= //

    if (debug_1)
        cout << "nLambda before beginning = " << gv::details->nLambda << endl;
    if (debug_1)
        cout << "nLambdaRot before beginning = " << gv::details->nLambdaRot << endl;

    // Gamma_112_module(mcen_curr, mopt_weight_curr, msys_err_opt_curr, mJobIDName_curr, mEpdHits_maker);

    if ((gv::details->nLambda == 0) && (gv::details->nLambdaRot == 0))
        return kStOK;

    fill_tree();

    if (debug_3) cout << "tree filled!" << endl;
    corr_tree->Fill();

    /////////////////////////////////////////////////////////
    return kStOK;
}

void StKFParticleAnalysisMaker::fill_tree(){
    
    if (debug_2) cout << "n_lam = " << (int) gv::details->nLambda/2 << endl;
    if (debug_2) cout << "n_dau = " << gv::details->nLambda << endl;
    if (debug_2) cout << "gv::details->num_trk = " << gv::details->num_trk << endl;
    if (debug_2) cout << "gv::details->n_Proton = " << gv::details->n_Proton << endl;
    if (debug_2) cout << "gv::details->n_Pion = " << gv::details->n_Pion << endl;
    
    dt_cent = gv::details->cent;
    dt_num_trk = gv::details->num_trk;
    n_lam = (int) gv::details->nLambda/2;
    n_dau = gv::details->nLambda;
    dt_Run = gv::details->Run;
    dt_TOFMult = gv::details->TOFMult;
    dt_RefMult = gv::details->RefMult;
    dt_n_Proton = gv::details->n_Proton;
    dt_n_Pion = gv::details->n_Pion;
    dt_EventID = gv::details->EventID;
    dt_VPDvz = gv::details->VPDvz;
    dt_PVtxz = gv::details->PVtxz;
    dt_PVtxx = gv::details->PVtxx;
    dt_PVtxy = gv::details->PVtxy;
    dt_Eweight = gv::details->Eweight;
    dt_EPD_EP1_east = gv::details->EPD_EP1_east;
    dt_EPD_EP1_west = gv::details->EPD_EP1_west;
    dt_EPD_EP_east = gv::details->EPD_EP_east;
    dt_EPD_EP_west = gv::details->EPD_EP_west;
    dt_Magn = gv::details->Magn;
    
    for (int pl = 0; pl < n_lam; pl++)
    {
        lambda_px->push_back(gv::p_lambda[pl].px);
        lambda_py->push_back(gv::p_lambda[pl].py);
        lambda_pz->push_back(gv::p_lambda[pl].pz);
        lambda_Charge->push_back(gv::p_lambda[pl].Charge);
        lambda_dcaglobal->push_back(gv::p_lambda[pl].dcaglobal);
        lambda_nsigma->push_back(gv::p_lambda[pl].nsigma);
        lambda_mass->push_back(gv::p_lambda[pl].mass);
        lambda_trk_id->push_back(gv::p_lambda[pl].trk_id);
        lambda_hits_ratio->push_back(gv::p_lambda[pl].hits_ratio);
        lambda_nhitsfit->push_back(gv::p_lambda[pl].nhitsfit);
        lambda_nhitsmax->push_back(gv::p_lambda[pl].nhitsmax);   
    }
    // for (int pl = 0; pl < dt_n_Proton; pl++)
    // {
    //     proton_px->push_back(gv::p_proton[pl].px);
    //     proton_py->push_back(gv::p_proton[pl].py);
    //     proton_pz->push_back(gv::p_proton[pl].pz);
    //     proton_Charge->push_back(gv::p_proton[pl].Charge);
    //     proton_dcaglobal->push_back(gv::p_proton[pl].dcaglobal);
    //     proton_nsigma->push_back(gv::p_proton[pl].nsigma);
    //     proton_mass->push_back(gv::p_proton[pl].mass);
    //     proton_trk_id->push_back(gv::p_proton[pl].trk_id);
    //     proton_hits_ratio->push_back(gv::p_proton[pl].hits_ratio);
    //     proton_nhitsfit->push_back(gv::p_proton[pl].nhitsfit);
    //     proton_nhitsmax->push_back(gv::p_proton[pl].nhitsmax);   
    // }
    // for (int pl = 0; pl < dt_n_Pion; pl++)
    // {
    //     pion_px->push_back(gv::p_pion[pl].px);
    //     pion_py->push_back(gv::p_pion[pl].py);
    //     pion_pz->push_back(gv::p_pion[pl].pz);
    //     pion_Charge->push_back(gv::p_pion[pl].Charge);
    //     pion_dcaglobal->push_back(gv::p_pion[pl].dcaglobal);
    //     pion_nsigma->push_back(gv::p_pion[pl].nsigma);
    //     pion_mass->push_back(gv::p_pion[pl].mass);
    //     pion_trk_id->push_back(gv::p_pion[pl].trk_id);
    //     pion_hits_ratio->push_back(gv::p_pion[pl].hits_ratio);
    //     pion_nhitsfit->push_back(gv::p_pion[pl].nhitsfit);
    //     pion_nhitsmax->push_back(gv::p_pion[pl].nhitsmax);   
    // }
    for (int pl = 0; pl < dt_num_trk; pl++)
    {
        all_px->push_back(gv::p_all[pl].px);
        all_py->push_back(gv::p_all[pl].py);
        all_pz->push_back(gv::p_all[pl].pz);
        all_Charge->push_back(gv::p_all[pl].Charge);
        all_dcaglobal->push_back(gv::p_all[pl].dcaglobal);
        all_nSigmaProton->push_back(gv::p_all[pl].nSigmaProton);
        all_nSigmaPion->push_back(gv::p_all[pl].nSigmaPion);
        all_isTofTrack->push_back(gv::p_all[pl].isTofTrack);
        all_trk_id->push_back(gv::p_all[pl].trk_id);
        all_is_pion->push_back(gv::p_all[pl].is_pion);
        all_is_proton->push_back(gv::p_all[pl].is_proton);
        all_is_all->push_back(gv::p_all[pl].is_all);
        all_nhitsmax->push_back(gv::p_all[pl].nhitsmax);
        all_nhitsfit->push_back(gv::p_all[pl].nhitsfit);
    }
    for (int pl = 0; pl < n_dau; pl++)
    {
        dau_px->push_back(gv::p_dau[pl].px);
        dau_py->push_back(gv::p_dau[pl].py);
        dau_pz->push_back(gv::p_dau[pl].pz);
        dau_dcaglobal->push_back(gv::p_dau[pl].dcaglobal);
        dau_nSigma->push_back(gv::p_dau[pl].nSigma);
        dau_trk_id->push_back(gv::p_dau[pl].trk_id);
        dau_nHitsFit->push_back(gv::p_dau[pl].nHitsFit);
        dau_nHitsMax->push_back(gv::p_dau[pl].nHitsMax);   
    }
}

//------------------------------------------------------------
bool StKFParticleAnalysisMaker::isGoodObs(double Obs)
{

    bool mGoodObs = true;
    mGoodObs = mGoodObs && !std::isnan(Obs);
    mGoodObs = mGoodObs && !std::isinf(Obs);
    // mGoodObs = mGoodObs && fabs(Obs)<2.1;

    return mGoodObs;
}

//------------------------------------------------------------
void StKFParticleAnalysisMaker::BookVertexPlots()
{
    KFParticleInterface = new StKFParticleInterface;
    bool storeMCHistograms = false;
    KFParticlePerformanceInterface = new StKFParticlePerformanceInterface(KFParticleInterface->GetTopoReconstructor(), storeMCHistograms);
}

//------------------------------------------------------------
void StKFParticleAnalysisMaker::SetupKFParticle()
{
    if (debug_3)
        cout << "SetupKFParticle" << endl;

    int maxGBTrackIndex = -1; // find max global track index
    for (unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++)
    {
        StPicoTrack *track = PicoDst->track(iTrack);
        if (!track)
            continue;
        if (track->id() > maxGBTrackIndex)
            maxGBTrackIndex = track->id();
    }
    vector<KFMCTrack> mcTracks(0);
    vector<int> triggeredTracks;
    vector<int> mcIndices(maxGBTrackIndex + 1);
    for (unsigned int iIndex = 0; iIndex < mcIndices.size(); iIndex++)
        mcIndices[iIndex] = -1;
    if (maxGBTrackIndex > 0)
        KFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex + 1);

    if (!KFParticleInterface->ProcessEvent(PicoDst, triggeredTracks))
        InterfaceCantProcessEvent = true;
    else
        InterfaceCantProcessEvent = false;

    trackMap.resize(maxGBTrackIndex + 1, -1); // make a map from trackID to track index in global track array

    //check up to here

    for (unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++)
    {
        //        if(debug_1) cout << "in all tracks loop" << endl;

        StPicoTrack *track = PicoDst->track(iTrack);
        if (!track)
            continue;
        int index = track->id();
        trackMap[index] = iTrack;

        //        if(debug_1) cout << "after !track cut" << endl;

        TVector3 primary_p;
        StPicoPhysicalHelix helix = track->helix(mEvent->bField()); // inner helix. good for dca to PV.
        if (debug_4)
            cout << "mEvent->bField() = " << mEvent->bField() << endl;
        primary_p = helix.momentum(mEvent->bField() * kilogauss); // momentum at origin
        if (debug_4)
            cout << "mEvent->bField() * kilogauss = " << mEvent->bField() * kilogauss << endl;
        if (debug_4)
            cout << "primary_p.Pt() = " << primary_p.Pt() << endl;
        double pathlength = helix.pathLength(mEvent->primaryVertex(), false); // do scan periods. NOTE: the default is false. this is not necessary for tracks with pt>0.15GeV/c
        TVector3 dca = helix.at(pathlength) - mEvent->primaryVertex();

        int index_btof = track->bTofPidTraitsIndex();
        int tofflag = -1;
        if (index_btof >= 0)
            tofflag = (PicoDst->btofPidTraits(index_btof))->btofMatchFlag();

        
        if (!(primary_p.Pt() <= 2.0) && (primary_p.Pt() >= 0.15)) continue;
        if(((float)track->nHitsFit()) / ((float)track->nHitsMax()) < 0.52) continue; // To eliminate track splitting effects
        if((track->nHitsFit() <= 15) || (abs(track->charge()) != 1)) continue;
        if (dca.Mag() > 3) continue;

        gv::p_all[gv::details->num_trk].px = primary_p.X();
        gv::p_all[gv::details->num_trk].py = primary_p.Y();
        gv::p_all[gv::details->num_trk].pz = primary_p.Z();
        gv::p_all[gv::details->num_trk].Charge = track->charge();
        gv::p_all[gv::details->num_trk].nSigmaProton = track->nSigmaProton();
        gv::p_all[gv::details->num_trk].nSigmaPion = track->nSigmaPion();
        gv::p_all[gv::details->num_trk].trk_id = track->id();
        gv::p_all[gv::details->num_trk].dcaglobal = dca.Mag();
        if (track->isTofTrack())
            gv::p_all[gv::details->num_trk].isTofTrack = true;
        else
            gv::p_all[gv::details->num_trk].isTofTrack = false;

        // gv::p_all[gv::details->num_trk].hits_ratio = ((float)track->nHitsFit()) / ((float)track->nHitsMax());
        gv::p_all[gv::details->num_trk].nhitsfit = (float)track->nHitsFit();
        gv::p_all[gv::details->num_trk].nhitsmax = (float)track->nHitsMax();

        gv::p_all[gv::details->num_trk].is_pion = false;
        gv::p_all[gv::details->num_trk].is_all = false;
        gv::p_all[gv::details->num_trk].is_proton = false;

        if ((primary_p.Eta() <= 1.0) && (primary_p.Eta() >= -1.0))
        {
            gv::p_all[gv::details->num_trk].is_all = true;
        }

        if (!((track->nHitsDedx() < 15) || (index_btof < 0))){
            double tof = (PicoDst->btofPidTraits(index_btof))->btof();
            double BtofYLocal = (PicoDst->btofPidTraits(index_btof))->btofYLocal();
            pTDist_b4_TOF->Fill(primary_p.Perp());

            if (!((tofflag < 1) || (tof <= 0) || (BtofYLocal < -1.8) || (BtofYLocal > 1.8))){
                pTDist_after_TOF->Fill(primary_p.Perp());
                
                double beta = (PicoDst->btofPidTraits(index_btof))->btofBeta();
                if (beta != 0){
                    float mass2 = primary_p.Mag2() * (1.0 / beta / beta - 1.0);
                    TVector3 pp = track->pMom();

                    double nsigmaproton = track->nSigmaProton();
                    double nsigmapion = track->nSigmaPion();

                    TVector3 gdca = track->gDCA(mEvent->primaryVertex());

                    if (debug_4)
                        cout << "nsigmaproton = " << nsigmaproton << endl;
                    if (debug_4)
                        cout << "dca.Mag() = " << dca.Mag() << endl;
                    if (debug_4)
                        cout << "mass2 = " << mass2 << endl;

                    m2_proton_vs_pT->Fill(primary_p.Perp(), mass2);

                    if (!((nsigmaproton < -2) || (nsigmaproton > 2) || (dca.Mag() > 1)))
                    {
                        if (!((mass2 < 0.8) || ((mass2 > 1.0))))
                        {
                            if (debug_4)
                                cout << "got proton!" << endl;

                            gv::p_all[gv::details->num_trk].is_proton = true;
                            gv::details->n_Proton++;
                            nproton++;
                        }
                    }

                    if (!((pp.Perp() < 0.2) || (pp.Mag() > 1.6)))
                    {
                        if (!((nsigmapion < -2) || (nsigmapion > 2) || (dca.Mag() > 1)))
                        {
                            if (!((mass2 < -0.01) || ((mass2 > 0.1))))
                            {
                                gv::details->n_Pion++;
                                npion++;
                                gv::p_all[gv::details->num_trk].is_pion = true;
                            }
                        }
                    }
                }
            }
        }

        if(gv::p_all[gv::details->num_trk].is_pion || gv::p_all[gv::details->num_trk].is_proton || gv::p_all[gv::details->num_trk].is_all) gv::details->num_trk++;
        // cout << "gv::details->num_trk = " << gv::details->num_trk << endl;
    }

    if (debug_3)
        cout << "details->num_trk = " << gv::details->num_trk << endl;

    KFParticlePerformanceInterface->SetMCTracks(mcTracks);
    KFParticlePerformanceInterface->SetMCIndexes(mcIndices);
    KFParticlePerformanceInterface->SetCentralityBin(-1);
    KFParticlePerformanceInterface->SetCentralityWeight(1.);
    KFParticlePerformanceInterface->SetPrintEffFrequency(100000);
    KFParticlePerformanceInterface->PerformanceAnalysis();
} // void SetupKFParticle

//------------------------------------------------------------
void StKFParticleAnalysisMaker::SetDaughterTrackPointers(int iKFParticle) // Get Daughter Tracks to calculate decay variables manually
{
    const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
    int upQ;
    if (particle.GetPDG() == LambdaPdg)
        upQ = 1;
    else if (particle.GetPDG() == -1 * LambdaPdg)
        upQ = -1;
    else
        upQ = 0;
    ProtonTrackIndex = -99999, PionTrackIndex = -99999;
    for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
    {
        const int daughterId = particle.DaughterIds()[iDaughter];
        const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
        const int globalTrackId = daughter.DaughterIds()[0];
        int trackIndex = trackMap[globalTrackId];
        if (daughter.GetPDG() == upQ * ProtonPdg)
            ProtonTrackIndex = trackIndex;
        else if (daughter.GetPDG() == upQ * PionPdg)
            PionTrackIndex = trackIndex;
    } // iDaughter
    ProtonTrack = PicoDst->track(ProtonTrackIndex);
    PionTrack = PicoDst->track(PionTrackIndex);
} // void SetDaughterTrackPointers
