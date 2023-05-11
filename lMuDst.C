#include "Riostream.h"
//#include "./StRoot/particle.h"
//#include "./StRoot/particle_all.h"
//#include "./StRoot/particle_dau.h"
//#include "./StRoot/event.h"
//
//#include "./StRoot/gv.h"

void lMuDst(Int_t opt = -2, const Char_t *input = "", const Char_t *copt = "StEvent,RMuDst,mysql,tpcDb,magF,nodefault,CorrX",const Char_t *tfile =  0) {
#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#else
//    gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/KFParticle4Lambda/StRoot/event.so");
//    gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/KFParticle4Lambda/StRoot/gv.so");
//    gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/KFParticle4Lambda/StRoot/particle.so");
//    gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/KFParticle4Lambda/StRoot/particle_all.so");
//    gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/KFParticle4Lambda/StRoot/particle_dau.so");

  gSystem->AddIncludePath(" -I$ROOTROOT/root/tmva/test");
  gSystem->Load("libTMVA");
    //gSystem->Load("StMyAna");
  gSystem->Load("StRefMultCorr");
  gROOT->LoadMacro("bfc.C");
  TString Chain(copt);
  if (TString(gSystem->Getenv("STAR_VERSION")) == ".DEV2") Chain += ",TMVARank";
  bfc(opt,Chain,input,0,tfile);
#endif
}
