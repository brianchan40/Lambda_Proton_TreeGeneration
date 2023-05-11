#include <TSystem>
//#include "RooFit.h"

#include <iostream>

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StEpdGeom;
class StEpdEpInfo;
class StEpdEpFinder;

//#include "./StRoot/StMyAna/particle.h"
//#include "./StRoot/StMyAna/particle_all.h"
//#include "./StRoot/StMyAna/particle_dau.h"
//#include "./StRoot/StMyAna/event.h"
//
// #include "./StRoot/StKFParticleAnalysisMaker/gv.h"


//StChain *chain;

//#ifdef __MAKECINT__
//#pragma link C++ class vector<particle>+;
//#pragma link C++ class vector<particle_all>+;
//#pragma link C++ class vector<particle_dau>+;
//#endif

//Global Variables
TFile *outputfile;
TFile *treefile;

//void init_output(TString JobIdName);

void readPicoDst(const Char_t *inputFile="test.list", const TString JobIdName = "1234", int run=11, float energy=200.0, Char_t *ListDir="./datalist/", const int cen_ana = 0, const TString lamtype_ana = "lam", const int opt_weight = 1, const int sys_err_opt = 0)
{
	// int iJob = jobindex;
	// cout << "current job id: " << iJob << endl;

	// Char_t *InputFileList = "./new.list";
 //    std::ofstream fout(InputFileList); //create a file to write
 //    std::ifstream fin(inputFile);

 //    string line;
 //    int file_count = 0;
 //    while(getline(fin, line)) //loop wiill run till end of file
 //    {
 //        //cout << line << endl;
 //        std::istringstream iss(line);
 //        std::string token;
 //        while(std::getline(iss, token, ' '))   // but we can specify a different one
 //        {
 //            if(file_count % 2 == 0)
 //            {
 //                fout << token << "\n";     //writing data to file
 //                cout << token << endl;
 //            }
 //            file_count++;
 //        }
 //    }

 //    fin.close();
 //    fout.close();

	TString StrOutName = TString::Format("sched%scen%d_output.root", JobIdName.Data(), cen_ana);
	const Char_t *outputFile = StrOutName.Data();

	Long64_t nEvents = 1000000000;
	//nEvents = 50000;

	//Load all the System libraries
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

 	gROOT->LoadMacro("lMuDst.C");
	TString outputKFParticleQA = TString::Format("sched%scen%d_KFParticleQA.root", JobIdName.Data(), cen_ana);
	lMuDst(-1,inputFile,"ry2016,RpicoDst,kfpInter,mysql,quiet,nodefault",outputKFParticleQA.Data());

	//gSystem->AddIncludePath("-lRooFitCore");
	//gSystem->Load("libRooFit");
	gSystem->Load("StUtilities");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("StPicoEvent");
	gSystem->Load("StPicoDstMaker");
	//gSystem->Load("StEpdUtil");
	gSystem->Load("StEpdUtil.so");
	//gSystem->Load("StLambdaDecayPair");
	// gSystem->Load("StKFParticleAnalysisMaker");
	//gSystem->Load("KFParticle");

	StKFParticleAnalysisMaker *anaMaker = new StKFParticleAnalysisMaker("ana", outputFile);
	anaMaker->setRunEnergyAndListDir(run,energy,ListDir);       
	anaMaker->set_cen_curr(cen_ana);
//	anaMaker->set_lambdatype_curr(lamtype_ana);
	anaMaker->set_JobIDName_curr(JobIdName);
	anaMaker->set_opt_weight_curr(opt_weight);
	anaMaker->set_sys_err_opt_curr(sys_err_opt);

	// Loop over the links in the chain
	Int_t iInit = chain -> Init() ;
	if (iInit) chain->Fatal(iInit,"on init");
	cout<<"chain->Init();"<<endl;

	// KFParticle
	StKFParticleInterface::instance()->CleanLowPVTrackEvents();
	//StKFParticleInterface::instance()->UseHFTTracksOnly();
	StKFParticleInterface::instance()->SetSoftKaonPIDMode();
	StKFParticleInterface::instance()->SetSoftTofPidMode();
	StKFParticleInterface::instance()->SetChiPrimaryCut(10);
	//Add decays to the reconstruction list
	StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122);
	StKFParticleInterface::instance()->AddDecayToReconstructionList(-3122);

	// StPicoDstMaker & chain
	StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
	if (! maker) return;
	maker->SetStatus("*",1);
	TChain *tree = maker->chain();
	Long64_t nentries = tree->GetEntries();
	Long64_t EventsToRun = TMath::Min(nEvents,nentries);
	cout << nentries << " events in chain " << EventsToRun << " will be read." << endl;

	unsigned int found;
    tree->SetBranchStatus("EpdHit*", 1, &found);
    cout << "StPicoDstMaker::SetStatus "<< found <<" to EpdHit" << endl; 
	TClonesArray *mEpdHits_pass = new TClonesArray("StPicoEpdHit");
	tree->SetBranchAddress("EpdHit", &mEpdHits_pass);
	anaMaker->SetEpdHitsArray(mEpdHits_pass);
    
    if (nentries <= 0) return;

	//init_output(JobIdName); // Create output file
    // init_tree();

	time_t time_start;
	time_t time_now;
	time(&time_start);

	Int_t istat = 0, i = 1;
	while(i<=EventsToRun && istat!=2) {
		//if(i%5000==0){cout << endl; cout << "== Event " << i << " start ==" << endl;}
		chain->Clear();
		istat = chain->Make(i);

		if(i%100==0) {
			cout << endl; 
			cout << "== Event " << i << " finish == " << flush;
			time(&time_now);
			int time_diff = (int)difftime(time_now, time_start);
			cout << time_diff/60 << "min " << time_diff%60 << "s: " << 1.0*time_diff/i << "s/event" << endl;
		}
		if (istat == 2)
			cout << "Last  event processed. Status = " << istat << endl;
		if (istat == 3)
			cout << "Error event processed. Status = " << istat << endl;
		i++;
	}

	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << EventsToRun << endl;
	cout << "****************************************** " << endl;

	delete chain;


}

//creating output files
/*void init_output(TString JobIdName)
{
    // Create output file for histograms
    TString Name = "sched";
    Name.Append(JobIdName);
    Name.Append("_plam.root") ;
    outputfile = new TFile(Name, "recreate");

    std::cout << Name << std::endl;

    // Create output file for tree;
    TString Name4 = "sched";
    Name4.Append(JobIdName);
    Name4.Append("_tree.root");
    treefile = new TFile(Name4, "recreate");
}*/
