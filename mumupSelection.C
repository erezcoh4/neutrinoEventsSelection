/*
 selection of muon - proton scattering events
 for GBDT proton classification study
 erez cohen, august 2016
 */

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "TProfile.h"
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TVector3.h>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <stdlib.h>

using namespace std;

float FVx = 256.35;
float FVy = 233;
float FVz = 1036.8;
float borderx = 20.;
float bordery = 20.;
float borderz = 10.;
float dEdxMCScale[3] = {1.118,1.308,1.194};
bool inFV(float x, float y, float z){
    if(x< (FVx - borderx) && (x>borderx) && (y<(FVy/2. - bordery)) && (y>(-FVy/2. + bordery)) && (z<(FVz-borderz)) && (z>borderz))
    return true;
    else
    return false;
}

float scaledEdx(float x, int plane, bool isdata){
    float dEdx = 1.63;
    float p0_data[3] = {1.927, 2.215, 1.936};
    float p1_data[3] = {0.001495, 0.0001655, 0.001169};
    float p0_mc[3] = {1.843, 1.904, 1.918};
    float p1_mc[3] = {-0.0008329, -0.001357, -0.0007563};
    if (isdata){
        return dEdx/(p0_data[plane]+x*p1_data[plane]);
    }
    else{
        return  dEdx/(p0_mc[plane]+x*p1_mc[plane]);
    }
}

void NeutrinoSelection(int is = -1){
    
    bool debug = false;
    
    const int nsamples = 4;
    const int ncuts = 10;
    int nevents[nsamples][ncuts];
    for (int i = 0; i<nsamples; ++i){
        for (int j = 0; j<ncuts; ++j){
            nevents[i][j] = 0;
        }
    }
    
    //  TChain *ch[nsamples];
    string filelists[nsamples] = {
        //"/uboone/app/users/tjyang/larsoft_dev/ccinclusive/SelectionII/beamon.list",
        "/pnfs/uboone/persistent/users/aschu/devel/v05_11_01/hadd/GOODBNBBEAM/filesana.list",
        "/pnfs/uboone/persistent/users/aschu/devel/v05_11_01/hadd/GOODEXTBNB/filesana.list",
        "/uboone/app/users/tjyang/larsoft_dev/ccinclusive/tjyang/bnb_nu_cosmic.txt",
        "/uboone/app/users/tjyang/larsoft_dev/ccinclusive/SelectionII/corsikaintime.txt"
    };
    
    string samplename[nsamples] = {
        "BeamOnData",
        "BeamOffData",
        "BNBCosmicMC",
        "CorsikaInTime"};
    
    //  const string tracklabel = "pandoraNuPMA";
    //  const string vertexlabel = "pmtrack";
    
    const string tracklabel = "pandoraNu";
    const string vertexlabel = "pandoraNu";
    
    //  const string tracklabel = "pmtrack";
    //  const string vertexlabel = "pmtrack";
    
    const int ncomps = 12;
    string outhistfile = "hist_"+tracklabel+"_"+vertexlabel+".root";
    if (is!=-1){
        outhistfile = "hist_"+tracklabel+"_"+vertexlabel+to_string(is)+".root";
    }
    cout<<outhistfile<<endl;
    //exit(0);
    TFile *OutHist = new TFile(outhistfile.c_str(),"recreate");
    //Final results
    TH1F *hTrackLength[nsamples][ncomps];
    TH1F *hTrackMul[nsamples][ncomps];
    TH1F *hTrackCosz[nsamples][ncomps];
    TH1F *hTrackPhi[nsamples][ncomps];
    TH1F *hVtxx[nsamples][ncomps];
    TH1F *hVtxy[nsamples][ncomps];
    TH1F *hVtxz[nsamples][ncomps];
    TH1F *hTrackStartx[nsamples][ncomps];
    TH1F *hTrackStarty[nsamples][ncomps];
    TH1F *hTrackStartz[nsamples][ncomps];
    TH1F *hTrackEndx[nsamples][ncomps];
    TH1F *hTrackEndy[nsamples][ncomps];
    TH1F *hTrackEndz[nsamples][ncomps];
    
    //For development
    TH1F *hMul[nsamples][ncomps];          //Track multiplicity before selection
    TH1F *hTrackLength1[nsamples][ncomps]; //Track length for mult=1
    TH1F *hTrackLength2[nsamples][ncomps]; //Minimal track length for mult=3
    TH1F *hTrackLength3[nsamples][ncomps]; //Minimal track length for mult=2
    TH1F *hCos1[nsamples][ncomps]; //cosine of angle between two tracks for mult=2
    TH1F *hCos2[nsamples][ncomps]; //cosine of angle between two longer tracks for mult=3
    TH2F *hSingleTrackLengthRatioVsdEdxRatio[nsamples][ncomps]; //Mult = 2, track length y vs dEdx ratio
    TH2F *hLongTrackdEdxStartEnd[nsamples][ncomps]; //Mult = 2, startdEdx vs enddEdx
    TH1F *hTrackEndy1[nsamples][ncomps]; //Mult = 2, track end of the longer track
    TH1F *hCos3[nsamples][ncomps]; //cosine of angle btween two longest tracks for mult>=2
    TH2F *hCos0VsLen1[nsamples][ncomps]; //cosy of longest track vs length of 2nd longest track
    TH2F *hCosVsLen[nsamples][ncomps]; //cosy of track vs length of track for mult=1
    TH2F *hdEdxVsX[nsamples][ncomps][3]; // dEdx vs X
    TH2F *hdEdxVsXCor[nsamples][ncomps][3]; // dEdx vs X
    
    for (int i = 0; i<nsamples; ++i){
        for (int j = 0; j<ncomps; ++j){
            hTrackLength[i][j] = new TH1F(Form("hTrackLength_%d_%d",i,j),";Track Length (cm);N events",20,0,1000);
            hTrackLength[i][j]->Sumw2();
            hTrackMul[i][j] = new TH1F(Form("hTrackMul_%d_%d",i,j),";Track Multiplicity;N events",10,0,10);
            hTrackMul[i][j]->Sumw2();
            hVtxx[i][j] = new TH1F(Form("hVtxx_%d_%d",i,j),";Vertex X (cm);N events", 20, 0, 260);
            hVtxx[i][j]->Sumw2();
            hVtxy[i][j] = new TH1F(Form("hVtxy_%d_%d",i,j),";Vertex Y (cm);N events", 20, -120, 120);
            hVtxy[i][j]->Sumw2();
            hVtxz[i][j] = new TH1F(Form("hVtxz_%d_%d",i,j),";Vertex Z (cm);N events", 20, 0, 1000);
            hVtxz[i][j]->Sumw2();
            hTrackStartx[i][j] = new TH1F(Form("hTrackStartx_%d_%d",i,j),";Track Start X (cm);N events", 20, 0, 260);
            hTrackStartx[i][j]->Sumw2();
            hTrackStarty[i][j] = new TH1F(Form("hTrackStarty_%d_%d",i,j),";Track Start Y (cm);N events", 20, -120, 120);
            hTrackStarty[i][j]->Sumw2();
            hTrackStartz[i][j] = new TH1F(Form("hTrackStartz_%d_%d",i,j),";Track Start Z (cm);N events", 20, 0, 1000);
            hTrackStartz[i][j]->Sumw2();
            hTrackEndx[i][j] = new TH1F(Form("hTrackEndx_%d_%d",i,j),";Track End X (cm);N events", 20, 0, 260);
            hTrackEndx[i][j]->Sumw2();
            hTrackEndy[i][j] = new TH1F(Form("hTrackEndy_%d_%d",i,j),";Track End Y (cm);N events", 20, -120, 120);
            hTrackEndy[i][j]->Sumw2();
            hTrackEndz[i][j] = new TH1F(Form("hTrackEndz_%d_%d",i,j),";Track End Z (cm);N events", 20, 0, 1000);
            hTrackEndz[i][j]->Sumw2();
            hTrackCosz[i][j] = new TH1F(Form("hTrackCosz_%d_%d",i,j),";Track cos(#theta);N events", 20, -1, 1);
            hTrackCosz[i][j]->Sumw2();
            hTrackPhi[i][j] = new TH1F(Form("hTrackPhi_%d_%d",i,j),";Track #phi;N events", 20, -TMath::Pi(), TMath::Pi());
            hTrackPhi[i][j]->Sumw2();
            
            hMul[i][j] = new TH1F(Form("hMul_%d_%d",i,j),"Before Selection;Track Multiplicity;N events",10,0,10);
            hMul[i][j]->Sumw2();
            hTrackLength1[i][j] = new TH1F(Form("hTrackLength1_%d_%d",i,j),"Single Tracks;Track Length (cm);N events",200,0,200);
            hTrackLength1[i][j]->Sumw2();
            hTrackLength2[i][j] = new TH1F(Form("hTrackLength2_%d_%d",i,j),"Track Multiplicity = 3;Track Length (cm);N events",200,0,200);
            hTrackLength2[i][j]->Sumw2();
            hTrackLength3[i][j] = new TH1F(Form("hTrackLength3_%d_%d",i,j),"Track Multiplicity = 2;Track Length (cm);N events",200,0,200);
            hTrackLength3[i][j]->Sumw2();
            hCos1[i][j] = new TH1F(Form("hCos1_%d_%d",i,j),"Track Multiplicity = 2;cos(angle);N events",100,0,1);
            hCos1[i][j]->Sumw2();
            hCos2[i][j] = new TH1F(Form("hCos2_%d_%d",i,j),"Track Multiplicity = 3;cos(angle);N events",100,0,1);
            hCos2[i][j]->Sumw2();
            hCos3[i][j] = new TH1F(Form("hCos3_%d_%d",i,j),"Track Multiplicity > 1;cos(angle);N events",100,0,1);
            hCos3[i][j]->Sumw2();
            hSingleTrackLengthRatioVsdEdxRatio[i][j] = new TH2F(Form("hSingleTrackLengthRatioVsdEdxRatio_%d_%d",i,j),"Track Multiplicity = 1;Projected y length (cm);dE/dx^{high}/dE/dx^{low}",200,0,200,100,0,10);
            hLongTrackdEdxStartEnd[i][j] = new TH2F(Form("hLongTrackdEdxStartEnd_%d_%d",i,j),"Track Multiplicity = 2;dE/dx^{start};dE/dx^{end}",100,0,10,100,0,10);
            hTrackEndy1[i][j] = new TH1F(Form("hTrackEndy1_%d_%d",i,j),"Track Multiplicity = 2;Track End Y (cm);N events", 240, -120, 120);
            hTrackEndy1[i][j]->Sumw2();
            hCos0VsLen1[i][j] = new TH2F(Form("hCos0VsLen1_%d_%d",i,j),"Track Multiplicity > 1;cosy_{track0};Length_{track1}",100,0,1,100,0,100);
            hCosVsLen[i][j] = new TH2F(Form("hCosVsLen_%d_%d",i,j),"Track Multiplicity = 1;Track cosy;Track length",100,0,1,100,0,100);
            for (int k = 0; k<3; ++k){
                hdEdxVsX[i][j][k] = new TH2F(Form("hdEdxVsX_%d_%d_%d",i,j,k),";X (cm);dE/dx (MeV/cm)",40,0,260,100,0,10);
                hdEdxVsX[i][j][k]->Sumw2();
                hdEdxVsXCor[i][j][k] = new TH2F(Form("hdEdxVsXCor_%d_%d_%d",i,j,k),";X (cm);dE/dx (MeV/cm)",40,0,260,100,0,10);
                hdEdxVsXCor[i][j][k]->Sumw2();
            }
        }
    }
    
    ofstream outputfile[nsamples];
    
    const int maxtrks=1000;
    const int maxvtx=1000;
    const int maxflash=1000;
    const int maxtruth=5;
    
    Int_t       run;
    Int_t       subrun;
    Int_t       event;
    Short_t     ntracks;
    Float_t     trklen[maxtrks];
    Float_t     trkstartx[maxtrks];
    Float_t     trkstarty[maxtrks];
    Float_t     trkstartz[maxtrks];
    Float_t     trkendx[maxtrks];
    Float_t     trkendy[maxtrks];
    Float_t     trkendz[maxtrks];
    Float_t     trkstartdcosx[maxtrks];
    Float_t     trkstartdcosy[maxtrks];
    Float_t     trkstartdcosz[maxtrks];
    Float_t     trkenddcosx[maxtrks];
    Float_t     trkenddcosy[maxtrks];
    Float_t     trkenddcosz[maxtrks];
    Float_t     trkphi[maxtrks];
    Int_t       trkorig[maxtrks];
    
    Short_t     ntrkhits[maxtrks][3];
    Float_t     *trkxyz = new Float_t[maxtrks*3*2000*3];
    Float_t     *trkdedx = new Float_t[maxtrks*3*2000];
    Float_t     *trkresrg = new Float_t[maxtrks*3*2000];
    Float_t     *trkthetayz[maxtrks];
    
    Short_t     nvtx;
    Float_t     vtxx[maxvtx];
    Float_t     vtxy[maxvtx];
    Float_t     vtxz[maxvtx];
    
    Int_t       no_flashes;
    Float_t     flash_time[maxflash];
    Float_t     flash_pe[maxflash];
    Float_t     flash_ycenter[maxflash];
    Float_t     flash_zcenter[maxflash];
    
    Int_t           nuPDG_truth[maxtruth];   //[mcevts_truth]
    Int_t           ccnc_truth[maxtruth];   //[mcevts_truth]
    Int_t           mode_truth[maxtruth];   //[mcevts_truth]
    Float_t         enu_truth[maxtruth];   //[mcevts_truth]
    Float_t         Q2_truth[maxtruth];   //[mcevts_truth]
    Float_t         W_truth[maxtruth];   //[mcevts_truth]
    Float_t         X_truth[maxtruth];   //[mcevts_truth]
    Float_t         Y_truth[maxtruth];   //[mcevts_truth]
    Float_t         nuvtxx_truth[maxtruth];   //[mcevts_truth]
    Float_t         nuvtxy_truth[maxtruth];   //[mcevts_truth]
    Float_t         nuvtxz_truth[maxtruth];   //[mcevts_truth]
    Float_t         lep_mom_truth[maxtruth];   //[mcevts_truth]
    Float_t         lep_dcosx_truth[maxtruth];   //[mcevts_truth]
    Float_t         lep_dcosy_truth[maxtruth];   //[mcevts_truth]
    Float_t         lep_dcosz_truth[maxtruth];   //[mcevts_truth]
    
    
    ifstream in;
    for (int isample = 0; isample<nsamples; ++isample){
        
        if (is!=-1 && isample != is) continue;
        
        //for (int isample = 0; isample<1; ++isample){
        TChain ch("analysistree/anatree");
        in.open(filelists[isample].c_str());
        TString filename;
        while(1){
            in>>filename;
            if (!in.good()) break;
            ch.Add(filename);
        }
        in.close();
        in.clear();
        ch.SetBranchStatus("*", 0);
        ch.SetBranchStatus("run", 1);
        ch.SetBranchStatus("subrun", 1);
        ch.SetBranchStatus("event", 1);
        ch.SetBranchStatus(Form("ntracks_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trklen_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkstartx_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkstarty_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkstartz_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkendx_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkendy_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkendz_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkstartdcosx_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkstartdcosy_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkstartdcosz_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkenddcosx_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkenddcosy_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkenddcosz_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkphi_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("ntrkhits_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkxyz_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkdedx_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkresrg_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("trkthetayz_%s",tracklabel.c_str()), 1);
        ch.SetBranchStatus(Form("nvtx_%s",vertexlabel.c_str()), 1);
        ch.SetBranchStatus(Form("vtxx_%s",vertexlabel.c_str()), 1);
        ch.SetBranchStatus(Form("vtxy_%s",vertexlabel.c_str()), 1);
        ch.SetBranchStatus(Form("vtxz_%s",vertexlabel.c_str()), 1);
        ch.SetBranchStatus("no_flashes", 1);
        ch.SetBranchStatus("flash_time", 1);
        ch.SetBranchStatus("flash_pe", 1);
        ch.SetBranchStatus("flash_ycenter", 1);
        ch.SetBranchStatus("flash_zcenter", 1);
        
        ch.SetBranchAddress("run", &run);
        ch.SetBranchAddress("subrun", &subrun);
        ch.SetBranchAddress("event", &event);
        ch.SetBranchAddress(Form("ntracks_%s",tracklabel.c_str()), &ntracks);
        ch.SetBranchAddress(Form("trklen_%s",tracklabel.c_str()), trklen);
        ch.SetBranchAddress(Form("trkstartx_%s",tracklabel.c_str()), trkstartx);
        ch.SetBranchAddress(Form("trkstarty_%s",tracklabel.c_str()), trkstarty);
        ch.SetBranchAddress(Form("trkstartz_%s",tracklabel.c_str()), trkstartz);
        ch.SetBranchAddress(Form("trkendx_%s",tracklabel.c_str()), trkendx);
        ch.SetBranchAddress(Form("trkendy_%s",tracklabel.c_str()), trkendy);
        ch.SetBranchAddress(Form("trkendz_%s",tracklabel.c_str()), trkendz);
        ch.SetBranchAddress(Form("trkstartdcosx_%s",tracklabel.c_str()), trkstartdcosx);
        ch.SetBranchAddress(Form("trkstartdcosy_%s",tracklabel.c_str()), trkstartdcosy);
        ch.SetBranchAddress(Form("trkstartdcosz_%s",tracklabel.c_str()), trkstartdcosz);
        ch.SetBranchAddress(Form("trkenddcosx_%s",tracklabel.c_str()), trkenddcosx);
        ch.SetBranchAddress(Form("trkenddcosy_%s",tracklabel.c_str()), trkenddcosy);
        ch.SetBranchAddress(Form("trkenddcosz_%s",tracklabel.c_str()), trkenddcosz);
        ch.SetBranchAddress(Form("trkphi_%s",tracklabel.c_str()), trkphi);
        ch.SetBranchAddress(Form("ntrkhits_%s",tracklabel.c_str()), ntrkhits);
        ch.SetBranchAddress(Form("trkxyz_%s",tracklabel.c_str()), trkxyz);
        ch.SetBranchAddress(Form("trkdedx_%s",tracklabel.c_str()), trkdedx);
        ch.SetBranchAddress(Form("trkresrg_%s",tracklabel.c_str()), trkresrg);
        ch.SetBranchAddress(Form("trkthetayz_%s",tracklabel.c_str()), trkthetayz);
        ch.SetBranchAddress(Form("nvtx_%s",vertexlabel.c_str()), &nvtx);
        ch.SetBranchAddress(Form("vtxx_%s",vertexlabel.c_str()), vtxx);
        ch.SetBranchAddress(Form("vtxy_%s",vertexlabel.c_str()), vtxy);
        ch.SetBranchAddress(Form("vtxz_%s",vertexlabel.c_str()), vtxz);
        ch.SetBranchAddress("no_flashes", &no_flashes );
        ch.SetBranchAddress("flash_time", flash_time);
        ch.SetBranchAddress("flash_pe", flash_pe);
        ch.SetBranchAddress("flash_ycenter", flash_ycenter);
        ch.SetBranchAddress("flash_zcenter", flash_zcenter);
        
        if (isample>1){
            ch.SetBranchStatus(Form("trkorig_%s",tracklabel.c_str()), 1);
            ch.SetBranchStatus("nuPDG_truth", 1);
            ch.SetBranchStatus("ccnc_truth", 1);
            ch.SetBranchStatus("mode_truth", 1);
            ch.SetBranchStatus("nuvtxx_truth", 1);
            ch.SetBranchStatus("nuvtxy_truth", 1);
            ch.SetBranchStatus("nuvtxz_truth", 1);
            ch.SetBranchAddress(Form("trkorig_%s",tracklabel.c_str()), trkorig);
            ch.SetBranchAddress("nuPDG_truth",nuPDG_truth);
            ch.SetBranchAddress("ccnc_truth",ccnc_truth);
            ch.SetBranchAddress("mode_truth",mode_truth);
            ch.SetBranchAddress("nuvtxx_truth",nuvtxx_truth);
            ch.SetBranchAddress("nuvtxy_truth",nuvtxy_truth);
            ch.SetBranchAddress("nuvtxz_truth",nuvtxz_truth);
        }
        
        outputfile[isample].open(Form("%s_%s_%s.txt",samplename[isample].c_str(), tracklabel.c_str(), vertexlabel.c_str()));
        
        cout<<"Processing sample "<<isample<<" "<<samplename[isample]<<endl;
        Long64_t nentries = ch.GetEntries();
        for (Long64_t jentry = 0; jentry<nentries; ++jentry){
            if (jentry%10000==0) cout<<"processing event "<<jentry<<"/"<<nentries<<endl;
            ch.GetEntry(jentry);
            ++nevents[isample][0];
            
            //check the flash info
            float FlashPEmax=0;
            int NuFlashID=-1;
            
            //Look for most energetic flash above 50 PE in the beam window
            float t0 = 0;
            float t1 = 0;
            if (isample == 0){
                t0 = 3.3;
                t1 = 4.9;
            }
            else if (isample == 1){
                t0 = 3.65;
                t1 = 5.25;
            }
            else if (isample == 2){
                t0 = 3.55;
                t1 = 5.15;
            }
            else if (isample == 3){
                t0 = 3.2;
                t1 = 4.3;
            }
            else {
                cout<<"Unknow sample, beam window not set."<<endl;
            }
            for (int i = 0; i<no_flashes; ++i){
                if (flash_pe[i]>50 && flash_time[i]>t0 && flash_time[i]<t1){
                    if (flash_pe[i]>FlashPEmax){
                        FlashPEmax = flash_pe[i];
                        NuFlashID = i;
                    }
                }
            }
            
            //Did not find the desired flash, continue
            if (NuFlashID == -1){
                if (debug) cout<<"Did not find flash in the beam window."<<endl;
                continue;
            }
            ++nevents[isample][1];
            
            //Match each track with the selected flash
            vector<bool> trackflashmatch(ntracks);
            bool foundtrackflashmatch = false;
            
            float TaggedFlashYCenter = flash_ycenter[NuFlashID];
            float TaggedFlashZCenter = flash_zcenter[NuFlashID];
            for (int i = 0; i<ntracks; ++i){
                float FlashTrackDis = 1e10;
                if ((trkstartz[i]<TaggedFlashZCenter && trkendz[i]>TaggedFlashZCenter)||
                    (trkstartz[i]>TaggedFlashZCenter && trkendz[i]<TaggedFlashZCenter)){
                    FlashTrackDis = 0;
                }
                else{
                    FlashTrackDis = std::min(std::abs(trkstartz[i] - TaggedFlashZCenter),
                                             std::abs(trkendz[i] - TaggedFlashZCenter));
                }
                if (FlashTrackDis<70){
                    trackflashmatch[i] = true;
                    foundtrackflashmatch = true;
                }
                else{
                    trackflashmatch[i] = false;
                }
            }
            if (!foundtrackflashmatch) {
                if (debug) cout<<"Did not find any tracks matching flash."<<endl;
                continue;
            }
            ++nevents[isample][2];
            
            //Match tracks with vertices
            if (!nvtx) {
                if (debug) cout<<"No vertex found"<<endl;
                continue; //No vertex found
            }
            vector<vector<int>> trkindex(nvtx);
            vector<vector<bool>> fliptrack(nvtx);
            for (int i = 0; i<nvtx; ++i){
                if (!inFV(vtxx[i], vtxy[i], vtxz[i])) continue;
                for (int j = 0; j<ntracks; ++j){
                    float vtxtrkStartDis = sqrt(pow(trkstartx[j]-vtxx[i],2)+
                                                pow(trkstarty[j]-vtxy[i],2)+
                                                pow(trkstartz[j]-vtxz[i],2));
                    float vtxtrkEndDis = sqrt(pow(trkendx[j]-vtxx[i],2)+
                                              pow(trkendy[j]-vtxy[i],2)+
                                              pow(trkendz[j]-vtxz[i],2));
                    float vtxtrkDis = min(vtxtrkStartDis, vtxtrkEndDis);
                    if (vtxtrkDis<3){
                        trkindex[i].push_back(j);
                        if (vtxtrkEndDis<vtxtrkStartDis){
                            fliptrack[i].push_back(true);
                        }
                        else{
                            fliptrack[i].push_back(false);
                        }
                    }
                }//Loope over all tracks
            }//Loop over all vertices
            
            //calculate average dE/dx near the track start and track end
            vector<float> trkStartdEdx(ntracks);
            vector<float> trkEnddEdx(ntracks);
            for (int i = 0; i<ntracks; ++i){
                
                int totalnhits0 = ntrkhits[i][0];
                int totalnhits1 = ntrkhits[i][1];
                int totalnhits2 = ntrkhits[i][2];
                int totalnhits=totalnhits0;
                int iplane=0;
                if(totalnhits<totalnhits1){
                    totalnhits=totalnhits1;
                    iplane=1;
                }
                if(totalnhits<totalnhits2){
                    totalnhits=totalnhits2;
                    iplane=2;
                }
                
                float sumdEdxStart=0;
                float sumdEdxEnd=0;
                
                int MaxHits=0;
                if(totalnhits>=20){
                    MaxHits=10;
                }
                else if(totalnhits>0){
                    MaxHits=totalnhits/2;}
                for(int ihit=0;ihit<MaxHits;ihit++){
                    sumdEdxStart += trkdedx[i*3*2000+iplane*2000+ihit]*scaledEdx(trkxyz[i*3*2000*3+iplane*2000*3+ihit*3+0], iplane, isample<2);
                    sumdEdxEnd += trkdedx[i*3*2000+iplane*2000+totalnhits-ihit-1]*scaledEdx(trkxyz[i*3*2000*3+iplane*2000*3+(totalnhits-ihit-1)*3+0], iplane, isample<2);
                }
                if (debug) cout<<trkxyz[i*3*2000*3+iplane*2000*3+0*3+0]<<" "
                <<trkxyz[i*3*2000*3+iplane*2000*3+0*3+1]<<" "
                <<trkxyz[i*3*2000*3+iplane*2000*3+0*3+2]<<" "
                <<trkxyz[i*3*2000*3+iplane*2000*3+(totalnhits-1)*3+0]<<" "
                <<trkxyz[i*3*2000*3+iplane*2000*3+(totalnhits-1)*3+1]<<" "
                <<trkxyz[i*3*2000*3+iplane*2000*3+(totalnhits-1)*3+2]<<" "
                <<trkstartx[i]<<" "<<trkstarty[i]<<" "<<trkstartz[i]<<" "
                <<trkendx[i]<<" "<<trkendy[i]<<" "<<trkendz[i]<<endl;
                if (sqrt(pow(trkxyz[i*3*2000*3+iplane*2000*3+0*3+0]-trkstartx[i],2)+
                         pow(trkxyz[i*3*2000*3+iplane*2000*3+0*3+1]-trkstarty[i],2)+
                         pow(trkxyz[i*3*2000*3+iplane*2000*3+0*3+2]-trkstartz[i],2))>
                    sqrt(pow(trkxyz[i*3*2000*3+iplane*2000*3+0*3+0]-trkendx[i],2)+
                         pow(trkxyz[i*3*2000*3+iplane*2000*3+0*3+1]-trkendy[i],2)+
                         pow(trkxyz[i*3*2000*3+iplane*2000*3+0*3+2]-trkendz[i],2))){
                        swap(sumdEdxEnd, sumdEdxStart);
                    }
                if(MaxHits>0 && sumdEdxEnd>0 && sumdEdxStart>0){
                    trkStartdEdx[i] = sumdEdxStart/MaxHits;
                    trkEnddEdx[i] = sumdEdxEnd/MaxHits;
                }
                if (debug) cout<<"Trkid = "<<i<<" MaxHits = "<<MaxHits<<" sumdEdxStart "<<sumdEdxStart<<" sumdEdxEnd "<<sumdEdxEnd<<" iplane "<<iplane<<endl;
            }//Loop over tracks
            
            //Examine tracks around each vertex and select neutrino candidates
            vector<bool> nuvtx(nvtx);
            vector<float> cosangle(nvtx); //angle between two longest tracks
            vector<float> dcosylong(nvtx); //dcosy of the longest track
            vector<float> trklen2nd(nvtx); //track length of the second longest track
            for (int i = 0; i<nvtx; ++i){
                
                //Check if there are tracks associated with this vertex
                if (!trkindex[i].size()) {
                    //no tracks associated with vertex
                    nuvtx[i] = false;
                    if (debug) cout<<"ivtx = "<<i<<" no tracks associated with this vertex."<<endl;
                    continue;
                }
                
                //Check if at least one track matches the flash
                bool flashmatch = false;
                for (size_t j = 0; j<trkindex[i].size(); ++j){
                    if (trackflashmatch[trkindex[i][j]]){
                        flashmatch = true;
                    }
                }
                if (!flashmatch){
                    //no tracks matched to the flash around this vertex
                    nuvtx[i] = false;
                    if (debug) cout<<"ivtx = "<<i<<" no tracks around the vertex matched to the flash."<<endl;
                    continue;
                }
                
                //flag for cosmic/neutrino, -1 is data, 1 is cosmic, 2 is neutrino
                int itype = -1;
                if (isample>1){
                    itype = 1; //cosmic
                    for (size_t j = 0; j<trkindex[i].size(); ++j){
                        if (trkorig[trkindex[i][j]]==1) itype = 2; //neutrino
                    }
                }
                
                //Fill track multiplicity information
                hMul[isample][0]->Fill(trkindex[i].size());
                if (itype!=-1) hMul[isample][itype]->Fill(trkindex[i].size());
                
                //study if mult>=2
                if (trkindex[i].size()>1){
                    //find two longest tracks
                    int j0 = -1;
                    int j1 = -1;
                    float trklen0 = -1;
                    float trklen1 = -1;
                    //find the highest track
                    int jhigh = -1;
                    float highy= -1000;
                    for (size_t j = 0; j<trkindex[i].size(); ++j){
                        if (trklen[trkindex[i][j]]>trklen0){
                            j1 = j0;
                            trklen1 = trklen0;
                            j0 = j;
                            trklen0 = trklen[trkindex[i][j]];
                        }
                        else if (trklen[trkindex[i][j]]>trklen1){
                            j1 = j;
                            trklen1 = trklen[trkindex[i][j]];
                        }
                        if (trkstarty[trkindex[i][j]]>highy){
                            highy = trkstarty[trkindex[i][j]];
                            jhigh = j;
                        }
                        if (trkendy[trkindex[i][j]]>highy){
                            highy = trkendy[trkindex[i][j]];
                            jhigh = j;
                        }
                    }
                    float dcosx0 = trkstartdcosx[trkindex[i][j0]];
                    float dcosy0 = trkstartdcosy[trkindex[i][j0]];
                    float dcosz0 = trkstartdcosz[trkindex[i][j0]];
                    float dcosx1 = trkstartdcosx[trkindex[i][j1]];
                    float dcosy1 = trkstartdcosy[trkindex[i][j1]];
                    float dcosz1 = trkstartdcosz[trkindex[i][j1]];
                    if (fliptrack[i][j0]){
                        dcosx0 = trkenddcosx[trkindex[i][j0]];
                        dcosy0 = trkenddcosy[trkindex[i][j0]];
                        dcosz0 = trkenddcosz[trkindex[i][j0]];
                    }
                    if (fliptrack[i][j1]){
                        dcosx1 = trkenddcosx[trkindex[i][j1]];
                        dcosy1 = trkenddcosy[trkindex[i][j1]];
                        dcosz1 = trkenddcosz[trkindex[i][j1]];
                    }
                    cosangle[i] = abs(dcosx0*dcosx1+dcosy0*dcosy1+dcosz0*dcosz1);
                    hCos3[isample][0]->Fill(cosangle[i]);
                    if (itype!=-1) hCos3[isample][itype]->Fill(cosangle[i]);
                    if (j0 == jhigh){
                        dcosylong[i] = fliptrack[i][j0]?abs(trkstartdcosy[trkindex[i][j0]]):abs(trkenddcosy[trkindex[i][j0]]);
                        trklen2nd[i] = trklen[trkindex[i][j1]];
                        hCos0VsLen1[isample][0]->Fill(dcosylong[i],trklen2nd[i]);
                        if (itype!=-1) hCos0VsLen1[isample][itype]->Fill(dcosylong[i],trklen2nd[i]);
                    }
                    if (cosangle[i]>0.9) {
                        nuvtx[i] = false;
                        continue;
                    }
                    if (j0 == jhigh){
                        if (dcosylong[i]>0.6&&trklen2nd[i]<30){
                            nuvtx[i] = false;
                            continue;
                        }
                    }
                }
                
                //If there are more than 3 tracks, accept the vertex
                if (trkindex[i].size()>3){
                    nuvtx[i] = true;
                    if (debug) cout<<"ivtx = "<<i<<" track multiplicity = "<<trkindex[i].size()<<endl;
                    continue;
                }
                
                //Single track
                if (trkindex[i].size()==1){
                    nuvtx[i] = false; //will update later
                    //Only select contained track
                    int itrk = trkindex[i][0];
                    //Track is fully contained
                    if (inFV(trkstartx[itrk],
                             trkstarty[itrk],
                             trkstartz[itrk])&&
                        inFV(trkendx[itrk],
                             trkendy[itrk],
                             trkendz[itrk])){
                            hTrackLength1[isample][0]->Fill(trklen[itrk]);
                            if (itype!=-1) hTrackLength1[isample][itype]->Fill(trklen[itrk]);
                            hCosVsLen[isample][0]->Fill(abs(trkstartdcosy[itrk]), trklen[itrk]);
                            if (itype!=-1) hCosVsLen[isample][itype]->Fill(abs(trkstartdcosy[itrk]), trklen[itrk]);
                            if (abs(trkstartdcosy[itrk])>0.7) {
                                nuvtx[i] = false;
                                continue;
                            }
                            //At least 40 cm
                            if (trklen[itrk]>40){
                                float TrackLengthYNu = trklen[itrk]*abs(trkstartdcosy[itrk]);
                                float LongSingleTrackdEdxRatio = -999;
                                if (trkstarty[itrk]>trkendy[itrk]){
                                    LongSingleTrackdEdxRatio = trkStartdEdx[itrk]/trkEnddEdx[itrk];
                                }
                                else{
                                    LongSingleTrackdEdxRatio = trkEnddEdx[itrk]/trkStartdEdx[itrk];
                                }
                                if (debug) cout<<"LongSingleTrackdEdxRatio = "<<LongSingleTrackdEdxRatio<<" TrackLengthYNu "<<TrackLengthYNu<<" trkStartdEdx "<<trkStartdEdx[itrk]<<" trkEnddEdx "<<trkEnddEdx[itrk]<<endl;
                                if(LongSingleTrackdEdxRatio>1.5 || (TrackLengthYNu<=25 && LongSingleTrackdEdxRatio<=1.5)){
                                    nuvtx[i] = true;
                                }
                                hSingleTrackLengthRatioVsdEdxRatio[isample][0]->Fill(TrackLengthYNu, LongSingleTrackdEdxRatio);
                                if (itype!=-1) hSingleTrackLengthRatioVsdEdxRatio[isample][itype]->Fill(TrackLengthYNu, LongSingleTrackdEdxRatio);
                            }//At least 40 cm
                        }//Track is contained
                    if (debug) cout<<"ivtx = "<<i<<" single track "<<nuvtx[i]<<endl;
                    continue;
                }//Single track
                
                //Multiplicity = 2
                if (trkindex[i].size()==2){
                    bool isMichel = false;
                    if (debug) cout<<trklen[trkindex[i][0]]<<" "<<fliptrack[i][0]<<" "<<trkStartdEdx[trkindex[i][0]]<<" "<<trkEnddEdx[trkindex[i][0]]<<" "<<trkstarty[trkindex[i][0]]<<" "<<trkendy[trkindex[i][0]]<<" "<<trklen[trkindex[i][1]]<<" "<<fliptrack[i][1]<<" "<<trkStartdEdx[trkindex[i][1]]<<" "<<trkEnddEdx[trkindex[i][1]]<<" "<<trkstarty[trkindex[i][1]]<<" "<<trkendy[trkindex[i][1]]<<endl;
                    float trkstartdedx0 = 0;
                    float trkenddedx0 = 0;
                    float trkendy0 = 0;
                    float trklen1 = 0;
                    if (trklen[trkindex[i][0]]>trklen[trkindex[i][1]]){//first track is longer
                        if (!fliptrack[i][0]){
                            trkstartdedx0 = trkStartdEdx[trkindex[i][0]];
                            trkenddedx0 = trkEnddEdx[trkindex[i][0]];
                            trkendy0 = trkendy[trkindex[i][0]];
                            trklen1 = trklen[trkindex[i][1]];
                        }
                        else{
                            trkstartdedx0 = trkEnddEdx[trkindex[i][0]];
                            trkenddedx0 = trkStartdEdx[trkindex[i][0]];
                            trkendy0 = trkstarty[trkindex[i][0]];
                            trklen1 = trklen[trkindex[i][1]];
                        }
                    }//first track is longer
                    else{//second track is longer
                        if (!fliptrack[i][1]){
                            trkstartdedx0 = trkStartdEdx[trkindex[i][1]];
                            trkenddedx0 = trkEnddEdx[trkindex[i][1]];
                            trkendy0 = trkendy[trkindex[i][1]];
                            trklen1 = trklen[trkindex[i][0]];
                        }
                        else{
                            trkstartdedx0 = trkEnddEdx[trkindex[i][1]];
                            trkenddedx0 = trkStartdEdx[trkindex[i][1]];
                            trkendy0 = trkstarty[trkindex[i][1]];
                            trklen1 = trklen[trkindex[i][0]];
                        }
                    }//second track is longer
                    hLongTrackdEdxStartEnd[isample][0]->Fill(trkstartdedx0,trkenddedx0);
                    hTrackLength3[isample][0]->Fill(trklen1);
                    hTrackEndy1[isample][0]->Fill(trkendy0);
                    if (itype!=-1){
                        hLongTrackdEdxStartEnd[isample][itype]->Fill(trkstartdedx0,trkenddedx0);
                        hTrackLength3[isample][itype]->Fill(trklen1);
                        hTrackEndy1[isample][itype]->Fill(trkendy0);
                    }
                    if (((trkstartdedx0>trkenddedx0&&
                          trkstartdedx0>2.5&&trkenddedx0<4)||
                         trkendy0>96.5)&&trklen1<30)
                    isMichel = true;
                    
                    hCos1[isample][0]->Fill(abs(cosangle[i]));
                    if (itype!=-1) hCos1[isample][itype]->Fill(abs(cosangle[i]));
                    
                    if (abs(cosangle[i])>0.9) isMichel = true;
                    if (isMichel){
                        nuvtx[i] = false;
                    }
                    else{
                        nuvtx[i] = true;
                    }
                    if (debug) cout<<"ivtx = "<<i<<" mul = 2 "<<nuvtx[i]<<endl;
                    continue;
                }//Multiplicity = 2
                
                //mult = 3
                if (trkindex[i].size() == 3){
                    //find the min track length
                    float mintrklen = 1e10;
                    int ishortest = -1;
                    for (size_t j = 0; j<trkindex[i].size(); ++j){
                        if (trklen[trkindex[i][j]]<mintrklen){
                            mintrklen = trklen[trkindex[i][j]];
                            ishortest = trkindex[i][j];
                        }
                    }
                    hTrackLength2[isample][0]->Fill(mintrklen);
                    if (itype!=-1) hTrackLength2[isample][itype]->Fill(mintrklen);
                    
                    hCos2[isample][0]->Fill(abs(cosangle[i]));
                    if (itype!=-1) hCos2[isample][itype]->Fill(abs(cosangle[i]));
                    
                    if (cosangle[i]<0.9) nuvtx[i] = true;
                    else nuvtx[i] = false;
                    continue;
                } // mult = 3
                
            }//Loop over all vertices
            
            //Find the longest track
            int ivtx = -1;
            int itrk = -1;
            float longesttracklength = -1;
            for (int i = 0; i<nvtx; ++i){
                if (!nuvtx[i]) continue;
                for (size_t j = 0; j<trkindex[i].size(); ++j){
                    if (debug) cout<<"ivtx = "<<i<<" trkid = "<<trkindex[i][j]<<" tracklen "<<trklen[trkindex[i][j]]<<" trackflashmatch "<<trackflashmatch[trkindex[i][j]]<<endl;
                    if (trklen[trkindex[i][j]]>longesttracklength&&
                        trackflashmatch[trkindex[i][j]]&&
                        trklen[trkindex[i][j]]>15){
                        longesttracklength = trklen[trkindex[i][j]];
                        ivtx = i;
                        itrk = j;
                    }
                }
            }//Loop over all vertices
            if (ivtx!=-1 && itrk!=-1){
                outputfile[isample]<<run<<" "<<subrun<<" "<<event<<" "<<ivtx<<" "<<trkindex[ivtx][itrk]<<" "<<trkindex[ivtx].size()<<endl;
                if (trkindex[ivtx].size()==1) ++nevents[isample][3];
                else ++nevents[isample][4];
                hTrackLength[isample][0]->Fill(trklen[trkindex[ivtx][itrk]]);
                hTrackMul[isample][0]->Fill(trkindex[ivtx].size());
                hVtxx[isample][0]->Fill(vtxx[ivtx]);
                hVtxy[isample][0]->Fill(vtxy[ivtx]);
                hVtxz[isample][0]->Fill(vtxz[ivtx]);
                if (!fliptrack[ivtx][itrk]){
                    hTrackStartx[isample][0]->Fill(trkstartx[trkindex[ivtx][itrk]]);
                    hTrackStarty[isample][0]->Fill(trkstarty[trkindex[ivtx][itrk]]);
                    hTrackStartz[isample][0]->Fill(trkstartz[trkindex[ivtx][itrk]]);
                    hTrackEndx[isample][0]->Fill(trkendx[trkindex[ivtx][itrk]]);
                    hTrackEndy[isample][0]->Fill(trkendy[trkindex[ivtx][itrk]]);
                    hTrackEndz[isample][0]->Fill(trkendz[trkindex[ivtx][itrk]]);
                    hTrackCosz[isample][0]->Fill(trkstartdcosz[trkindex[ivtx][itrk]]);
                    hTrackPhi[isample][0]->Fill(trkphi[trkindex[ivtx][itrk]]);
                }
                else{
                    hTrackStartx[isample][0]->Fill(trkendx[trkindex[ivtx][itrk]]);
                    hTrackStarty[isample][0]->Fill(trkendy[trkindex[ivtx][itrk]]);
                    hTrackStartz[isample][0]->Fill(trkendz[trkindex[ivtx][itrk]]);
                    hTrackEndx[isample][0]->Fill(trkstartx[trkindex[ivtx][itrk]]);
                    hTrackEndy[isample][0]->Fill(trkstarty[trkindex[ivtx][itrk]]);
                    hTrackEndz[isample][0]->Fill(trkstartz[trkindex[ivtx][itrk]]);
                    hTrackCosz[isample][0]->Fill(-trkenddcosz[trkindex[ivtx][itrk]]);
                    TVector3 v(-trkenddcosx[trkindex[ivtx][itrk]],-trkenddcosy[trkindex[ivtx][itrk]],-trkenddcosz[trkindex[ivtx][itrk]]);
                    hTrackPhi[isample][0]->Fill(v.Phi());
                }
                int iplane = 0;
                if (ntrkhits[trkindex[ivtx][itrk]][1]>ntrkhits[trkindex[ivtx][itrk]][iplane]) iplane = 1;
                if (ntrkhits[trkindex[ivtx][itrk]][2]>ntrkhits[trkindex[ivtx][itrk]][iplane]) iplane = 2;
                for (int l = 0; l<ntrkhits[trkindex[ivtx][itrk]][iplane]; ++l){
                    hdEdxVsX[isample][0][iplane]->Fill(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0],
                                                       trkdedx[trkindex[ivtx][itrk]*3*2000+iplane*2000+l]);
                    hdEdxVsXCor[isample][0][iplane]->Fill(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0],
                                                          trkdedx[trkindex[ivtx][itrk]*3*2000+iplane*2000+l]*scaledEdx(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0], iplane, isample<2));
                }
                if (trkindex[ivtx].size()==1){
                    hTrackLength[isample][10]->Fill(trklen[trkindex[ivtx][itrk]]);
                    hTrackMul[isample][10]->Fill(trkindex[ivtx].size());
                    hVtxx[isample][10]->Fill(vtxx[ivtx]);
                    hVtxy[isample][10]->Fill(vtxy[ivtx]);
                    hVtxz[isample][10]->Fill(vtxz[ivtx]);
                    for (int l = 0; l<ntrkhits[trkindex[ivtx][itrk]][iplane]; ++l){
                        hdEdxVsX[isample][10][iplane]->Fill(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0],
                                                            trkdedx[trkindex[ivtx][itrk]*3*2000+iplane*2000+l]);
                        hdEdxVsXCor[isample][10][iplane]->Fill(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0],
                                                               trkdedx[trkindex[ivtx][itrk]*3*2000+iplane*2000+l]*scaledEdx(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0], iplane, isample<2));
                    }
                    
                    if (!fliptrack[ivtx][itrk]){
                        hTrackStartx[isample][10]->Fill(trkstartx[trkindex[ivtx][itrk]]);
                        hTrackStarty[isample][10]->Fill(trkstarty[trkindex[ivtx][itrk]]);
                        hTrackStartz[isample][10]->Fill(trkstartz[trkindex[ivtx][itrk]]);
                        hTrackEndx[isample][10]->Fill(trkendx[trkindex[ivtx][itrk]]);
                        hTrackEndy[isample][10]->Fill(trkendy[trkindex[ivtx][itrk]]);
                        hTrackEndz[isample][10]->Fill(trkendz[trkindex[ivtx][itrk]]);
                        hTrackCosz[isample][10]->Fill(trkstartdcosz[trkindex[ivtx][itrk]]);
                        hTrackPhi[isample][10]->Fill(trkphi[trkindex[ivtx][itrk]]);
                    }
                    else{
                        hTrackStartx[isample][10]->Fill(trkendx[trkindex[ivtx][itrk]]);
                        hTrackStarty[isample][10]->Fill(trkendy[trkindex[ivtx][itrk]]);
                        hTrackStartz[isample][10]->Fill(trkendz[trkindex[ivtx][itrk]]);
                        hTrackEndx[isample][10]->Fill(trkstartx[trkindex[ivtx][itrk]]);
                        hTrackEndy[isample][10]->Fill(trkstarty[trkindex[ivtx][itrk]]);
                        hTrackEndz[isample][10]->Fill(trkstartz[trkindex[ivtx][itrk]]);
                        hTrackCosz[isample][10]->Fill(-trkenddcosz[trkindex[ivtx][itrk]]);
                        TVector3 v(-trkenddcosx[trkindex[ivtx][itrk]],-trkenddcosy[trkindex[ivtx][itrk]],-trkenddcosz[trkindex[ivtx][itrk]]);
                        hTrackPhi[isample][10]->Fill(v.Phi());
                    }
                }          
                else{
                    hTrackLength[isample][11]->Fill(trklen[trkindex[ivtx][itrk]]);
                    hTrackMul[isample][11]->Fill(trkindex[ivtx].size());
                    hVtxx[isample][11]->Fill(vtxx[ivtx]);
                    hVtxy[isample][11]->Fill(vtxy[ivtx]);
                    hVtxz[isample][11]->Fill(vtxz[ivtx]);
                    for (int l = 0; l<ntrkhits[trkindex[ivtx][itrk]][iplane]; ++l){
                        hdEdxVsX[isample][11][iplane]->Fill(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0],
                                                            trkdedx[trkindex[ivtx][itrk]*3*2000+iplane*2000+l]);
                        hdEdxVsXCor[isample][11][iplane]->Fill(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0],
                                                               trkdedx[trkindex[ivtx][itrk]*3*2000+iplane*2000+l]*scaledEdx(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0], iplane, isample<2));
                    }
                    if (!fliptrack[ivtx][itrk]){
                        hTrackStartx[isample][11]->Fill(trkstartx[trkindex[ivtx][itrk]]);
                        hTrackStarty[isample][11]->Fill(trkstarty[trkindex[ivtx][itrk]]);
                        hTrackStartz[isample][11]->Fill(trkstartz[trkindex[ivtx][itrk]]);
                        hTrackEndx[isample][11]->Fill(trkendx[trkindex[ivtx][itrk]]);
                        hTrackEndy[isample][11]->Fill(trkendy[trkindex[ivtx][itrk]]);
                        hTrackEndz[isample][11]->Fill(trkendz[trkindex[ivtx][itrk]]);
                        hTrackCosz[isample][11]->Fill(trkstartdcosz[trkindex[ivtx][itrk]]);
                        hTrackPhi[isample][11]->Fill(trkphi[trkindex[ivtx][itrk]]);
                    }
                    else{
                        hTrackStartx[isample][11]->Fill(trkendx[trkindex[ivtx][itrk]]);
                        hTrackStarty[isample][11]->Fill(trkendy[trkindex[ivtx][itrk]]);
                        hTrackStartz[isample][11]->Fill(trkendz[trkindex[ivtx][itrk]]);
                        hTrackEndx[isample][11]->Fill(trkstartx[trkindex[ivtx][itrk]]);
                        hTrackEndy[isample][11]->Fill(trkstarty[trkindex[ivtx][itrk]]);
                        hTrackEndz[isample][11]->Fill(trkstartz[trkindex[ivtx][itrk]]);
                        hTrackCosz[isample][11]->Fill(-trkenddcosz[trkindex[ivtx][itrk]]);
                        TVector3 v(-trkenddcosx[trkindex[ivtx][itrk]],-trkenddcosy[trkindex[ivtx][itrk]],-trkenddcosz[trkindex[ivtx][itrk]]);
                        hTrackPhi[isample][11]->Fill(v.Phi());
                    }
                }
                if (isample>1){ //MC
                    int itype = -1;
                    if (trkorig[trkindex[ivtx][itrk]]!=1) itype = 1; //cosmic background
                    else if (!inFV(nuvtxx_truth[0], nuvtxy_truth[0], nuvtxz_truth[0])) itype = 2; //out of FV neutrino
                    else if (ccnc_truth[0] == 1) itype = 3; //NC background
                    else if (abs(nuPDG_truth[0]) == 12) itype = 4; //nue
                    else if (nuPDG_truth[0] == -14) itype = 5; //anti-numu
                    else if (mode_truth[0] == 0) itype = 6; //QE
                    else if (mode_truth[0] == 1) itype = 7; //RES
                    else if (mode_truth[0] == 2) itype = 8; //DIS
                    else if (mode_truth[0] == 3) itype = 9; //Coh
                    if (itype != -1){
                        hTrackLength[isample][itype]->Fill(trklen[trkindex[ivtx][itrk]]);
                        hTrackMul[isample][itype]->Fill(trkindex[ivtx].size());
                        hVtxx[isample][itype]->Fill(vtxx[ivtx]);
                        hVtxy[isample][itype]->Fill(vtxy[ivtx]);
                        hVtxz[isample][itype]->Fill(vtxz[ivtx]);
                        for (int l = 0; l<ntrkhits[trkindex[ivtx][itrk]][iplane]; ++l){
                            hdEdxVsX[isample][itype][iplane]->Fill(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0],
                                                                   trkdedx[trkindex[ivtx][itrk]*3*2000+iplane*2000+l]);
                            hdEdxVsXCor[isample][itype][iplane]->Fill(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0],
                                                                      trkdedx[trkindex[ivtx][itrk]*3*2000+iplane*2000+l]*scaledEdx(trkxyz[trkindex[ivtx][itrk]*3*2000*3+iplane*2000*3+l*3+0], iplane, isample<2));
                        }
                        if (!fliptrack[ivtx][itrk]){
                            hTrackStartx[isample][itype]->Fill(trkstartx[trkindex[ivtx][itrk]]);
                            hTrackStarty[isample][itype]->Fill(trkstarty[trkindex[ivtx][itrk]]);
                            hTrackStartz[isample][itype]->Fill(trkstartz[trkindex[ivtx][itrk]]);
                            hTrackEndx[isample][itype]->Fill(trkendx[trkindex[ivtx][itrk]]);
                            hTrackEndy[isample][itype]->Fill(trkendy[trkindex[ivtx][itrk]]);
                            hTrackEndz[isample][itype]->Fill(trkendz[trkindex[ivtx][itrk]]);
                            hTrackCosz[isample][itype]->Fill(trkstartdcosz[trkindex[ivtx][itrk]]);
                            hTrackPhi[isample][itype]->Fill(trkphi[trkindex[ivtx][itrk]]);
                        }
                        else{
                            hTrackStartx[isample][itype]->Fill(trkendx[trkindex[ivtx][itrk]]);
                            hTrackStarty[isample][itype]->Fill(trkendy[trkindex[ivtx][itrk]]);
                            hTrackStartz[isample][itype]->Fill(trkendz[trkindex[ivtx][itrk]]);
                            hTrackEndx[isample][itype]->Fill(trkstartx[trkindex[ivtx][itrk]]);
                            hTrackEndy[isample][itype]->Fill(trkstarty[trkindex[ivtx][itrk]]);
                            hTrackEndz[isample][itype]->Fill(trkstartz[trkindex[ivtx][itrk]]);
                            hTrackCosz[isample][itype]->Fill(-trkenddcosz[trkindex[ivtx][itrk]]);
                            TVector3 v(-trkenddcosx[trkindex[ivtx][itrk]],-trkenddcosy[trkindex[ivtx][itrk]],-trkenddcosz[trkindex[ivtx][itrk]]);
                            hTrackPhi[isample][itype]->Fill(v.Phi());
                        }
                    }
                }
            }
        }//Loop over all entries
        outputfile[isample].close();
    }//Loop over all samples
    
    OutHist->Write();
    OutHist->Close();
    for (int i = 0; i<nsamples; ++i){
        for (int j = 0; j<ncuts; ++j){
            cout<<"isample = "<<i<<" icut = "<<j<<" "<<nevents[i][j]<<endl;
        }
    }
    
}
