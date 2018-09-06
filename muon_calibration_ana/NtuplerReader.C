#define NtuplerReader_cxx
#include "NtuplerReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>

using namespace std;

void NtuplerReader::Loop( TString path, TString num, int nvtx_sel )
{

   TString nvtx_sel_Str;
   nvtx_sel_Str.Form("%d",nvtx_sel);
   TFile* outputFile=TFile::Open(path+"/tree_Mmumu_nvtx_sel"+nvtx_sel_Str+"_split"+num+".root","recreate");

   TTree *tree = new TTree("tree","tree");

   double chg, ecal1x1, ecal3x3,ecal5x5, ecal15x15,ecal25x25, weight, ISOR04, p, energy, pt, eta, eEcal_o, phi, pt2, eta2, phi2, e1dx, e2dx, e3dx, e4dx, e5dx, e6dx, e7dx, e1, e2, e3, e4, e5, e6, e7, e1_o, e2_o, e3_o, e4_o, e5_o, e6_o, e7_o;

   double e1r, e12r, e123r, e23r, e1s, e15, e2r, e2s, e25, e3r, e3s, e35, e4r, e4s, e45, e5r, e5s, e55, e6r, e6s, e65, e7r, e7s, e75;

   double e1r_o, e12r_o, e123r_o, e23r_o, e1s_o, e15_o, e2r_o, e2s_o, e25_o, e3r_o, e3s_o, e35_o, e4r_o, e4s_o, e45_o, e5r_o, e5s_o, e55_o, e6r_o, e6s_o, e65_o, e7r_o, e7s_o, e75_o;

   int Ecal_label_o, index1, index2, run, event, n_nvtx_good, HCal_cellHot_o, HCal_cellHot, tight1, tight2, Trkmatch, HCal_ieta, HCal_iphi, HCal_ieta_o, HCal_iphi_o;

   e1s = 0, e2s = 0, e3s = 0, e4s = 0, e5s = 0, e6s = 0, e7s = 0;
   tree->Branch("run", &run, "run/I");
   tree->Branch("nvtx_good", &n_nvtx_good, "nvtx_good/I");
   tree->Branch("event", &event, "event/I");
   tree->Branch("weight", &weight, "weight/D"); 
   tree->Branch("tight1", &tight1, "tight1/I");
   tree->Branch("HCal_cellHot", &HCal_cellHot, "HCal_cellHot/I");
   tree->Branch("HCal_cellHot_o", &HCal_cellHot_o, "HCal_cellHot_o/I"); 
   tree->Branch("Trkmatch", &Trkmatch, "Trkmatch/I"); //front and back proprogated track is matched.
   tree->Branch("ecal1x1", &ecal1x1, "ecal1x1/D");
   tree->Branch("ecal3x3", &ecal3x3, "ecal3x3/D");
   tree->Branch("ecal5x5", &ecal5x5, "ecal5x5/D");
   tree->Branch("ecal15x15", &ecal15x15, "ecal15x15/D");
   tree->Branch("ecal25x25", &ecal25x25, "ecal25x25/D");
   tree->Branch("eEcal_o", &eEcal_o, "eEcal_o/D");
   tree->Branch("Ecal_label_o", &Ecal_label_o, "Ecal_label_o/D");
   tree->Branch("ISOR04", &ISOR04, "ISOR04/D");
   tree->Branch("HCal_ieta", &HCal_ieta, "HCal_ieta/I");
   tree->Branch("HCal_iphi", &HCal_iphi, "HCal_iphi/I");
   tree->Branch("HCal_ieta_o", &HCal_ieta_o, "HCal_ieta_o/I");
   tree->Branch("HCal_iphi_o", &HCal_iphi_o, "HCal_iphi_o/I");
   tree->Branch("eta", &eta, "eta/D");
   tree->Branch("pt", &pt, "pt/D");
   tree->Branch("energy", &energy, "energy/D");
   tree->Branch("phi", &phi, "phi/D");
   tree->Branch("charge", &chg, "charge/D");
   tree->Branch("p", &p, "p/D");

   tree->Branch("e1", &e1, "e1/D");
   tree->Branch("e2", &e2, "e2/D");
   tree->Branch("e3", &e3, "e3/D");
   tree->Branch("e4", &e4, "e4/D");
   tree->Branch("e5", &e5, "e5/D");
   tree->Branch("e6", &e6, "e6/D");
   tree->Branch("e7", &e7, "e7/D");
    
   tree->Branch("e1r", &e1r, "e1r/D");
   tree->Branch("e12r", &e12r, "e12r/D");
   tree->Branch("e123r", &e123r, "e123r/D");
   tree->Branch("e2r", &e2r, "e2r/D");
   tree->Branch("e3r", &e3r, "e3r/D");
   tree->Branch("e4r", &e4r, "e4r/D");
   tree->Branch("e5r", &e5r, "e5r/D");
   tree->Branch("e6r", &e6r, "e6r/D");
   tree->Branch("e7r", &e7r, "e7r/D");
   tree->Branch("e1s", &e1s, "e1s/D");
   tree->Branch("e2s", &e2s, "e2s/D");
   tree->Branch("e3s", &e3s, "e3s/D");
   tree->Branch("e4s", &e4s, "e4s/D");
   tree->Branch("e5s", &e5s, "e5s/D");
   tree->Branch("e6s", &e6s, "e6s/D");
   tree->Branch("e7s", &e7s, "e7s/D");

   tree->Branch("e1_o", &e1_o, "e1_o/D");
   tree->Branch("e2_o", &e2_o, "e2_o/D");
   tree->Branch("e3_o", &e3_o, "e3_o/D");
   tree->Branch("e4_o", &e4_o, "e4_o/D");
   tree->Branch("e5_o", &e5_o, "e5_o/D");
   tree->Branch("e6_o", &e6_o, "e6_o/D");
   tree->Branch("e7_o", &e7_o, "e7_o/D");

   tree->Branch("e1r_o", &e1r_o, "e1r_o/D");
   tree->Branch("e12r_o", &e12r_o, "e12r_o/D");
   tree->Branch("e123r_o", &e123r_o, "e123r_o/D");
   tree->Branch("e2r_o", &e2r_o, "e2r_o/D");
   tree->Branch("e3r_o", &e3r_o, "e3r_o/D");
   tree->Branch("e4r_o", &e4r_o, "e4r_o/D");
   tree->Branch("e5r_o", &e5r_o, "e5r_o/D");
   tree->Branch("e6r_o", &e6r_o, "e6r_o/D");
   tree->Branch("e7r_o", &e7r_o, "e7r_o/D");
   tree->Branch("e1s_o", &e1s_o, "e1s_o/D");
   tree->Branch("e2s_o", &e2s_o, "e2s_o/D");
   tree->Branch("e3s_o", &e3s_o, "e3s_o/D");
   tree->Branch("e4s_o", &e4s_o, "e4s_o/D");
   tree->Branch("e5s_o", &e5s_o, "e5s_o/D");
   tree->Branch("e6s_o", &e6s_o, "e6s_o/D");
   tree->Branch("e7s_o", &e7s_o, "e7s_o/D");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //Long64_t nentries = 500;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%5000==0) std::cout <<"events: "<<jentry<<std::endl;

      int nmuon = IsolationR04->size();

      for(int k = 0; k<nmuon; k ++){
      
           if(nvtx_sel==-1){
              if((!(hcal_cellHot->at(k))) || (!(matchedId->at(k))) || (MuonIsTight->at(k)==0) || (IsolationR04->at(k)>=0.15) || (!(hcal_cellHot_o->at(k))) || (pt_of_muon->at(k))<20 ) continue;
           }
           else{
               if((!(hcal_cellHot->at(k))) || (!(matchedId->at(k))) || (GoodVertex!= nvtx_sel) || (MuonIsTight->at(k)==0) || (IsolationR04->at(k)>=0.15) || (!(hcal_cellHot_o->at(k))) || (pt_of_muon->at(k))<20 ) continue;
           }
           bool match1 = hcal_depthMatch1->at(k);
           bool match2 = hcal_depthMatch2->at(k);
           bool match3 = hcal_depthMatch3->at(k);
           bool match4 = hcal_depthMatch4->at(k);
           bool match5 = hcal_depthMatch5->at(k);
           bool match6 = hcal_depthMatch6->at(k);
           bool match7 = hcal_depthMatch7->at(k);
           if(!(match1 && match2 && match3 && match4 && match5 && match6 && match7)) continue;

           run = Run_No;
           event = Event_No;
           n_nvtx_good = GoodVertex;
           weight = 1.0;
           HCal_cellHot = hcal_cellHot->at(k);
           HCal_cellHot_o = hcal_cellHot_o->at(k);
           tight1 = MuonIsTight->at(k); //istight->at(k);
           ISOR04 = IsolationR04->at(k);
           pt = pt_of_muon->at(k); 
           phi = phi_of_muon->at(k); 
           eta = eta_of_muon->at(k);
           chg = charge_of_muon->at(k);
           energy = energy_of_muon->at(k);
           p = p_of_muon->at(k);
           ecal1x1 = ecal_1x1->at(k);
           ecal3x3 = ecal_3x3->at(k);
           ecal5x5 = ecal_5x5->at(k);
           ecal15x15 = ecal_15x15->at(k);
           ecal25x25 = ecal_25x25->at(k);
           eEcal_o  = eEcal->at(k);
           Ecal_label_o  = Ecal_label->at(k);
           HCal_iphi = hcal_iphi->at(k);
           HCal_ieta = hcal_ieta->at(k);
           HCal_iphi_o = hcal_iphi_o->at(k);
           HCal_ieta_o = hcal_ieta_o->at(k);

           e1 = hcal_edepth1->at(k);
           e2 = hcal_edepth2->at(k);
           e3 = hcal_edepth3->at(k);
           e4 = hcal_edepth4->at(k);
           e5 = hcal_edepth5->at(k);
           e6 = hcal_edepth6->at(k);
           e7 = hcal_edepth7->at(k);

           e1s  = (hcal_edepth11 ->at(k) + hcal_edepth21 ->at(k) + hcal_edepth31 ->at(k) + hcal_edepth41 ->at(k) + hcal_edepth51 ->at(k) + hcal_edepth61 ->at(k) + hcal_edepth71 ->at(k) + hcal_edepth81 ->at(k));
           e2s  = (hcal_edepth12 ->at(k) + hcal_edepth22 ->at(k) + hcal_edepth32 ->at(k) + hcal_edepth42 ->at(k) + hcal_edepth52 ->at(k) + hcal_edepth62 ->at(k) + hcal_edepth72 ->at(k) + hcal_edepth82 ->at(k));
           e3s  = (hcal_edepth13 ->at(k) + hcal_edepth23 ->at(k) + hcal_edepth33 ->at(k) + hcal_edepth43 ->at(k) + hcal_edepth53 ->at(k) + hcal_edepth63 ->at(k) + hcal_edepth73 ->at(k) + hcal_edepth83 ->at(k));
           e4s  = (hcal_edepth14 ->at(k) + hcal_edepth24 ->at(k) + hcal_edepth34 ->at(k) + hcal_edepth44 ->at(k) + hcal_edepth54 ->at(k) + hcal_edepth64 ->at(k) + hcal_edepth74 ->at(k) + hcal_edepth84 ->at(k));
           e5s  = (hcal_edepth15 ->at(k) + hcal_edepth25 ->at(k) + hcal_edepth35 ->at(k) + hcal_edepth45 ->at(k) + hcal_edepth55 ->at(k) + hcal_edepth65 ->at(k) + hcal_edepth75 ->at(k) + hcal_edepth85 ->at(k));
           e6s  = (hcal_edepth16 ->at(k) + hcal_edepth26 ->at(k) + hcal_edepth36 ->at(k) + hcal_edepth46 ->at(k) + hcal_edepth56 ->at(k) + hcal_edepth66 ->at(k) + hcal_edepth76 ->at(k) + hcal_edepth86 ->at(k));
           e7s  = (hcal_edepth17 ->at(k) + hcal_edepth27 ->at(k) + hcal_edepth37 ->at(k) + hcal_edepth47 ->at(k) + hcal_edepth57 ->at(k) + hcal_edepth67 ->at(k) + hcal_edepth77 ->at(k) + hcal_edepth87 ->at(k));

           e1r  = e1s/e1 ;
           e12r  =( e1s+e2s)/e1 ;
           e123r  =( e1s+e2s+e3s)/e1 ;

           e2r  = e2s/e2 ;
           e3r  = e3s/e3 ;
           e4r  = e4s/e4 ;
           e5r  = e5s/e5 ;
           e6r  = e6s/e6 ;
           e7r  = e7s/e7 ;

           e1_o = hcal_edepth1_o->at(k);
           e2_o = hcal_edepth2_o->at(k);
           e3_o = hcal_edepth3_o->at(k);
           e4_o = hcal_edepth4_o->at(k);
           e5_o = hcal_edepth5_o->at(k);
           e6_o = hcal_edepth6_o->at(k);
           e7_o = hcal_edepth7_o->at(k);

           e1s_o = (hcal_edepth11_o->at(k) + hcal_edepth21_o->at(k) + hcal_edepth31_o->at(k) + hcal_edepth41_o->at(k) + hcal_edepth51_o->at(k) + hcal_edepth61_o->at(k) + hcal_edepth71_o->at(k) + hcal_edepth81_o->at(k));
           e2s_o = (hcal_edepth12_o->at(k) + hcal_edepth22_o->at(k) + hcal_edepth32_o->at(k) + hcal_edepth42_o->at(k) + hcal_edepth52_o->at(k) + hcal_edepth62_o->at(k) + hcal_edepth72_o->at(k) + hcal_edepth82_o->at(k));
           e3s_o = (hcal_edepth13_o->at(k) + hcal_edepth23_o->at(k) + hcal_edepth33_o->at(k) + hcal_edepth43_o->at(k) + hcal_edepth53_o->at(k) + hcal_edepth63_o->at(k) + hcal_edepth73_o->at(k) + hcal_edepth83_o->at(k));
           e4s_o = (hcal_edepth14_o->at(k) + hcal_edepth24_o->at(k) + hcal_edepth34_o->at(k) + hcal_edepth44_o->at(k) + hcal_edepth54_o->at(k) + hcal_edepth64_o->at(k) + hcal_edepth74_o->at(k) + hcal_edepth84_o->at(k));
           e5s_o = (hcal_edepth15_o->at(k) + hcal_edepth25_o->at(k) + hcal_edepth35_o->at(k) + hcal_edepth45_o->at(k) + hcal_edepth55_o->at(k) + hcal_edepth65_o->at(k) + hcal_edepth75_o->at(k) + hcal_edepth85_o->at(k));
           e6s_o = (hcal_edepth16_o->at(k) + hcal_edepth26_o->at(k) + hcal_edepth36_o->at(k) + hcal_edepth46_o->at(k) + hcal_edepth56_o->at(k) + hcal_edepth66_o->at(k) + hcal_edepth76_o->at(k) + hcal_edepth86_o->at(k));
           e7s_o = (hcal_edepth17_o->at(k) + hcal_edepth27_o->at(k) + hcal_edepth37_o->at(k) + hcal_edepth47_o->at(k) + hcal_edepth57_o->at(k) + hcal_edepth67_o->at(k) + hcal_edepth77_o->at(k) + hcal_edepth87_o->at(k));

           e1r_o = e1s_o/e1_o;
           e12r_o =( e1s_o+e2s_o)/e1_o;
           e123r_o =( e1s_o+e2s_o+e3s_o)/e1_o;

           e2r_o = e2s_o/e2_o;
           e3r_o = e3s_o/e3_o;
           e4r_o = e4s_o/e4_o;
           e5r_o = e5s_o/e5_o;
           e6r_o = e6s_o/e6_o;
           e7r_o = e7s_o/e7_o;

           tree->Fill();
     }
   }
   outputFile->cd();
   tree->Write();
   outputFile->Close();
}
