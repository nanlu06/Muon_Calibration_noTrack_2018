#define NtuplerReader_cxx
#include "NtuplerReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>

using namespace std;

void NtuplerReader::Loop( TString path, TString prefix, int nvtx_sel )
{

   TString nvtx_sel_Str;
   nvtx_sel_Str.Form("%d",nvtx_sel);
   TFile* outputFile=TFile::Open(path+"/"+prefix+"_nvtx_sel"+nvtx_sel_Str+".root","recreate");

   TTree *tree = new TTree("tree","tree");

   double chg= -999., ecal1x1= -999., ecal3x3= -999.,ecal5x5= -999., ecal15x15= -999.,ecal25x25= -999., weight= -999., ISOR04= -999., p= -999., energy= -999., pt= -999., eta= -999., eEcal_o= -999., phi= -999., pt2= -999., eta2= -999., phi2= -999., e1dx= -999., e2dx= -999., e3dx= -999., e4dx= -999., e5dx= -999., e6dx= -999., e7dx= -999., e1= -999., e2= -999., e3= -999., e4= -999., e5= -999., e6= -999., e7= -999., e1_o = -999., e2_o = -999., e3_o = -999. , e4_o = -999., e5_o = -999., e6_o= -999. , e7_o = -999. ;

   vector<double> pt_genMuon,eta_genMuon,phi_genMuon,energy_genMuon;
   vector<double> pt_probe,eta_probe,phi_probe,energy_probe,charge_probe,p_probe,ISOR04_probe;
   vector<bool>MuonIsMed_probe;
   int nMuon_probe=-999;
   double e1r = -999., e12r= -999., e123r= -999., e23r= -999., e1s= -999., e15= -999., e2r= -999., e2s= -999., e25= -999., e3r= -999., e3s= -999., e35= -999., e4r= -999., e4s= -999., e45= -999., e5r= -999., e5s= -999., e55= -999., e6r= -999., e6s= -999., e65= -999., e7r= -999., e7s= -999., e75= -999.;

   double e1r_o = -999., e12r_o= -999., e123r_o= -999., e23r_o= -999., e1s_o= -999., e15_o= -999., e2r_o= -999., e2s_o= -999., e25_o= -999., e3r_o= -999., e3s_o= -999., e35_o= -999., e4r_o= -999., e4s_o= -999., e45_o= -999., e5r_o= -999., e5s_o= -999., e55_o= -999., e6r_o= -999., e6s_o= -999., e65_o= -999., e7r_o= -999., e7s_o= -999., e75_o= -999.;

   int Ecal_label_o = -999 , index1= -999, index2 = -999, run = -999, event = -999, n_nvtx_good = -999, HCal_cellHot_o = -999, HCal_cellHot = -999, tight1 = -999, tight2 = -999, Trkmatch = -999, HCal_ieta = -999, HCal_iphi = -999, HCal_ieta_o = -999, HCal_iphi_o = -999;

   bool IsMuonRec = false;
   
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
   tree->Branch("eta_genMuon", &eta_genMuon);
   tree->Branch("pt_genMuon", &pt_genMuon);
   tree->Branch("energy_genMuon", &energy_genMuon);
   tree->Branch("phi_genMuon", &phi_genMuon);
   tree->Branch("IsMuonRec", &IsMuonRec);
   tree->Branch("eta_probe", &eta_probe);
   tree->Branch("pt_probe", &pt_probe);
   tree->Branch("energy_probe", &energy_probe);
   tree->Branch("phi_probe", &phi_probe);
   tree->Branch("charge_probe", &charge_probe);
   tree->Branch("p_probe", &p_probe);
   tree->Branch("ISOR04_probe", &ISOR04_probe);
   tree->Branch("MuonIsMed_probe", &MuonIsMed_probe);
   tree->Branch("nMuon_probe",&nMuon_probe);
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

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //if(Event_No!=25221) continue;
      //else std::cout <<"find event==25221 "<<std::endl;
      pt_probe.clear();
      eta_probe.clear();
      phi_probe.clear();
      energy_probe.clear();
      charge_probe.clear();
      p_probe.clear();
      ISOR04_probe.clear();
      MuonIsMed_probe.clear();
      pt_genMuon.clear();
      phi_genMuon.clear();
      eta_genMuon.clear();
      energy_genMuon.clear();
      if(jentry%5000==0) std::cout <<"events: "<<jentry<<std::endl;
      int nmuon = IsolationR04->size();
      int nGenMuon=pt_of_genMuon->size();
      for(int k = 0; k<nGenMuon; k ++){
	pt_genMuon.push_back(pt_of_genMuon->at(k));
	phi_genMuon.push_back(phi_of_genMuon->at(k));
	eta_genMuon.push_back(eta_of_genMuon->at(k));
	energy_genMuon.push_back(energy_of_genMuon->at(k));
      }
      //cout<<jentry<<" : "<<nmuon<<endl;
      
      
      for(int k = 0; k<nmuon; k ++){
	if(nvtx_sel==-1){
	  if((!(hcal_cellHot->at(k))) || (!(matchedId->at(k))) || (MuonIsTight->at(k)==0) || (IsolationR04->at(k)>=0.15) || /*(!(hcal_cellHot_o->at(k))) ||*/ (pt_of_muon->at(k))<20 ) continue;
	}
	else{
	  if((!(hcal_cellHot->at(k))) || (!(matchedId->at(k))) || (GoodVertex!= nvtx_sel) || (MuonIsTight->at(k)==0) || (IsolationR04->at(k)>=0.15) || /*(!(hcal_cellHot_o->at(k))) ||*/ (pt_of_muon->at(k))<20 ) continue;
	}
	bool match1 = hcal_depthMatch1->at(k);
	bool match2 = hcal_depthMatch2->at(k);
           bool match3 = hcal_depthMatch3->at(k);
           bool match4 = hcal_depthMatch4->at(k);
           bool match5 = hcal_depthMatch5->at(k);
           bool match6 = hcal_depthMatch6->at(k);
           bool match7 = hcal_depthMatch7->at(k);
           if(!(match1 && match2 && match3 && match4 && match5 && match6 && match7)) continue;

	   //for(int j = 0; j<nGenMuon; j ++){
	   //if((pt_of_genMuon-pt_of_muon
	   //}
           run = Run_No;
	   
           event = Event_No;
           //cout<<event<<endl;
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
	   IsMuonRec = isMuonRec->at(k);
	   nMuon_probe=n_muon_probe->at(k);
	   double sum_probe=0.;
	   for(int s =0; s<k;s++) sum_probe+=n_muon_probe->at(s);
	   for(int pr=sum_probe;pr<sum_probe+n_muon_probe->at(k);pr++){
	     if(pt_of_muon_probe->at(pr)!=pt && isolationR04_probe->at(pr)<0.25){
	       //cout<<Event_No<<" "<<pt_of_muon_probe->at(pr)<<" "<<pt<<endl;
	       pt_probe.push_back(pt_of_muon_probe->at(pr));
	       phi_probe.push_back(phi_of_muon_probe->at(pr));
	       eta_probe.push_back(eta_of_muon_probe->at(pr));
	       energy_probe.push_back(energy_of_muon_probe->at(pr));
	       charge_probe.push_back(charge_of_muon_probe->at(pr));
	       p_probe.push_back(p_of_muon_probe->at(pr));
	       ISOR04_probe.push_back(isolationR04_probe->at(pr));
	       MuonIsMed_probe.push_back(MuonIsMedium_probe->at(pr));
	     }
	   }
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
	   //cout<<jentry<<": Filling event "<<event<<endl;
           tree->Fill();
	   //cout<<"Finished\n";
      }
   }
   outputFile->cd();
   tree->Write();
   outputFile->Close();
}
