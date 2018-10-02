//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Sep  1 18:56:00 2018 by ROOT version 5.34/38
// from TTree TREE/TREE
// found on file: Validation_209.root
//////////////////////////////////////////////////////////

#ifndef NtuplerReader_h
#define NtuplerReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtuplerReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          Event_No;
   UInt_t          Run_No;
   UInt_t          LumiNumber;
   UInt_t          BXNumber;
   UInt_t          GoodVertex;
   Float_t         LumiScaler;
   Float_t         BunchLumi;
   Float_t         pileup;
   Float_t         pileupRMS;
   vector<bool>    *PF_Muon;
   vector<bool>    *Global_Muon;
   vector<bool>    *Tracker_muon;
   vector<bool>    *MuonIsTight;
   vector<bool>    *MuonIsMedium;
   vector<double>  *charge_of_muon;
   vector<double>  *pt_of_muon;
   vector<double>  *eta_of_muon;
   vector<double>  *phi_of_muon;
   vector<double>  *energy_of_muon;
   vector<double>  *p_of_muon;
   vector<float>   *muon_trkKink;
   vector<float>   *muon_chi2LocalPosition;
   vector<float>   *muon_segComp;
   vector<int>     *TrackerLayer;
   vector<int>     *NumPixelLayers;
   vector<int>     *InnerTrackPixelHits;
   vector<bool>    *innerTrack;
   vector<double>  *chiTracker;
   vector<double>  *DxyTracker;
   vector<double>  *DzTracker;
   vector<double>  *innerTrackpt;
   vector<double>  *innerTracketa;
   vector<double>  *innerTrackphi;
   vector<double>  *tight_validFraction;
   vector<bool>    *OuterTrack;
   vector<double>  *OuterTrackChi;
   vector<double>  *OuterTrackPt;
   vector<double>  *OuterTrackEta;
   vector<double>  *OuterTrackPhi;
   vector<int>     *OuterTrackHits;
   vector<int>     *OuterTrackRHits;
   vector<bool>    *GlobalTrack;
   vector<double>  *GlobalTrckPt;
   vector<double>  *GlobalTrckEta;
   vector<double>  *GlobalTrckPhi;
   vector<int>     *Global_Muon_Hits;
   vector<int>     *MatchedStations;
   vector<double>  *GlobTrack_Chi;
   vector<double>  *Tight_LongitudinalImpactparameter;
   vector<double>  *Tight_TransImpactparameter;
   vector<double>  *IsolationR04;
   vector<double>  *IsolationR03;
   vector<double>  *ecal_3into3;
   vector<double>  *hcal_3into3;
   vector<double>  *tracker_3into3;
   vector<bool>    *matchedId;
   vector<bool>    *hcal_cellHot;
   vector<bool>    *hcal_cellHot_o;
   vector<double>  *pt_of_genMuon;
   vector<double>  *eta_of_genMuon;
   vector<double>  *phi_of_genMuon;
   vector<double>  *energy_of_genMuon;
   vector<double>  *ecal_1x1;
   vector<double>  *ecal_3x3;
   vector<double>  *ecal_5x5;
   vector<double>  *ecal_15x15;
   vector<double>  *ecal_25x25;
   vector<double>  *hcal_1x1;
   vector<unsigned int> *ecal_detID;
   vector<unsigned int> *hcal_detID;
   vector<unsigned int> *ehcal_detID;
   vector<int>     *hcal_ieta;
   vector<int>     *hcal_iphi;
   vector<double>  *hcal_edepthuncalibrated1;
   vector<double>  *hcal_ndepth1;
   vector<double>  *hcal_edepth11;
   vector<double>  *hcal_ndepth11;
   vector<double>  *hcal_edepth21;
   vector<double>  *hcal_ndepth21;
   vector<double>  *hcal_edepth31;
   vector<double>  *hcal_ndepth31;
   vector<double>  *hcal_edepth41;
   vector<double>  *hcal_ndepth41;
   vector<double>  *hcal_edepth51;
   vector<double>  *hcal_ndepth51;
   vector<double>  *hcal_edepth61;
   vector<double>  *hcal_ndepth61;
   vector<double>  *hcal_edepth71;
   vector<double>  *hcal_ndepth71;
   vector<double>  *hcal_edepth81;
   vector<double>  *hcal_ndepth81;
   vector<double>  *hcal_activeL1;
   vector<double>  *hcal_edepth1;
   vector<bool>    *hcal_depthMatch1;
   vector<double>  *hcal_edepthuncalibrated2;
   vector<double>  *hcal_ndepth2;
   vector<double>  *hcal_edepth12;
   vector<double>  *hcal_ndepth12;
   vector<double>  *hcal_edepth22;
   vector<double>  *hcal_ndepth22;
   vector<double>  *hcal_edepth32;
   vector<double>  *hcal_ndepth32;
   vector<double>  *hcal_edepth42;
   vector<double>  *hcal_ndepth42;
   vector<double>  *hcal_edepth52;
   vector<double>  *hcal_ndepth52;
   vector<double>  *hcal_edepth62;
   vector<double>  *hcal_ndepth62;
   vector<double>  *hcal_edepth72;
   vector<double>  *hcal_ndepth72;
   vector<double>  *hcal_edepth82;
   vector<double>  *hcal_ndepth82;
   vector<double>  *hcal_activeL2;
   vector<double>  *hcal_edepth2;
   vector<bool>    *hcal_depthMatch2;
   vector<double>  *hcal_edepthuncalibrated3;
   vector<double>  *hcal_ndepth3;
   vector<double>  *hcal_edepth13;
   vector<double>  *hcal_ndepth13;
   vector<double>  *hcal_edepth23;
   vector<double>  *hcal_ndepth23;
   vector<double>  *hcal_edepth33;
   vector<double>  *hcal_ndepth33;
   vector<double>  *hcal_edepth43;
   vector<double>  *hcal_ndepth43;
   vector<double>  *hcal_edepth53;
   vector<double>  *hcal_ndepth53;
   vector<double>  *hcal_edepth63;
   vector<double>  *hcal_ndepth63;
   vector<double>  *hcal_edepth73;
   vector<double>  *hcal_ndepth73;
   vector<double>  *hcal_edepth83;
   vector<double>  *hcal_ndepth83;
   vector<double>  *hcal_activeL3;
   vector<double>  *hcal_edepth3;
   vector<bool>    *hcal_depthMatch3;
   vector<double>  *hcal_edepthuncalibrated4;
   vector<double>  *hcal_ndepth4;
   vector<double>  *hcal_edepth14;
   vector<double>  *hcal_ndepth14;
   vector<double>  *hcal_edepth24;
   vector<double>  *hcal_ndepth24;
   vector<double>  *hcal_edepth34;
   vector<double>  *hcal_ndepth34;
   vector<double>  *hcal_edepth44;
   vector<double>  *hcal_ndepth44;
   vector<double>  *hcal_edepth54;
   vector<double>  *hcal_ndepth54;
   vector<double>  *hcal_edepth64;
   vector<double>  *hcal_ndepth64;
   vector<double>  *hcal_edepth74;
   vector<double>  *hcal_ndepth74;
   vector<double>  *hcal_edepth84;
   vector<double>  *hcal_ndepth84;
   vector<double>  *hcal_activeL4;
   vector<double>  *hcal_edepth4;
   vector<bool>    *hcal_depthMatch4;
   vector<double>  *hcal_edepthuncalibrated5;
   vector<double>  *hcal_ndepth5;
   vector<double>  *hcal_edepth15;
   vector<double>  *hcal_ndepth15;
   vector<double>  *hcal_edepth25;
   vector<double>  *hcal_ndepth25;
   vector<double>  *hcal_edepth35;
   vector<double>  *hcal_ndepth35;
   vector<double>  *hcal_edepth45;
   vector<double>  *hcal_ndepth45;
   vector<double>  *hcal_edepth55;
   vector<double>  *hcal_ndepth55;
   vector<double>  *hcal_edepth65;
   vector<double>  *hcal_ndepth65;
   vector<double>  *hcal_edepth75;
   vector<double>  *hcal_ndepth75;
   vector<double>  *hcal_edepth85;
   vector<double>  *hcal_ndepth85;
   vector<double>  *hcal_activeL5;
   vector<double>  *hcal_edepth5;
   vector<bool>    *hcal_depthMatch5;
   vector<double>  *hcal_edepthuncalibrated6;
   vector<double>  *hcal_ndepth6;
   vector<double>  *hcal_edepth16;
   vector<double>  *hcal_ndepth16;
   vector<double>  *hcal_edepth26;
   vector<double>  *hcal_ndepth26;
   vector<double>  *hcal_edepth36;
   vector<double>  *hcal_ndepth36;
   vector<double>  *hcal_edepth46;
   vector<double>  *hcal_ndepth46;
   vector<double>  *hcal_edepth56;
   vector<double>  *hcal_ndepth56;
   vector<double>  *hcal_edepth66;
   vector<double>  *hcal_ndepth66;
   vector<double>  *hcal_edepth76;
   vector<double>  *hcal_ndepth76;
   vector<double>  *hcal_edepth86;
   vector<double>  *hcal_ndepth86;
   vector<double>  *hcal_activeL6;
   vector<double>  *hcal_edepth6;
   vector<bool>    *hcal_depthMatch6;
   vector<double>  *hcal_edepthuncalibrated7;
   vector<double>  *hcal_ndepth7;
   vector<double>  *hcal_edepth17;
   vector<double>  *hcal_ndepth17;
   vector<double>  *hcal_edepth27;
   vector<double>  *hcal_ndepth27;
   vector<double>  *hcal_edepth37;
   vector<double>  *hcal_ndepth37;
   vector<double>  *hcal_edepth47;
   vector<double>  *hcal_ndepth47;
   vector<double>  *hcal_edepth57;
   vector<double>  *hcal_ndepth57;
   vector<double>  *hcal_edepth67;
   vector<double>  *hcal_ndepth67;
   vector<double>  *hcal_edepth77;
   vector<double>  *hcal_ndepth77;
   vector<double>  *hcal_edepth87;
   vector<double>  *hcal_ndepth87;
   vector<double>  *hcal_activeL7;
   vector<double>  *hcal_edepth7;
   vector<bool>    *hcal_depthMatch7;
   vector<double>  *activeLength;
   vector<double>  *eEcal;
   vector<double>  *Ecal_label;
   vector<int>     *hcal_ieta_o;
   vector<int>     *hcal_iphi_o;
   vector<double>  *hcal_edepth1_o;
   vector<double>  *hcal_edepth2_o;
   vector<double>  *hcal_edepth3_o;
   vector<double>  *hcal_edepth4_o;
   vector<double>  *hcal_edepth11_o;
   vector<double>  *hcal_edepth12_o;
   vector<double>  *hcal_edepth13_o;
   vector<double>  *hcal_edepth14_o;
   vector<double>  *hcal_edepth21_o;
   vector<double>  *hcal_edepth22_o;
   vector<double>  *hcal_edepth23_o;
   vector<double>  *hcal_edepth24_o;
   vector<double>  *hcal_edepth31_o;
   vector<double>  *hcal_edepth32_o;
   vector<double>  *hcal_edepth33_o;
   vector<double>  *hcal_edepth34_o;
   vector<double>  *hcal_edepth41_o;
   vector<double>  *hcal_edepth42_o;
   vector<double>  *hcal_edepth43_o;
   vector<double>  *hcal_edepth44_o;
   vector<double>  *hcal_edepth51_o;
   vector<double>  *hcal_edepth52_o;
   vector<double>  *hcal_edepth53_o;
   vector<double>  *hcal_edepth54_o;
   vector<double>  *hcal_edepth61_o;
   vector<double>  *hcal_edepth62_o;
   vector<double>  *hcal_edepth63_o;
   vector<double>  *hcal_edepth64_o;
   vector<double>  *hcal_edepth71_o;
   vector<double>  *hcal_edepth72_o;
   vector<double>  *hcal_edepth73_o;
   vector<double>  *hcal_edepth74_o;
   vector<double>  *hcal_edepth81_o;
   vector<double>  *hcal_edepth82_o;
   vector<double>  *hcal_edepth83_o;
   vector<double>  *hcal_edepth84_o;
   vector<double>  *hcal_edepth5_o;
   vector<double>  *hcal_edepth15_o;
   vector<double>  *hcal_edepth25_o;
   vector<double>  *hcal_edepth35_o;
   vector<double>  *hcal_edepth45_o;
   vector<double>  *hcal_edepth55_o;
   vector<double>  *hcal_edepth65_o;
   vector<double>  *hcal_edepth75_o;
   vector<double>  *hcal_edepth85_o;
   vector<double>  *hcal_edepth6_o;
   vector<double>  *hcal_edepth16_o;
   vector<double>  *hcal_edepth26_o;
   vector<double>  *hcal_edepth36_o;
   vector<double>  *hcal_edepth46_o;
   vector<double>  *hcal_edepth56_o;
   vector<double>  *hcal_edepth66_o;
   vector<double>  *hcal_edepth76_o;
   vector<double>  *hcal_edepth86_o;
   vector<double>  *hcal_edepth7_o;
   vector<double>  *hcal_edepth17_o;
   vector<double>  *hcal_edepth27_o;
   vector<double>  *hcal_edepth37_o;
   vector<double>  *hcal_edepth47_o;
   vector<double>  *hcal_edepth57_o;
   vector<double>  *hcal_edepth67_o;
   vector<double>  *hcal_edepth77_o;
   vector<double>  *hcal_edepth87_o;
   vector<int>     *hltresults;
   vector<string>  *all_triggers;

   // List of branches
   TBranch        *b_Event_No;   //!
   TBranch        *b_Run_No;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_BXNumber;   //!
   TBranch        *b_GoodVertex;   //!
   TBranch        *b_LumiScaler;   //!
   TBranch        *b_BunchLumi;   //!
   TBranch        *b_pileup;   //!
   TBranch        *b_pileupRMS;   //!
   TBranch        *b_PF_Muon;   //!
   TBranch        *b_Global_Muon;   //!
   TBranch        *b_Tracker_muon;   //!
   TBranch        *b_MuonIsTight;   //!
   TBranch        *b_MuonIsMedium;   //!
   TBranch        *b_charge_of_muon;   //!
   TBranch        *b_pt_of_muon;   //!
   TBranch        *b_eta_of_muon;   //!
   TBranch        *b_phi_of_muon;   //!
   TBranch        *b_energy_of_muon;   //!
   TBranch        *b_p_of_muon;   //!
   TBranch        *b_muon_trkKink;   //!
   TBranch        *b_muon_chi2LocalPosition;   //!
   TBranch        *b_muon_segComp;   //!
   TBranch        *b_TrackerLayer;   //!
   TBranch        *b_NumPixelLayers;   //!
   TBranch        *b_InnerTrackPixelHits;   //!
   TBranch        *b_innerTrack;   //!
   TBranch        *b_chiTracker;   //!
   TBranch        *b_DxyTracker;   //!
   TBranch        *b_DzTracker;   //!
   TBranch        *b_innerTrackpt;   //!
   TBranch        *b_innerTracketa;   //!
   TBranch        *b_innerTrackphi;   //!
   TBranch        *b_tight_validFraction;   //!
   TBranch        *b_OuterTrack;   //!
   TBranch        *b_OuterTrackChi;   //!
   TBranch        *b_OuterTrackPt;   //!
   TBranch        *b_OuterTrackEta;   //!
   TBranch        *b_OuterTrackPhi;   //!
   TBranch        *b_OuterTrackHits;   //!
   TBranch        *b_OuterTrackRHits;   //!
   TBranch        *b_GlobalTrack;   //!
   TBranch        *b_GlobalTrckPt;   //!
   TBranch        *b_GlobalTrckEta;   //!
   TBranch        *b_GlobalTrckPhi;   //!
   TBranch        *b_Global_Muon_Hits;   //!
   TBranch        *b_MatchedStations;   //!
   TBranch        *b_GlobTrack_Chi;   //!
   TBranch        *b_Tight_LongitudinalImpactparameter;   //!
   TBranch        *b_Tight_TransImpactparameter;   //!
   TBranch        *b_IsolationR04;   //!
   TBranch        *b_IsolationR03;   //!
   TBranch        *b_ecal_3into3;   //!
   TBranch        *b_hcal_3into3;   //!
   TBranch        *b_tracker_3into3;   //!
   TBranch        *b_matchedId;   //!
   TBranch        *b_hcal_cellHot;   //!
   TBranch        *b_hcal_cellHot_o;   //!
   TBranch        *b_pt_of_genMuon;//!
   TBranch        *b_eta_of_genMuon;//!
   TBranch        *b_phi_of_genMuon;//!
   TBranch        *b_energy_of_genMuon;//!
   TBranch        *b_ecal_1x1;   //!
   TBranch        *b_ecal_3x3;   //!
   TBranch        *b_ecal_5x5;   //!
   TBranch        *b_ecal_15x15;   //!
   TBranch        *b_ecal_25x25;   //!
   TBranch        *b_hcal_1x1;   //!
   TBranch        *b_ecal_detID;   //!
   TBranch        *b_hcal_detID;   //!
   TBranch        *b_ehcal_detID;   //!
   TBranch        *b_hcal_ieta;   //!
   TBranch        *b_hcal_iphi;   //!
   TBranch        *b_hcal_edepthuncalibrated1;   //!
   TBranch        *b_hcal_ndepth1;   //!
   TBranch        *b_hcal_edepth11;   //!
   TBranch        *b_hcal_ndepth11;   //!
   TBranch        *b_hcal_edepth21;   //!
   TBranch        *b_hcal_ndepth21;   //!
   TBranch        *b_hcal_edepth31;   //!
   TBranch        *b_hcal_ndepth31;   //!
   TBranch        *b_hcal_edepth41;   //!
   TBranch        *b_hcal_ndepth41;   //!
   TBranch        *b_hcal_edepth51;   //!
   TBranch        *b_hcal_ndepth51;   //!
   TBranch        *b_hcal_edepth61;   //!
   TBranch        *b_hcal_ndepth61;   //!
   TBranch        *b_hcal_edepth71;   //!
   TBranch        *b_hcal_ndepth71;   //!
   TBranch        *b_hcal_edepth81;   //!
   TBranch        *b_hcal_ndepth81;   //!
   TBranch        *b_hcal_activeL1;   //!
   TBranch        *b_hcal_edepth1;   //!
   TBranch        *b_hcal_depthMatch1;   //!
   TBranch        *b_hcal_edepthuncalibrated2;   //!
   TBranch        *b_hcal_ndepth2;   //!
   TBranch        *b_hcal_edepth12;   //!
   TBranch        *b_hcal_ndepth12;   //!
   TBranch        *b_hcal_edepth22;   //!
   TBranch        *b_hcal_ndepth22;   //!
   TBranch        *b_hcal_edepth32;   //!
   TBranch        *b_hcal_ndepth32;   //!
   TBranch        *b_hcal_edepth42;   //!
   TBranch        *b_hcal_ndepth42;   //!
   TBranch        *b_hcal_edepth52;   //!
   TBranch        *b_hcal_ndepth52;   //!
   TBranch        *b_hcal_edepth62;   //!
   TBranch        *b_hcal_ndepth62;   //!
   TBranch        *b_hcal_edepth72;   //!
   TBranch        *b_hcal_ndepth72;   //!
   TBranch        *b_hcal_edepth82;   //!
   TBranch        *b_hcal_ndepth82;   //!
   TBranch        *b_hcal_activeL2;   //!
   TBranch        *b_hcal_edepth2;   //!
   TBranch        *b_hcal_depthMatch2;   //!
   TBranch        *b_hcal_edepthuncalibrated3;   //!
   TBranch        *b_hcal_ndepth3;   //!
   TBranch        *b_hcal_edepth13;   //!
   TBranch        *b_hcal_ndepth13;   //!
   TBranch        *b_hcal_edepth23;   //!
   TBranch        *b_hcal_ndepth23;   //!
   TBranch        *b_hcal_edepth33;   //!
   TBranch        *b_hcal_ndepth33;   //!
   TBranch        *b_hcal_edepth43;   //!
   TBranch        *b_hcal_ndepth43;   //!
   TBranch        *b_hcal_edepth53;   //!
   TBranch        *b_hcal_ndepth53;   //!
   TBranch        *b_hcal_edepth63;   //!
   TBranch        *b_hcal_ndepth63;   //!
   TBranch        *b_hcal_edepth73;   //!
   TBranch        *b_hcal_ndepth73;   //!
   TBranch        *b_hcal_edepth83;   //!
   TBranch        *b_hcal_ndepth83;   //!
   TBranch        *b_hcal_activeL3;   //!
   TBranch        *b_hcal_edepth3;   //!
   TBranch        *b_hcal_depthMatch3;   //!
   TBranch        *b_hcal_edepthuncalibrated4;   //!
   TBranch        *b_hcal_ndepth4;   //!
   TBranch        *b_hcal_edepth14;   //!
   TBranch        *b_hcal_ndepth14;   //!
   TBranch        *b_hcal_edepth24;   //!
   TBranch        *b_hcal_ndepth24;   //!
   TBranch        *b_hcal_edepth34;   //!
   TBranch        *b_hcal_ndepth34;   //!
   TBranch        *b_hcal_edepth44;   //!
   TBranch        *b_hcal_ndepth44;   //!
   TBranch        *b_hcal_edepth54;   //!
   TBranch        *b_hcal_ndepth54;   //!
   TBranch        *b_hcal_edepth64;   //!
   TBranch        *b_hcal_ndepth64;   //!
   TBranch        *b_hcal_edepth74;   //!
   TBranch        *b_hcal_ndepth74;   //!
   TBranch        *b_hcal_edepth84;   //!
   TBranch        *b_hcal_ndepth84;   //!
   TBranch        *b_hcal_activeL4;   //!
   TBranch        *b_hcal_edepth4;   //!
   TBranch        *b_hcal_depthMatch4;   //!
   TBranch        *b_hcal_edepthuncalibrated5;   //!
   TBranch        *b_hcal_ndepth5;   //!
   TBranch        *b_hcal_edepth15;   //!
   TBranch        *b_hcal_ndepth15;   //!
   TBranch        *b_hcal_edepth25;   //!
   TBranch        *b_hcal_ndepth25;   //!
   TBranch        *b_hcal_edepth35;   //!
   TBranch        *b_hcal_ndepth35;   //!
   TBranch        *b_hcal_edepth45;   //!
   TBranch        *b_hcal_ndepth45;   //!
   TBranch        *b_hcal_edepth55;   //!
   TBranch        *b_hcal_ndepth55;   //!
   TBranch        *b_hcal_edepth65;   //!
   TBranch        *b_hcal_ndepth65;   //!
   TBranch        *b_hcal_edepth75;   //!
   TBranch        *b_hcal_ndepth75;   //!
   TBranch        *b_hcal_edepth85;   //!
   TBranch        *b_hcal_ndepth85;   //!
   TBranch        *b_hcal_activeL5;   //!
   TBranch        *b_hcal_edepth5;   //!
   TBranch        *b_hcal_depthMatch5;   //!
   TBranch        *b_hcal_edepthuncalibrated6;   //!
   TBranch        *b_hcal_ndepth6;   //!
   TBranch        *b_hcal_edepth16;   //!
   TBranch        *b_hcal_ndepth16;   //!
   TBranch        *b_hcal_edepth26;   //!
   TBranch        *b_hcal_ndepth26;   //!
   TBranch        *b_hcal_edepth36;   //!
   TBranch        *b_hcal_ndepth36;   //!
   TBranch        *b_hcal_edepth46;   //!
   TBranch        *b_hcal_ndepth46;   //!
   TBranch        *b_hcal_edepth56;   //!
   TBranch        *b_hcal_ndepth56;   //!
   TBranch        *b_hcal_edepth66;   //!
   TBranch        *b_hcal_ndepth66;   //!
   TBranch        *b_hcal_edepth76;   //!
   TBranch        *b_hcal_ndepth76;   //!
   TBranch        *b_hcal_edepth86;   //!
   TBranch        *b_hcal_ndepth86;   //!
   TBranch        *b_hcal_activeL6;   //!
   TBranch        *b_hcal_edepth6;   //!
   TBranch        *b_hcal_depthMatch6;   //!
   TBranch        *b_hcal_edepthuncalibrated7;   //!
   TBranch        *b_hcal_ndepth7;   //!
   TBranch        *b_hcal_edepth17;   //!
   TBranch        *b_hcal_ndepth17;   //!
   TBranch        *b_hcal_edepth27;   //!
   TBranch        *b_hcal_ndepth27;   //!
   TBranch        *b_hcal_edepth37;   //!
   TBranch        *b_hcal_ndepth37;   //!
   TBranch        *b_hcal_edepth47;   //!
   TBranch        *b_hcal_ndepth47;   //!
   TBranch        *b_hcal_edepth57;   //!
   TBranch        *b_hcal_ndepth57;   //!
   TBranch        *b_hcal_edepth67;   //!
   TBranch        *b_hcal_ndepth67;   //!
   TBranch        *b_hcal_edepth77;   //!
   TBranch        *b_hcal_ndepth77;   //!
   TBranch        *b_hcal_edepth87;   //!
   TBranch        *b_hcal_ndepth87;   //!
   TBranch        *b_hcal_activeL7;   //!
   TBranch        *b_hcal_edepth7;   //!
   TBranch        *b_hcal_depthMatch7;   //!
   TBranch        *b_activeLength;   //!
   TBranch        *b_eEcal;   //!
   TBranch        *b_Ecal_label;   //!
   TBranch        *b_hcal_ieta_o;   //!
   TBranch        *b_hcal_iphi_o;   //!
   TBranch        *b_hcal_edepth1_o;   //!
   TBranch        *b_hcal_edepth2_o;   //!
   TBranch        *b_hcal_edepth3_o;   //!
   TBranch        *b_hcal_edepth4_o;   //!
   TBranch        *b_hcal_edepth11_o;   //!
   TBranch        *b_hcal_edepth12_o;   //!
   TBranch        *b_hcal_edepth13_o;   //!
   TBranch        *b_hcal_edepth14_o;   //!
   TBranch        *b_hcal_edepth21_o;   //!
   TBranch        *b_hcal_edepth22_o;   //!
   TBranch        *b_hcal_edepth23_o;   //!
   TBranch        *b_hcal_edepth24_o;   //!
   TBranch        *b_hcal_edepth31_o;   //!
   TBranch        *b_hcal_edepth32_o;   //!
   TBranch        *b_hcal_edepth33_o;   //!
   TBranch        *b_hcal_edepth34_o;   //!
   TBranch        *b_hcal_edepth41_o;   //!
   TBranch        *b_hcal_edepth42_o;   //!
   TBranch        *b_hcal_edepth43_o;   //!
   TBranch        *b_hcal_edepth44_o;   //!
   TBranch        *b_hcal_edepth51_o;   //!
   TBranch        *b_hcal_edepth52_o;   //!
   TBranch        *b_hcal_edepth53_o;   //!
   TBranch        *b_hcal_edepth54_o;   //!
   TBranch        *b_hcal_edepth61_o;   //!
   TBranch        *b_hcal_edepth62_o;   //!
   TBranch        *b_hcal_edepth63_o;   //!
   TBranch        *b_hcal_edepth64_o;   //!
   TBranch        *b_hcal_edepth71_o;   //!
   TBranch        *b_hcal_edepth72_o;   //!
   TBranch        *b_hcal_edepth73_o;   //!
   TBranch        *b_hcal_edepth74_o;   //!
   TBranch        *b_hcal_edepth81_o;   //!
   TBranch        *b_hcal_edepth82_o;   //!
   TBranch        *b_hcal_edepth83_o;   //!
   TBranch        *b_hcal_edepth84_o;   //!
   TBranch        *b_hcal_edepth5_o;   //!
   TBranch        *b_hcal_edepth15_o;   //!
   TBranch        *b_hcal_edepth25_o;   //!
   TBranch        *b_hcal_edepth35_o;   //!
   TBranch        *b_hcal_edepth45_o;   //!
   TBranch        *b_hcal_edepth55_o;   //!
   TBranch        *b_hcal_edepth65_o;   //!
   TBranch        *b_hcal_edepth75_o;   //!
   TBranch        *b_hcal_edepth85_o;   //!
   TBranch        *b_hcal_edepth6_o;   //!
   TBranch        *b_hcal_edepth16_o;   //!
   TBranch        *b_hcal_edepth26_o;   //!
   TBranch        *b_hcal_edepth36_o;   //!
   TBranch        *b_hcal_edepth46_o;   //!
   TBranch        *b_hcal_edepth56_o;   //!
   TBranch        *b_hcal_edepth66_o;   //!
   TBranch        *b_hcal_edepth76_o;   //!
   TBranch        *b_hcal_edepth86_o;   //!
   TBranch        *b_hcal_edepth7_o;   //!
   TBranch        *b_hcal_edepth17_o;   //!
   TBranch        *b_hcal_edepth27_o;   //!
   TBranch        *b_hcal_edepth37_o;   //!
   TBranch        *b_hcal_edepth47_o;   //!
   TBranch        *b_hcal_edepth57_o;   //!
   TBranch        *b_hcal_edepth67_o;   //!
   TBranch        *b_hcal_edepth77_o;   //!
   TBranch        *b_hcal_edepth87_o;   //!
   TBranch        *b_hltresults;   //!
   TBranch        *b_all_triggers;   //!

   NtuplerReader(TTree *tree=0);
   virtual ~NtuplerReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop( TString path, TString num,  int nvtx_sel);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtuplerReader_cxx
NtuplerReader::NtuplerReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Validation_209.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Validation_209.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Validation_209.root:/hcalHBHEMuon");
      dir->GetObject("TREE",tree);

   }
   Init(tree);
}

NtuplerReader::~NtuplerReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NtuplerReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtuplerReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NtuplerReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PF_Muon = 0;
   Global_Muon = 0;
   Tracker_muon = 0;
   MuonIsTight = 0;
   MuonIsMedium = 0;
   charge_of_muon = 0;
   pt_of_muon = 0;
   eta_of_muon = 0;
   phi_of_muon = 0;
   energy_of_muon = 0;
   p_of_muon = 0;
   muon_trkKink = 0;
   muon_chi2LocalPosition = 0;
   muon_segComp = 0;
   TrackerLayer = 0;
   NumPixelLayers = 0;
   InnerTrackPixelHits = 0;
   innerTrack = 0;
   chiTracker = 0;
   DxyTracker = 0;
   DzTracker = 0;
   innerTrackpt = 0;
   innerTracketa = 0;
   innerTrackphi = 0;
   tight_validFraction = 0;
   OuterTrack = 0;
   OuterTrackChi = 0;
   OuterTrackPt = 0;
   OuterTrackEta = 0;
   OuterTrackPhi = 0;
   OuterTrackHits = 0;
   OuterTrackRHits = 0;
   GlobalTrack = 0;
   GlobalTrckPt = 0;
   GlobalTrckEta = 0;
   GlobalTrckPhi = 0;
   Global_Muon_Hits = 0;
   MatchedStations = 0;
   GlobTrack_Chi = 0;
   Tight_LongitudinalImpactparameter = 0;
   Tight_TransImpactparameter = 0;
   IsolationR04 = 0;
   IsolationR03 = 0;
   ecal_3into3 = 0;
   hcal_3into3 = 0;
   tracker_3into3 = 0;
   matchedId = 0;
   hcal_cellHot = 0;
   hcal_cellHot_o = 0;
   pt_of_genMuon = 0;//!
   eta_of_genMuon = 0;//!
   phi_of_genMuon = 0;//!
   energy_of_genMuon = 0;//!
   ecal_1x1 = 0;
   ecal_3x3 = 0;
   ecal_5x5 = 0;
   ecal_15x15 = 0;
   ecal_25x25 = 0;
   hcal_1x1 = 0;
   ecal_detID = 0;
   hcal_detID = 0;
   ehcal_detID = 0;
   hcal_ieta = 0;
   hcal_iphi = 0;
   hcal_edepthuncalibrated1 = 0;
   hcal_ndepth1 = 0;
   hcal_edepth11 = 0;
   hcal_ndepth11 = 0;
   hcal_edepth21 = 0;
   hcal_ndepth21 = 0;
   hcal_edepth31 = 0;
   hcal_ndepth31 = 0;
   hcal_edepth41 = 0;
   hcal_ndepth41 = 0;
   hcal_edepth51 = 0;
   hcal_ndepth51 = 0;
   hcal_edepth61 = 0;
   hcal_ndepth61 = 0;
   hcal_edepth71 = 0;
   hcal_ndepth71 = 0;
   hcal_edepth81 = 0;
   hcal_ndepth81 = 0;
   hcal_activeL1 = 0;
   hcal_edepth1 = 0;
   hcal_depthMatch1 = 0;
   hcal_edepthuncalibrated2 = 0;
   hcal_ndepth2 = 0;
   hcal_edepth12 = 0;
   hcal_ndepth12 = 0;
   hcal_edepth22 = 0;
   hcal_ndepth22 = 0;
   hcal_edepth32 = 0;
   hcal_ndepth32 = 0;
   hcal_edepth42 = 0;
   hcal_ndepth42 = 0;
   hcal_edepth52 = 0;
   hcal_ndepth52 = 0;
   hcal_edepth62 = 0;
   hcal_ndepth62 = 0;
   hcal_edepth72 = 0;
   hcal_ndepth72 = 0;
   hcal_edepth82 = 0;
   hcal_ndepth82 = 0;
   hcal_activeL2 = 0;
   hcal_edepth2 = 0;
   hcal_depthMatch2 = 0;
   hcal_edepthuncalibrated3 = 0;
   hcal_ndepth3 = 0;
   hcal_edepth13 = 0;
   hcal_ndepth13 = 0;
   hcal_edepth23 = 0;
   hcal_ndepth23 = 0;
   hcal_edepth33 = 0;
   hcal_ndepth33 = 0;
   hcal_edepth43 = 0;
   hcal_ndepth43 = 0;
   hcal_edepth53 = 0;
   hcal_ndepth53 = 0;
   hcal_edepth63 = 0;
   hcal_ndepth63 = 0;
   hcal_edepth73 = 0;
   hcal_ndepth73 = 0;
   hcal_edepth83 = 0;
   hcal_ndepth83 = 0;
   hcal_activeL3 = 0;
   hcal_edepth3 = 0;
   hcal_depthMatch3 = 0;
   hcal_edepthuncalibrated4 = 0;
   hcal_ndepth4 = 0;
   hcal_edepth14 = 0;
   hcal_ndepth14 = 0;
   hcal_edepth24 = 0;
   hcal_ndepth24 = 0;
   hcal_edepth34 = 0;
   hcal_ndepth34 = 0;
   hcal_edepth44 = 0;
   hcal_ndepth44 = 0;
   hcal_edepth54 = 0;
   hcal_ndepth54 = 0;
   hcal_edepth64 = 0;
   hcal_ndepth64 = 0;
   hcal_edepth74 = 0;
   hcal_ndepth74 = 0;
   hcal_edepth84 = 0;
   hcal_ndepth84 = 0;
   hcal_activeL4 = 0;
   hcal_edepth4 = 0;
   hcal_depthMatch4 = 0;
   hcal_edepthuncalibrated5 = 0;
   hcal_ndepth5 = 0;
   hcal_edepth15 = 0;
   hcal_ndepth15 = 0;
   hcal_edepth25 = 0;
   hcal_ndepth25 = 0;
   hcal_edepth35 = 0;
   hcal_ndepth35 = 0;
   hcal_edepth45 = 0;
   hcal_ndepth45 = 0;
   hcal_edepth55 = 0;
   hcal_ndepth55 = 0;
   hcal_edepth65 = 0;
   hcal_ndepth65 = 0;
   hcal_edepth75 = 0;
   hcal_ndepth75 = 0;
   hcal_edepth85 = 0;
   hcal_ndepth85 = 0;
   hcal_activeL5 = 0;
   hcal_edepth5 = 0;
   hcal_depthMatch5 = 0;
   hcal_edepthuncalibrated6 = 0;
   hcal_ndepth6 = 0;
   hcal_edepth16 = 0;
   hcal_ndepth16 = 0;
   hcal_edepth26 = 0;
   hcal_ndepth26 = 0;
   hcal_edepth36 = 0;
   hcal_ndepth36 = 0;
   hcal_edepth46 = 0;
   hcal_ndepth46 = 0;
   hcal_edepth56 = 0;
   hcal_ndepth56 = 0;
   hcal_edepth66 = 0;
   hcal_ndepth66 = 0;
   hcal_edepth76 = 0;
   hcal_ndepth76 = 0;
   hcal_edepth86 = 0;
   hcal_ndepth86 = 0;
   hcal_activeL6 = 0;
   hcal_edepth6 = 0;
   hcal_depthMatch6 = 0;
   hcal_edepthuncalibrated7 = 0;
   hcal_ndepth7 = 0;
   hcal_edepth17 = 0;
   hcal_ndepth17 = 0;
   hcal_edepth27 = 0;
   hcal_ndepth27 = 0;
   hcal_edepth37 = 0;
   hcal_ndepth37 = 0;
   hcal_edepth47 = 0;
   hcal_ndepth47 = 0;
   hcal_edepth57 = 0;
   hcal_ndepth57 = 0;
   hcal_edepth67 = 0;
   hcal_ndepth67 = 0;
   hcal_edepth77 = 0;
   hcal_ndepth77 = 0;
   hcal_edepth87 = 0;
   hcal_ndepth87 = 0;
   hcal_activeL7 = 0;
   hcal_edepth7 = 0;
   hcal_depthMatch7 = 0;
   activeLength = 0;
   eEcal = 0;
   Ecal_label = 0;
   hcal_ieta_o = 0;
   hcal_iphi_o = 0;
   hcal_edepth1_o = 0;
   hcal_edepth2_o = 0;
   hcal_edepth3_o = 0;
   hcal_edepth4_o = 0;
   hcal_edepth11_o = 0;
   hcal_edepth12_o = 0;
   hcal_edepth13_o = 0;
   hcal_edepth14_o = 0;
   hcal_edepth21_o = 0;
   hcal_edepth22_o = 0;
   hcal_edepth23_o = 0;
   hcal_edepth24_o = 0;
   hcal_edepth31_o = 0;
   hcal_edepth32_o = 0;
   hcal_edepth33_o = 0;
   hcal_edepth34_o = 0;
   hcal_edepth41_o = 0;
   hcal_edepth42_o = 0;
   hcal_edepth43_o = 0;
   hcal_edepth44_o = 0;
   hcal_edepth51_o = 0;
   hcal_edepth52_o = 0;
   hcal_edepth53_o = 0;
   hcal_edepth54_o = 0;
   hcal_edepth61_o = 0;
   hcal_edepth62_o = 0;
   hcal_edepth63_o = 0;
   hcal_edepth64_o = 0;
   hcal_edepth71_o = 0;
   hcal_edepth72_o = 0;
   hcal_edepth73_o = 0;
   hcal_edepth74_o = 0;
   hcal_edepth81_o = 0;
   hcal_edepth82_o = 0;
   hcal_edepth83_o = 0;
   hcal_edepth84_o = 0;
   hcal_edepth5_o = 0;
   hcal_edepth15_o = 0;
   hcal_edepth25_o = 0;
   hcal_edepth35_o = 0;
   hcal_edepth45_o = 0;
   hcal_edepth55_o = 0;
   hcal_edepth65_o = 0;
   hcal_edepth75_o = 0;
   hcal_edepth85_o = 0;
   hcal_edepth6_o = 0;
   hcal_edepth16_o = 0;
   hcal_edepth26_o = 0;
   hcal_edepth36_o = 0;
   hcal_edepth46_o = 0;
   hcal_edepth56_o = 0;
   hcal_edepth66_o = 0;
   hcal_edepth76_o = 0;
   hcal_edepth86_o = 0;
   hcal_edepth7_o = 0;
   hcal_edepth17_o = 0;
   hcal_edepth27_o = 0;
   hcal_edepth37_o = 0;
   hcal_edepth47_o = 0;
   hcal_edepth57_o = 0;
   hcal_edepth67_o = 0;
   hcal_edepth77_o = 0;
   hcal_edepth87_o = 0;
   hltresults = 0;
   all_triggers = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event_No", &Event_No, &b_Event_No);
   fChain->SetBranchAddress("Run_No", &Run_No, &b_Run_No);
   fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
   fChain->SetBranchAddress("BXNumber", &BXNumber, &b_BXNumber);
   fChain->SetBranchAddress("GoodVertex", &GoodVertex, &b_GoodVertex);
   fChain->SetBranchAddress("LumiScaler", &LumiScaler, &b_LumiScaler);
   fChain->SetBranchAddress("BunchLumi", &BunchLumi, &b_BunchLumi);
   fChain->SetBranchAddress("pileup", &pileup, &b_pileup);
   fChain->SetBranchAddress("pileupRMS", &pileupRMS, &b_pileupRMS);
   fChain->SetBranchAddress("PF_Muon", &PF_Muon, &b_PF_Muon);
   fChain->SetBranchAddress("Global_Muon", &Global_Muon, &b_Global_Muon);
   fChain->SetBranchAddress("Tracker_muon", &Tracker_muon, &b_Tracker_muon);
   fChain->SetBranchAddress("MuonIsTight", &MuonIsTight, &b_MuonIsTight);
   fChain->SetBranchAddress("MuonIsMedium", &MuonIsMedium, &b_MuonIsMedium);
   fChain->SetBranchAddress("charge_of_muon", &charge_of_muon, &b_charge_of_muon);
   fChain->SetBranchAddress("pt_of_muon", &pt_of_muon, &b_pt_of_muon);
   fChain->SetBranchAddress("eta_of_muon", &eta_of_muon, &b_eta_of_muon);
   fChain->SetBranchAddress("phi_of_muon", &phi_of_muon, &b_phi_of_muon);
   fChain->SetBranchAddress("energy_of_muon", &energy_of_muon, &b_energy_of_muon);
   fChain->SetBranchAddress("p_of_muon", &p_of_muon, &b_p_of_muon);
   fChain->SetBranchAddress("muon_trkKink", &muon_trkKink, &b_muon_trkKink);
   fChain->SetBranchAddress("muon_chi2LocalPosition", &muon_chi2LocalPosition, &b_muon_chi2LocalPosition);
   fChain->SetBranchAddress("muon_segComp", &muon_segComp, &b_muon_segComp);
   fChain->SetBranchAddress("TrackerLayer", &TrackerLayer, &b_TrackerLayer);
   fChain->SetBranchAddress("NumPixelLayers", &NumPixelLayers, &b_NumPixelLayers);
   fChain->SetBranchAddress("InnerTrackPixelHits", &InnerTrackPixelHits, &b_InnerTrackPixelHits);
   fChain->SetBranchAddress("innerTrack", &innerTrack, &b_innerTrack);
   fChain->SetBranchAddress("chiTracker", &chiTracker, &b_chiTracker);
   fChain->SetBranchAddress("DxyTracker", &DxyTracker, &b_DxyTracker);
   fChain->SetBranchAddress("DzTracker", &DzTracker, &b_DzTracker);
   fChain->SetBranchAddress("innerTrackpt", &innerTrackpt, &b_innerTrackpt);
   fChain->SetBranchAddress("innerTracketa", &innerTracketa, &b_innerTracketa);
   fChain->SetBranchAddress("innerTrackphi", &innerTrackphi, &b_innerTrackphi);
   fChain->SetBranchAddress("tight_validFraction", &tight_validFraction, &b_tight_validFraction);
   fChain->SetBranchAddress("OuterTrack", &OuterTrack, &b_OuterTrack);
   fChain->SetBranchAddress("OuterTrackChi", &OuterTrackChi, &b_OuterTrackChi);
   fChain->SetBranchAddress("OuterTrackPt", &OuterTrackPt, &b_OuterTrackPt);
   fChain->SetBranchAddress("OuterTrackEta", &OuterTrackEta, &b_OuterTrackEta);
   fChain->SetBranchAddress("OuterTrackPhi", &OuterTrackPhi, &b_OuterTrackPhi);
   fChain->SetBranchAddress("OuterTrackHits", &OuterTrackHits, &b_OuterTrackHits);
   fChain->SetBranchAddress("OuterTrackRHits", &OuterTrackRHits, &b_OuterTrackRHits);
   fChain->SetBranchAddress("GlobalTrack", &GlobalTrack, &b_GlobalTrack);
   fChain->SetBranchAddress("GlobalTrckPt", &GlobalTrckPt, &b_GlobalTrckPt);
   fChain->SetBranchAddress("GlobalTrckEta", &GlobalTrckEta, &b_GlobalTrckEta);
   fChain->SetBranchAddress("GlobalTrckPhi", &GlobalTrckPhi, &b_GlobalTrckPhi);
   fChain->SetBranchAddress("Global_Muon_Hits", &Global_Muon_Hits, &b_Global_Muon_Hits);
   fChain->SetBranchAddress("MatchedStations", &MatchedStations, &b_MatchedStations);
   fChain->SetBranchAddress("GlobTrack_Chi", &GlobTrack_Chi, &b_GlobTrack_Chi);
   fChain->SetBranchAddress("Tight_LongitudinalImpactparameter", &Tight_LongitudinalImpactparameter, &b_Tight_LongitudinalImpactparameter);
   fChain->SetBranchAddress("Tight_TransImpactparameter", &Tight_TransImpactparameter, &b_Tight_TransImpactparameter);
   fChain->SetBranchAddress("IsolationR04", &IsolationR04, &b_IsolationR04);
   fChain->SetBranchAddress("IsolationR03", &IsolationR03, &b_IsolationR03);
   fChain->SetBranchAddress("ecal_3into3", &ecal_3into3, &b_ecal_3into3);
   fChain->SetBranchAddress("hcal_3into3", &hcal_3into3, &b_hcal_3into3);
   fChain->SetBranchAddress("tracker_3into3", &tracker_3into3, &b_tracker_3into3);
   fChain->SetBranchAddress("matchedId", &matchedId, &b_matchedId);
   fChain->SetBranchAddress("hcal_cellHot", &hcal_cellHot, &b_hcal_cellHot);
   fChain->SetBranchAddress("hcal_cellHot_o", &hcal_cellHot_o, &b_hcal_cellHot_o);
   fChain->SetBranchAddress("pt_of_genMuon", &pt_of_genMuon, &b_pt_of_genMuon);
   fChain->SetBranchAddress("eta_of_genMuon", &eta_of_genMuon, &b_eta_of_genMuon);
   fChain->SetBranchAddress("phi_of_genMuon", &phi_of_genMuon, &b_phi_of_genMuon);
   fChain->SetBranchAddress("energy_of_genMuon", &energy_of_genMuon, &b_energy_of_genMuon);
   fChain->SetBranchAddress("ecal_1x1", &ecal_1x1, &b_ecal_1x1);
   fChain->SetBranchAddress("ecal_3x3", &ecal_3x3, &b_ecal_3x3);
   fChain->SetBranchAddress("ecal_5x5", &ecal_5x5, &b_ecal_5x5);
   fChain->SetBranchAddress("ecal_15x15", &ecal_15x15, &b_ecal_15x15);
   fChain->SetBranchAddress("ecal_25x25", &ecal_25x25, &b_ecal_25x25);
   fChain->SetBranchAddress("hcal_1x1", &hcal_1x1, &b_hcal_1x1);
   fChain->SetBranchAddress("ecal_detID", &ecal_detID, &b_ecal_detID);
   fChain->SetBranchAddress("hcal_detID", &hcal_detID, &b_hcal_detID);
   fChain->SetBranchAddress("ehcal_detID", &ehcal_detID, &b_ehcal_detID);
   fChain->SetBranchAddress("hcal_ieta", &hcal_ieta, &b_hcal_ieta);
   fChain->SetBranchAddress("hcal_iphi", &hcal_iphi, &b_hcal_iphi);
   fChain->SetBranchAddress("hcal_edepthuncalibrated1", &hcal_edepthuncalibrated1, &b_hcal_edepthuncalibrated1);
   fChain->SetBranchAddress("hcal_ndepth1", &hcal_ndepth1, &b_hcal_ndepth1);
   fChain->SetBranchAddress("hcal_edepth11", &hcal_edepth11, &b_hcal_edepth11);
   fChain->SetBranchAddress("hcal_ndepth11", &hcal_ndepth11, &b_hcal_ndepth11);
   fChain->SetBranchAddress("hcal_edepth21", &hcal_edepth21, &b_hcal_edepth21);
   fChain->SetBranchAddress("hcal_ndepth21", &hcal_ndepth21, &b_hcal_ndepth21);
   fChain->SetBranchAddress("hcal_edepth31", &hcal_edepth31, &b_hcal_edepth31);
   fChain->SetBranchAddress("hcal_ndepth31", &hcal_ndepth31, &b_hcal_ndepth31);
   fChain->SetBranchAddress("hcal_edepth41", &hcal_edepth41, &b_hcal_edepth41);
   fChain->SetBranchAddress("hcal_ndepth41", &hcal_ndepth41, &b_hcal_ndepth41);
   fChain->SetBranchAddress("hcal_edepth51", &hcal_edepth51, &b_hcal_edepth51);
   fChain->SetBranchAddress("hcal_ndepth51", &hcal_ndepth51, &b_hcal_ndepth51);
   fChain->SetBranchAddress("hcal_edepth61", &hcal_edepth61, &b_hcal_edepth61);
   fChain->SetBranchAddress("hcal_ndepth61", &hcal_ndepth61, &b_hcal_ndepth61);
   fChain->SetBranchAddress("hcal_edepth71", &hcal_edepth71, &b_hcal_edepth71);
   fChain->SetBranchAddress("hcal_ndepth71", &hcal_ndepth71, &b_hcal_ndepth71);
   fChain->SetBranchAddress("hcal_edepth81", &hcal_edepth81, &b_hcal_edepth81);
   fChain->SetBranchAddress("hcal_ndepth81", &hcal_ndepth81, &b_hcal_ndepth81);
   fChain->SetBranchAddress("hcal_activeL1", &hcal_activeL1, &b_hcal_activeL1);
   fChain->SetBranchAddress("hcal_edepth1", &hcal_edepth1, &b_hcal_edepth1);
   fChain->SetBranchAddress("hcal_depthMatch1", &hcal_depthMatch1, &b_hcal_depthMatch1);
   fChain->SetBranchAddress("hcal_edepthuncalibrated2", &hcal_edepthuncalibrated2, &b_hcal_edepthuncalibrated2);
   fChain->SetBranchAddress("hcal_ndepth2", &hcal_ndepth2, &b_hcal_ndepth2);
   fChain->SetBranchAddress("hcal_edepth12", &hcal_edepth12, &b_hcal_edepth12);
   fChain->SetBranchAddress("hcal_ndepth12", &hcal_ndepth12, &b_hcal_ndepth12);
   fChain->SetBranchAddress("hcal_edepth22", &hcal_edepth22, &b_hcal_edepth22);
   fChain->SetBranchAddress("hcal_ndepth22", &hcal_ndepth22, &b_hcal_ndepth22);
   fChain->SetBranchAddress("hcal_edepth32", &hcal_edepth32, &b_hcal_edepth32);
   fChain->SetBranchAddress("hcal_ndepth32", &hcal_ndepth32, &b_hcal_ndepth32);
   fChain->SetBranchAddress("hcal_edepth42", &hcal_edepth42, &b_hcal_edepth42);
   fChain->SetBranchAddress("hcal_ndepth42", &hcal_ndepth42, &b_hcal_ndepth42);
   fChain->SetBranchAddress("hcal_edepth52", &hcal_edepth52, &b_hcal_edepth52);
   fChain->SetBranchAddress("hcal_ndepth52", &hcal_ndepth52, &b_hcal_ndepth52);
   fChain->SetBranchAddress("hcal_edepth62", &hcal_edepth62, &b_hcal_edepth62);
   fChain->SetBranchAddress("hcal_ndepth62", &hcal_ndepth62, &b_hcal_ndepth62);
   fChain->SetBranchAddress("hcal_edepth72", &hcal_edepth72, &b_hcal_edepth72);
   fChain->SetBranchAddress("hcal_ndepth72", &hcal_ndepth72, &b_hcal_ndepth72);
   fChain->SetBranchAddress("hcal_edepth82", &hcal_edepth82, &b_hcal_edepth82);
   fChain->SetBranchAddress("hcal_ndepth82", &hcal_ndepth82, &b_hcal_ndepth82);
   fChain->SetBranchAddress("hcal_activeL2", &hcal_activeL2, &b_hcal_activeL2);
   fChain->SetBranchAddress("hcal_edepth2", &hcal_edepth2, &b_hcal_edepth2);
   fChain->SetBranchAddress("hcal_depthMatch2", &hcal_depthMatch2, &b_hcal_depthMatch2);
   fChain->SetBranchAddress("hcal_edepthuncalibrated3", &hcal_edepthuncalibrated3, &b_hcal_edepthuncalibrated3);
   fChain->SetBranchAddress("hcal_ndepth3", &hcal_ndepth3, &b_hcal_ndepth3);
   fChain->SetBranchAddress("hcal_edepth13", &hcal_edepth13, &b_hcal_edepth13);
   fChain->SetBranchAddress("hcal_ndepth13", &hcal_ndepth13, &b_hcal_ndepth13);
   fChain->SetBranchAddress("hcal_edepth23", &hcal_edepth23, &b_hcal_edepth23);
   fChain->SetBranchAddress("hcal_ndepth23", &hcal_ndepth23, &b_hcal_ndepth23);
   fChain->SetBranchAddress("hcal_edepth33", &hcal_edepth33, &b_hcal_edepth33);
   fChain->SetBranchAddress("hcal_ndepth33", &hcal_ndepth33, &b_hcal_ndepth33);
   fChain->SetBranchAddress("hcal_edepth43", &hcal_edepth43, &b_hcal_edepth43);
   fChain->SetBranchAddress("hcal_ndepth43", &hcal_ndepth43, &b_hcal_ndepth43);
   fChain->SetBranchAddress("hcal_edepth53", &hcal_edepth53, &b_hcal_edepth53);
   fChain->SetBranchAddress("hcal_ndepth53", &hcal_ndepth53, &b_hcal_ndepth53);
   fChain->SetBranchAddress("hcal_edepth63", &hcal_edepth63, &b_hcal_edepth63);
   fChain->SetBranchAddress("hcal_ndepth63", &hcal_ndepth63, &b_hcal_ndepth63);
   fChain->SetBranchAddress("hcal_edepth73", &hcal_edepth73, &b_hcal_edepth73);
   fChain->SetBranchAddress("hcal_ndepth73", &hcal_ndepth73, &b_hcal_ndepth73);
   fChain->SetBranchAddress("hcal_edepth83", &hcal_edepth83, &b_hcal_edepth83);
   fChain->SetBranchAddress("hcal_ndepth83", &hcal_ndepth83, &b_hcal_ndepth83);
   fChain->SetBranchAddress("hcal_activeL3", &hcal_activeL3, &b_hcal_activeL3);
   fChain->SetBranchAddress("hcal_edepth3", &hcal_edepth3, &b_hcal_edepth3);
   fChain->SetBranchAddress("hcal_depthMatch3", &hcal_depthMatch3, &b_hcal_depthMatch3);
   fChain->SetBranchAddress("hcal_edepthuncalibrated4", &hcal_edepthuncalibrated4, &b_hcal_edepthuncalibrated4);
   fChain->SetBranchAddress("hcal_ndepth4", &hcal_ndepth4, &b_hcal_ndepth4);
   fChain->SetBranchAddress("hcal_edepth14", &hcal_edepth14, &b_hcal_edepth14);
   fChain->SetBranchAddress("hcal_ndepth14", &hcal_ndepth14, &b_hcal_ndepth14);
   fChain->SetBranchAddress("hcal_edepth24", &hcal_edepth24, &b_hcal_edepth24);
   fChain->SetBranchAddress("hcal_ndepth24", &hcal_ndepth24, &b_hcal_ndepth24);
   fChain->SetBranchAddress("hcal_edepth34", &hcal_edepth34, &b_hcal_edepth34);
   fChain->SetBranchAddress("hcal_ndepth34", &hcal_ndepth34, &b_hcal_ndepth34);
   fChain->SetBranchAddress("hcal_edepth44", &hcal_edepth44, &b_hcal_edepth44);
   fChain->SetBranchAddress("hcal_ndepth44", &hcal_ndepth44, &b_hcal_ndepth44);
   fChain->SetBranchAddress("hcal_edepth54", &hcal_edepth54, &b_hcal_edepth54);
   fChain->SetBranchAddress("hcal_ndepth54", &hcal_ndepth54, &b_hcal_ndepth54);
   fChain->SetBranchAddress("hcal_edepth64", &hcal_edepth64, &b_hcal_edepth64);
   fChain->SetBranchAddress("hcal_ndepth64", &hcal_ndepth64, &b_hcal_ndepth64);
   fChain->SetBranchAddress("hcal_edepth74", &hcal_edepth74, &b_hcal_edepth74);
   fChain->SetBranchAddress("hcal_ndepth74", &hcal_ndepth74, &b_hcal_ndepth74);
   fChain->SetBranchAddress("hcal_edepth84", &hcal_edepth84, &b_hcal_edepth84);
   fChain->SetBranchAddress("hcal_ndepth84", &hcal_ndepth84, &b_hcal_ndepth84);
   fChain->SetBranchAddress("hcal_activeL4", &hcal_activeL4, &b_hcal_activeL4);
   fChain->SetBranchAddress("hcal_edepth4", &hcal_edepth4, &b_hcal_edepth4);
   fChain->SetBranchAddress("hcal_depthMatch4", &hcal_depthMatch4, &b_hcal_depthMatch4);
   fChain->SetBranchAddress("hcal_edepthuncalibrated5", &hcal_edepthuncalibrated5, &b_hcal_edepthuncalibrated5);
   fChain->SetBranchAddress("hcal_ndepth5", &hcal_ndepth5, &b_hcal_ndepth5);
   fChain->SetBranchAddress("hcal_edepth15", &hcal_edepth15, &b_hcal_edepth15);
   fChain->SetBranchAddress("hcal_ndepth15", &hcal_ndepth15, &b_hcal_ndepth15);
   fChain->SetBranchAddress("hcal_edepth25", &hcal_edepth25, &b_hcal_edepth25);
   fChain->SetBranchAddress("hcal_ndepth25", &hcal_ndepth25, &b_hcal_ndepth25);
   fChain->SetBranchAddress("hcal_edepth35", &hcal_edepth35, &b_hcal_edepth35);
   fChain->SetBranchAddress("hcal_ndepth35", &hcal_ndepth35, &b_hcal_ndepth35);
   fChain->SetBranchAddress("hcal_edepth45", &hcal_edepth45, &b_hcal_edepth45);
   fChain->SetBranchAddress("hcal_ndepth45", &hcal_ndepth45, &b_hcal_ndepth45);
   fChain->SetBranchAddress("hcal_edepth55", &hcal_edepth55, &b_hcal_edepth55);
   fChain->SetBranchAddress("hcal_ndepth55", &hcal_ndepth55, &b_hcal_ndepth55);
   fChain->SetBranchAddress("hcal_edepth65", &hcal_edepth65, &b_hcal_edepth65);
   fChain->SetBranchAddress("hcal_ndepth65", &hcal_ndepth65, &b_hcal_ndepth65);
   fChain->SetBranchAddress("hcal_edepth75", &hcal_edepth75, &b_hcal_edepth75);
   fChain->SetBranchAddress("hcal_ndepth75", &hcal_ndepth75, &b_hcal_ndepth75);
   fChain->SetBranchAddress("hcal_edepth85", &hcal_edepth85, &b_hcal_edepth85);
   fChain->SetBranchAddress("hcal_ndepth85", &hcal_ndepth85, &b_hcal_ndepth85);
   fChain->SetBranchAddress("hcal_activeL5", &hcal_activeL5, &b_hcal_activeL5);
   fChain->SetBranchAddress("hcal_edepth5", &hcal_edepth5, &b_hcal_edepth5);
   fChain->SetBranchAddress("hcal_depthMatch5", &hcal_depthMatch5, &b_hcal_depthMatch5);
   fChain->SetBranchAddress("hcal_edepthuncalibrated6", &hcal_edepthuncalibrated6, &b_hcal_edepthuncalibrated6);
   fChain->SetBranchAddress("hcal_ndepth6", &hcal_ndepth6, &b_hcal_ndepth6);
   fChain->SetBranchAddress("hcal_edepth16", &hcal_edepth16, &b_hcal_edepth16);
   fChain->SetBranchAddress("hcal_ndepth16", &hcal_ndepth16, &b_hcal_ndepth16);
   fChain->SetBranchAddress("hcal_edepth26", &hcal_edepth26, &b_hcal_edepth26);
   fChain->SetBranchAddress("hcal_ndepth26", &hcal_ndepth26, &b_hcal_ndepth26);
   fChain->SetBranchAddress("hcal_edepth36", &hcal_edepth36, &b_hcal_edepth36);
   fChain->SetBranchAddress("hcal_ndepth36", &hcal_ndepth36, &b_hcal_ndepth36);
   fChain->SetBranchAddress("hcal_edepth46", &hcal_edepth46, &b_hcal_edepth46);
   fChain->SetBranchAddress("hcal_ndepth46", &hcal_ndepth46, &b_hcal_ndepth46);
   fChain->SetBranchAddress("hcal_edepth56", &hcal_edepth56, &b_hcal_edepth56);
   fChain->SetBranchAddress("hcal_ndepth56", &hcal_ndepth56, &b_hcal_ndepth56);
   fChain->SetBranchAddress("hcal_edepth66", &hcal_edepth66, &b_hcal_edepth66);
   fChain->SetBranchAddress("hcal_ndepth66", &hcal_ndepth66, &b_hcal_ndepth66);
   fChain->SetBranchAddress("hcal_edepth76", &hcal_edepth76, &b_hcal_edepth76);
   fChain->SetBranchAddress("hcal_ndepth76", &hcal_ndepth76, &b_hcal_ndepth76);
   fChain->SetBranchAddress("hcal_edepth86", &hcal_edepth86, &b_hcal_edepth86);
   fChain->SetBranchAddress("hcal_ndepth86", &hcal_ndepth86, &b_hcal_ndepth86);
   fChain->SetBranchAddress("hcal_activeL6", &hcal_activeL6, &b_hcal_activeL6);
   fChain->SetBranchAddress("hcal_edepth6", &hcal_edepth6, &b_hcal_edepth6);
   fChain->SetBranchAddress("hcal_depthMatch6", &hcal_depthMatch6, &b_hcal_depthMatch6);
   fChain->SetBranchAddress("hcal_edepthuncalibrated7", &hcal_edepthuncalibrated7, &b_hcal_edepthuncalibrated7);
   fChain->SetBranchAddress("hcal_ndepth7", &hcal_ndepth7, &b_hcal_ndepth7);
   fChain->SetBranchAddress("hcal_edepth17", &hcal_edepth17, &b_hcal_edepth17);
   fChain->SetBranchAddress("hcal_ndepth17", &hcal_ndepth17, &b_hcal_ndepth17);
   fChain->SetBranchAddress("hcal_edepth27", &hcal_edepth27, &b_hcal_edepth27);
   fChain->SetBranchAddress("hcal_ndepth27", &hcal_ndepth27, &b_hcal_ndepth27);
   fChain->SetBranchAddress("hcal_edepth37", &hcal_edepth37, &b_hcal_edepth37);
   fChain->SetBranchAddress("hcal_ndepth37", &hcal_ndepth37, &b_hcal_ndepth37);
   fChain->SetBranchAddress("hcal_edepth47", &hcal_edepth47, &b_hcal_edepth47);
   fChain->SetBranchAddress("hcal_ndepth47", &hcal_ndepth47, &b_hcal_ndepth47);
   fChain->SetBranchAddress("hcal_edepth57", &hcal_edepth57, &b_hcal_edepth57);
   fChain->SetBranchAddress("hcal_ndepth57", &hcal_ndepth57, &b_hcal_ndepth57);
   fChain->SetBranchAddress("hcal_edepth67", &hcal_edepth67, &b_hcal_edepth67);
   fChain->SetBranchAddress("hcal_ndepth67", &hcal_ndepth67, &b_hcal_ndepth67);
   fChain->SetBranchAddress("hcal_edepth77", &hcal_edepth77, &b_hcal_edepth77);
   fChain->SetBranchAddress("hcal_ndepth77", &hcal_ndepth77, &b_hcal_ndepth77);
   fChain->SetBranchAddress("hcal_edepth87", &hcal_edepth87, &b_hcal_edepth87);
   fChain->SetBranchAddress("hcal_ndepth87", &hcal_ndepth87, &b_hcal_ndepth87);
   fChain->SetBranchAddress("hcal_activeL7", &hcal_activeL7, &b_hcal_activeL7);
   fChain->SetBranchAddress("hcal_edepth7", &hcal_edepth7, &b_hcal_edepth7);
   fChain->SetBranchAddress("hcal_depthMatch7", &hcal_depthMatch7, &b_hcal_depthMatch7);
   fChain->SetBranchAddress("activeLength", &activeLength, &b_activeLength);
   fChain->SetBranchAddress("eEcal", &eEcal, &b_eEcal);
   fChain->SetBranchAddress("Ecal_label", &Ecal_label, &b_Ecal_label);
   fChain->SetBranchAddress("hcal_ieta_o", &hcal_ieta_o, &b_hcal_ieta_o);
   fChain->SetBranchAddress("hcal_iphi_o", &hcal_iphi_o, &b_hcal_iphi_o);
   fChain->SetBranchAddress("hcal_edepth1_o", &hcal_edepth1_o, &b_hcal_edepth1_o);
   fChain->SetBranchAddress("hcal_edepth2_o", &hcal_edepth2_o, &b_hcal_edepth2_o);
   fChain->SetBranchAddress("hcal_edepth3_o", &hcal_edepth3_o, &b_hcal_edepth3_o);
   fChain->SetBranchAddress("hcal_edepth4_o", &hcal_edepth4_o, &b_hcal_edepth4_o);
   fChain->SetBranchAddress("hcal_edepth11_o", &hcal_edepth11_o, &b_hcal_edepth11_o);
   fChain->SetBranchAddress("hcal_edepth12_o", &hcal_edepth12_o, &b_hcal_edepth12_o);
   fChain->SetBranchAddress("hcal_edepth13_o", &hcal_edepth13_o, &b_hcal_edepth13_o);
   fChain->SetBranchAddress("hcal_edepth14_o", &hcal_edepth14_o, &b_hcal_edepth14_o);
   fChain->SetBranchAddress("hcal_edepth21_o", &hcal_edepth21_o, &b_hcal_edepth21_o);
   fChain->SetBranchAddress("hcal_edepth22_o", &hcal_edepth22_o, &b_hcal_edepth22_o);
   fChain->SetBranchAddress("hcal_edepth23_o", &hcal_edepth23_o, &b_hcal_edepth23_o);
   fChain->SetBranchAddress("hcal_edepth24_o", &hcal_edepth24_o, &b_hcal_edepth24_o);
   fChain->SetBranchAddress("hcal_edepth31_o", &hcal_edepth31_o, &b_hcal_edepth31_o);
   fChain->SetBranchAddress("hcal_edepth32_o", &hcal_edepth32_o, &b_hcal_edepth32_o);
   fChain->SetBranchAddress("hcal_edepth33_o", &hcal_edepth33_o, &b_hcal_edepth33_o);
   fChain->SetBranchAddress("hcal_edepth34_o", &hcal_edepth34_o, &b_hcal_edepth34_o);
   fChain->SetBranchAddress("hcal_edepth41_o", &hcal_edepth41_o, &b_hcal_edepth41_o);
   fChain->SetBranchAddress("hcal_edepth42_o", &hcal_edepth42_o, &b_hcal_edepth42_o);
   fChain->SetBranchAddress("hcal_edepth43_o", &hcal_edepth43_o, &b_hcal_edepth43_o);
   fChain->SetBranchAddress("hcal_edepth44_o", &hcal_edepth44_o, &b_hcal_edepth44_o);
   fChain->SetBranchAddress("hcal_edepth51_o", &hcal_edepth51_o, &b_hcal_edepth51_o);
   fChain->SetBranchAddress("hcal_edepth52_o", &hcal_edepth52_o, &b_hcal_edepth52_o);
   fChain->SetBranchAddress("hcal_edepth53_o", &hcal_edepth53_o, &b_hcal_edepth53_o);
   fChain->SetBranchAddress("hcal_edepth54_o", &hcal_edepth54_o, &b_hcal_edepth54_o);
   fChain->SetBranchAddress("hcal_edepth61_o", &hcal_edepth61_o, &b_hcal_edepth61_o);
   fChain->SetBranchAddress("hcal_edepth62_o", &hcal_edepth62_o, &b_hcal_edepth62_o);
   fChain->SetBranchAddress("hcal_edepth63_o", &hcal_edepth63_o, &b_hcal_edepth63_o);
   fChain->SetBranchAddress("hcal_edepth64_o", &hcal_edepth64_o, &b_hcal_edepth64_o);
   fChain->SetBranchAddress("hcal_edepth71_o", &hcal_edepth71_o, &b_hcal_edepth71_o);
   fChain->SetBranchAddress("hcal_edepth72_o", &hcal_edepth72_o, &b_hcal_edepth72_o);
   fChain->SetBranchAddress("hcal_edepth73_o", &hcal_edepth73_o, &b_hcal_edepth73_o);
   fChain->SetBranchAddress("hcal_edepth74_o", &hcal_edepth74_o, &b_hcal_edepth74_o);
   fChain->SetBranchAddress("hcal_edepth81_o", &hcal_edepth81_o, &b_hcal_edepth81_o);
   fChain->SetBranchAddress("hcal_edepth82_o", &hcal_edepth82_o, &b_hcal_edepth82_o);
   fChain->SetBranchAddress("hcal_edepth83_o", &hcal_edepth83_o, &b_hcal_edepth83_o);
   fChain->SetBranchAddress("hcal_edepth84_o", &hcal_edepth84_o, &b_hcal_edepth84_o);
   fChain->SetBranchAddress("hcal_edepth5_o", &hcal_edepth5_o, &b_hcal_edepth5_o);
   fChain->SetBranchAddress("hcal_edepth15_o", &hcal_edepth15_o, &b_hcal_edepth15_o);
   fChain->SetBranchAddress("hcal_edepth25_o", &hcal_edepth25_o, &b_hcal_edepth25_o);
   fChain->SetBranchAddress("hcal_edepth35_o", &hcal_edepth35_o, &b_hcal_edepth35_o);
   fChain->SetBranchAddress("hcal_edepth45_o", &hcal_edepth45_o, &b_hcal_edepth45_o);
   fChain->SetBranchAddress("hcal_edepth55_o", &hcal_edepth55_o, &b_hcal_edepth55_o);
   fChain->SetBranchAddress("hcal_edepth65_o", &hcal_edepth65_o, &b_hcal_edepth65_o);
   fChain->SetBranchAddress("hcal_edepth75_o", &hcal_edepth75_o, &b_hcal_edepth75_o);
   fChain->SetBranchAddress("hcal_edepth85_o", &hcal_edepth85_o, &b_hcal_edepth85_o);
   fChain->SetBranchAddress("hcal_edepth6_o", &hcal_edepth6_o, &b_hcal_edepth6_o);
   fChain->SetBranchAddress("hcal_edepth16_o", &hcal_edepth16_o, &b_hcal_edepth16_o);
   fChain->SetBranchAddress("hcal_edepth26_o", &hcal_edepth26_o, &b_hcal_edepth26_o);
   fChain->SetBranchAddress("hcal_edepth36_o", &hcal_edepth36_o, &b_hcal_edepth36_o);
   fChain->SetBranchAddress("hcal_edepth46_o", &hcal_edepth46_o, &b_hcal_edepth46_o);
   fChain->SetBranchAddress("hcal_edepth56_o", &hcal_edepth56_o, &b_hcal_edepth56_o);
   fChain->SetBranchAddress("hcal_edepth66_o", &hcal_edepth66_o, &b_hcal_edepth66_o);
   fChain->SetBranchAddress("hcal_edepth76_o", &hcal_edepth76_o, &b_hcal_edepth76_o);
   fChain->SetBranchAddress("hcal_edepth86_o", &hcal_edepth86_o, &b_hcal_edepth86_o);
   fChain->SetBranchAddress("hcal_edepth7_o", &hcal_edepth7_o, &b_hcal_edepth7_o);
   fChain->SetBranchAddress("hcal_edepth17_o", &hcal_edepth17_o, &b_hcal_edepth17_o);
   fChain->SetBranchAddress("hcal_edepth27_o", &hcal_edepth27_o, &b_hcal_edepth27_o);
   fChain->SetBranchAddress("hcal_edepth37_o", &hcal_edepth37_o, &b_hcal_edepth37_o);
   fChain->SetBranchAddress("hcal_edepth47_o", &hcal_edepth47_o, &b_hcal_edepth47_o);
   fChain->SetBranchAddress("hcal_edepth57_o", &hcal_edepth57_o, &b_hcal_edepth57_o);
   fChain->SetBranchAddress("hcal_edepth67_o", &hcal_edepth67_o, &b_hcal_edepth67_o);
   fChain->SetBranchAddress("hcal_edepth77_o", &hcal_edepth77_o, &b_hcal_edepth77_o);
   fChain->SetBranchAddress("hcal_edepth87_o", &hcal_edepth87_o, &b_hcal_edepth87_o);
   fChain->SetBranchAddress("hltresults", &hltresults, &b_hltresults);
   fChain->SetBranchAddress("all_triggers", &all_triggers, &b_all_triggers);
   Notify();
}

Bool_t NtuplerReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtuplerReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtuplerReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NtuplerReader_cxx
