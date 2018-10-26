#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include "TPRegexp.h"
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//////////////trigger info////////////////////////////////////

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"

#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"
#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"

#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h" 
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h" 
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "Calibration/IsolatedParticles/interface/CaloPropagateTrack.h"
#include "Calibration/IsolatedParticles/interface/eECALMatrix.h" 
#include "Calibration/IsolatedParticles/interface/eHCALMatrix.h" 

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//#define EDM_ML_DEBUG

class HcalHBHEMuonAnalyzer :  public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {

public:
  explicit HcalHBHEMuonAnalyzer(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;
  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void endRun(edm::Run const&, edm::EventSetup const&) override ;
  void   clearVectors();
  int    matchId(const HcalDetId&, const HcalDetId&);
  double activeLength(const DetId&);
  bool   isGoodVertex(const reco::Vertex& vtx);
  double respCorr(const DetId& id);
  double gainFactor(const edm::ESHandle<HcalDbService>&, const HcalDetId& id);
  int    depth16HE(int ieta, int iphi);
  bool   goodCell(const HcalDetId& hcid, const reco::Track* pTrack, 
		  const CaloGeometry* geo, const MagneticField* bField);

  // ----------member data ---------------------------
  HLTConfigProvider          hltConfig_;
  edm::Service<TFileService> fs;
  const edm::InputTag        hlTriggerResults_;
  const edm::InputTag        labelEBRecHit_, labelEERecHit_, labelHBHERecHit_, labelLumiScalers_;
  const std::string          labelVtx_, labelMuon_, labelGenPart_,fileInCorr_;
  const std::vector<std::string> triggers_;
  const int                  verbosity_, useRaw_;
  const bool                 unCorrect_, collapseDepth_, isItPlan1_;
  const bool                 ignoreHECorr_, isItPreRecHit_;
  const bool                 getCharge_, writeRespCorr_;
  bool                       mergedDepth_, useMyCorr_;
  int                        maxDepth_, kount_;

  const HcalDDDRecConstants *hdc_;
  const HcalTopology        *theHBHETopology_;
  HcalRespCorrs             *respCorrs_;

  edm::EDGetTokenT<edm::TriggerResults>                   tok_trigRes_;
  edm::EDGetTokenT<reco::VertexCollection>                tok_Vtx_;
  edm::EDGetTokenT<EcalRecHitCollection>                  tok_EB_;
  edm::EDGetTokenT<EcalRecHitCollection>                  tok_EE_;
  edm::EDGetTokenT<HBHERecHitCollection>                  tok_HBHE_;
  edm::EDGetTokenT<reco::MuonCollection>                  tok_Muon_;
  edm::EDGetTokenT<LumiScalersCollection>                 lumiScalersSrc_; 
  //edm::EDGetTokenT<reco::GenParticleCollection>           tok_GenPart_;
  //////////////////////////////////////////////////////
  static const int          depthMax_ = 7;
  TTree                    *tree_;
  unsigned int              runNumber_, eventNumber_ , lumiNumber_, bxNumber_;
  unsigned int              goodVertex_;
  float                     LumiScaler_, BunchLumi_, pileup_, pileupRMS_;
  std::vector<bool>         muon_is_good_, muon_global_, muon_tracker_;
  std::vector<bool>         muon_is_tight_, muon_is_medium_, isMuonRec_;
  std::vector<double>       cGlob_, ptGlob_, etaGlob_, phiGlob_, energyMuon_, pMuon_;
  std::vector<double>       isolationR04_;
  std::vector<double>       ecalEnergy_, hcalEnergy_, hoEnergy_;
  std::vector<bool>         matchedId_, hcalHot_, o_hcalHot_;
  std::vector<double>       genMuon_pt_,genMuon_eta_,genMuon_phi_,genMuon_energy_;
  std::vector<double>       eEcal_, Ecal_label_, ecal1x1Energy_, ecal3x3Energy_, ecal5x5Energy_, ecal15x15Energy_, ecal25x25Energy_, hcal1x1Energy_;
  std::vector<unsigned int> ecalDetId_, hcalDetId_, ehcalDetId_;
  std::vector<int>          hcal_ieta_, hcal_iphi_, o_hcal_ieta_, o_hcal_iphi_;
  std::vector<double>       hcalDepthEnergy_[depthMax_], hcalDepthEnergy1_[depthMax_], hcalDepthEnergy2_[depthMax_], hcalDepthEnergy3_[depthMax_], hcalDepthEnergy4_[depthMax_], hcalDepthEnergy5_[depthMax_], hcalDepthEnergy6_[depthMax_], hcalDepthEnergy7_[depthMax_], hcalDepthEnergy8_[depthMax_];
  std::vector<double>       hcalDepthNumHits_[depthMax_], hcalDepthNumHits1_[depthMax_], hcalDepthNumHits2_[depthMax_], hcalDepthNumHits3_[depthMax_], hcalDepthNumHits4_[depthMax_], hcalDepthNumHits5_[depthMax_], hcalDepthNumHits6_[depthMax_], hcalDepthNumHits7_[depthMax_], hcalDepthNumHits8_[depthMax_];
  std::vector<double>       hcalDepthActiveLength_[depthMax_];
  std::vector<double>       hcalDepthEnergyHot_[depthMax_];
  std::vector<double>       hcalDepthActiveLengthHot_[depthMax_];
  std::vector<double>       hcalDepthChargeHot_[depthMax_];
  std::vector<double>       hcalDepthChargeHotBG_[depthMax_];
  std::vector<double>       hcalDepthEnergyCorr_[depthMax_];
  std::vector<double>       hcalDepthEnergyHotCorr_[depthMax_];
  std::vector<bool>         hcalDepthMatch_[depthMax_];
  std::vector<bool>         hcalDepthMatchHot_[depthMax_];

  std::vector<double>       o_hcalDepth1Energy_, o_hcalDepth1Energy1_, o_hcalDepth1Energy2_, o_hcalDepth1Energy3_, o_hcalDepth1Energy4_, o_hcalDepth1Energy5_, o_hcalDepth1Energy6_, o_hcalDepth1Energy7_, o_hcalDepth1Energy8_;
  std::vector<double>       o_hcalDepth2Energy_, o_hcalDepth2Energy1_, o_hcalDepth2Energy2_, o_hcalDepth2Energy3_, o_hcalDepth2Energy4_, o_hcalDepth2Energy5_, o_hcalDepth2Energy6_, o_hcalDepth2Energy7_, o_hcalDepth2Energy8_;
  std::vector<double>       o_hcalDepth3Energy_, o_hcalDepth3Energy1_, o_hcalDepth3Energy2_, o_hcalDepth3Energy3_, o_hcalDepth3Energy4_, o_hcalDepth3Energy5_, o_hcalDepth3Energy6_, o_hcalDepth3Energy7_, o_hcalDepth3Energy8_;
  std::vector<double>       o_hcalDepth4Energy_, o_hcalDepth4Energy1_, o_hcalDepth4Energy2_, o_hcalDepth4Energy3_, o_hcalDepth4Energy4_, o_hcalDepth4Energy5_, o_hcalDepth4Energy6_, o_hcalDepth4Energy7_, o_hcalDepth4Energy8_;
  std::vector<double>       o_hcalDepth5Energy_, o_hcalDepth5Energy1_, o_hcalDepth5Energy2_, o_hcalDepth5Energy3_, o_hcalDepth5Energy4_, o_hcalDepth5Energy5_, o_hcalDepth5Energy6_, o_hcalDepth5Energy7_, o_hcalDepth5Energy8_;
  std::vector<double>       o_hcalDepth6Energy_, o_hcalDepth6Energy1_, o_hcalDepth6Energy2_, o_hcalDepth6Energy3_, o_hcalDepth6Energy4_, o_hcalDepth6Energy5_, o_hcalDepth6Energy6_, o_hcalDepth6Energy7_, o_hcalDepth6Energy8_;
  std::vector<double>       o_hcalDepth7Energy_, o_hcalDepth7Energy1_, o_hcalDepth7Energy2_, o_hcalDepth7Energy3_, o_hcalDepth7Energy4_, o_hcalDepth7Energy5_, o_hcalDepth7Energy6_, o_hcalDepth7Energy7_, o_hcalDepth7Energy8_;

  std::vector<double>       hcalActiveLength_; //hcalActiveLengthHot_;
  std::vector<std::string>  all_triggers_;
  std::vector<int>          hltresults_;

  std::vector<HcalDDDRecConstants::HcalActiveLength> actHB, actHE;
  std::map<DetId,double>    corrValue_;
  ////////////////////////////////////////////////////////////
  TH1F *h_cutflow = new TH1F("h_cutflow","Cutflow distribution",10,0,11);;
};
HcalHBHEMuonAnalyzer::HcalHBHEMuonAnalyzer(const edm::ParameterSet& iConfig) :
  hlTriggerResults_(iConfig.getParameter<edm::InputTag>("hlTriggerResults")),
  labelEBRecHit_(iConfig.getParameter<edm::InputTag>("labelEBRecHit")),
  labelEERecHit_(iConfig.getParameter<edm::InputTag>("labelEERecHit")),
  labelHBHERecHit_(iConfig.getParameter<edm::InputTag>("labelHBHERecHit")),
  labelLumiScalers_(iConfig.getParameter<edm::InputTag>("labelLumiScalers")),
  labelVtx_(iConfig.getParameter<std::string>("labelVertex")),
  labelMuon_(iConfig.getParameter<std::string>("labelMuon")),
  //labelGenPart_(iConfig.getParameter<std::string>("labelGenPart")),
  fileInCorr_(iConfig.getUntrackedParameter<std::string>("fileInCorr","")),
  triggers_(iConfig.getParameter<std::vector<std::string>>("triggers")),
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity",0)),
  useRaw_(iConfig.getParameter<int>("useRaw")),
  unCorrect_(iConfig.getParameter<bool>("unCorrect")),
  collapseDepth_(iConfig.getParameter<bool>("collapseDepth")),
  isItPlan1_(iConfig.getParameter<bool>("isItPlan1")),
  ignoreHECorr_(iConfig.getUntrackedParameter<bool>("ignoreHECorr",false)),
  isItPreRecHit_(iConfig.getUntrackedParameter<bool>("isItPreRecHit",false)),
  getCharge_(iConfig.getParameter<bool>("getCharge")),
  writeRespCorr_(iConfig.getUntrackedParameter<bool>("writeRespCorr",false)),
  hdc_(nullptr), theHBHETopology_(nullptr), respCorrs_(nullptr) {
  
  usesResource(TFileService::kSharedResource);
  //now do what ever initialization is needed
  kount_            = 0;
  maxDepth_         = iConfig.getUntrackedParameter<int>("maxDepth",4);
  if      (maxDepth_ > depthMax_) maxDepth_ = depthMax_;
  else if (maxDepth_ < 1)         maxDepth_ = 4;
  std::string modnam = iConfig.getUntrackedParameter<std::string>("moduleName","");
  std::string procnm = iConfig.getUntrackedParameter<std::string>("processName","");

  mergedDepth_  = (!isItPreRecHit_) || (collapseDepth_);
  tok_trigRes_  = consumes<edm::TriggerResults>(hlTriggerResults_);
  tok_EB_       = consumes<EcalRecHitCollection>(labelEBRecHit_);
  tok_EE_       = consumes<EcalRecHitCollection>(labelEERecHit_);
  tok_HBHE_     = consumes<HBHERecHitCollection>(labelHBHERecHit_);
  lumiScalersSrc_ = consumes<LumiScalersCollection>(labelLumiScalers_);
  if (modnam.empty()) {
    tok_Vtx_      = consumes<reco::VertexCollection>(labelVtx_);
    tok_Muon_     = consumes<reco::MuonCollection>(labelMuon_);
    //tok_GenPart_     = consumes<reco::GenParticleCollection>(labelGenPart_);
    edm::LogVerbatim("HBHEMuon")  << "Labels used: Trig " << hlTriggerResults_
				  << " Vtx " << labelVtx_ << " EB " 
				  << labelEBRecHit_ << " EE "
				  << labelEERecHit_ << " HBHE " 
				  << labelHBHERecHit_ << " MU " << labelMuon_; //<< "GenPart " << labelGenPart_;
  } else {
    tok_Vtx_      = consumes<reco::VertexCollection>(edm::InputTag(modnam,labelVtx_,procnm));
    tok_Muon_     = consumes<reco::MuonCollection>(edm::InputTag(modnam,labelMuon_,procnm));
    //tok_GenPart_     = consumes<reco::GenParticleCollection>(edm::InputTag(modnam,labelGenPart_,procnm));
    edm::LogVerbatim("HBHEMuon")   << "Labels used Trig " << hlTriggerResults_
				   << "\n  Vtx  " << edm::InputTag(modnam,labelVtx_,procnm)
				   << "\n  EB   " << labelEBRecHit_
				   << "\n  EE   " << labelEERecHit_
				   << "\n  HBHE " << labelHBHERecHit_
				   << "\n  MU   " << edm::InputTag(modnam,labelMuon_,procnm);
				   //<< "\n GenPart "<<edm::InputTag(modnam,labelGenPart_,procnm);
  }

  if (!fileInCorr_.empty()) {
    std::ifstream infile(fileInCorr_.c_str());
    if (infile.is_open()) {
      while (true) {
	unsigned int id;
	double       cfac;
	infile >> id >> cfac;
	if (!infile.good()) break;
	corrValue_[DetId(id)] = cfac;
      }
      infile.close();
    }
  }
  useMyCorr_ = (!corrValue_.empty());
  edm::LogVerbatim("HBHEMuon")   << "Flags used: UseRaw " << useRaw_ 
				 << " GetCharge " << getCharge_ << " UnCorrect "
				 << unCorrect_ << " IgnoreHECorr "
				 << ignoreHECorr_ << " CollapseDepth "
				 << collapseDepth_ << ":" << mergedDepth_
				 << " IsItPlan1 " << isItPlan1_
				 << " IsItPreRecHit " << isItPreRecHit_ 
				 << " UseMyCorr " << useMyCorr_;
}

//
// member functions
//

// ------------ method called for each event  ------------
void HcalHBHEMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  ++kount_;
  clearVectors();
  runNumber_   = iEvent.id().run();
  eventNumber_ = iEvent.id().event();
  lumiNumber_  = iEvent.id().luminosityBlock();
  bxNumber_    = iEvent.bunchCrossing();
#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("HBHEMuon") << "Run " << runNumber_ << " Event "
			       << eventNumber_ << " Lumi " << lumiNumber_ 
			       << " BX " << bxNumber_ << std::endl;
#endif  
  edm::Handle<edm::TriggerResults> _Triggers;
  iEvent.getByToken(tok_trigRes_, _Triggers); 
#ifdef EDM_ML_DEBUG
  if ((verbosity_/10000)%10>0) 
    edm::LogVerbatim("HBHEMuon") << "Size of all triggers "  
				 << all_triggers_.size() << std::endl;
#endif
  int Ntriggers = all_triggers_.size();
#ifdef EDM_ML_DEBUG
  if ((verbosity_/10000)%10>0) 
    edm::LogVerbatim("HBHEMuon") << "Size of HLT MENU: " << _Triggers->size()
				 << std::endl;
#endif
  if (_Triggers.isValid()) {
    const edm::TriggerNames &triggerNames_ = iEvent.triggerNames(*_Triggers);
    std::vector<int> index;
    for (int i=0; i<Ntriggers; i++) {
      index.push_back(triggerNames_.triggerIndex(all_triggers_[i]));
      int triggerSize = int( _Triggers->size());
#ifdef EDM_ML_DEBUG
      if ((verbosity_/10000)%10>0) 
	edm::LogVerbatim("HBHEMuon") << "outside loop " << index[i]
				     << "\ntriggerSize " << triggerSize
				     << std::endl;
#endif
      if (index[i] < triggerSize) {
	hltresults_.push_back(_Triggers->accept(index[i]));
#ifdef EDM_ML_DEBUG
	if ((verbosity_/10000)%10>0) 
	  edm::LogVerbatim("HBHEMuon") << "Trigger_info " << triggerSize
				       << " triggerSize " << index[i]
				       << " trigger_index " << hltresults_.at(i)
				       << " hltresult" << std::endl;
#endif
      } else {
	if ((verbosity_/10000)%10>0) 
	  edm::LogVerbatim("HBHEMuon") << "Requested HLT path \"" 
				       << "\" does not exist\n";
      }
    }
  }

  // get handles to calogeometry and calotopology
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* geo = pG.product();
  
  edm::ESHandle<MagneticField> bFieldH;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldH);
  const MagneticField* bField = bFieldH.product();
  
  edm::ESHandle<EcalChannelStatus> ecalChStatus;
  iSetup.get<EcalChannelStatusRcd>().get(ecalChStatus);
  const EcalChannelStatus* theEcalChStatus = ecalChStatus.product();

  edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);

  edm::ESHandle<CaloTopology> theCaloTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopology);
  const CaloTopology *caloTopology = theCaloTopology.product();

  edm::ESHandle<HcalDbService> conditions;
  iSetup.get<HcalDbRecord>().get(conditions);

  edm::ESHandle<CaloTowerConstituentsMap> ct;
  iSetup.get<CaloGeometryRecord>().get(ct);
  const CaloTowerConstituentsMap* ctmap = ct.product();

  // Relevant blocks from iEvent
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByToken(tok_Vtx_, vtx);

  edm::Handle<LumiScalersCollection> hscalers;
  iEvent.getByToken(lumiScalersSrc_, hscalers);
  const LumiScalersCollection& lumiScalers = *hscalers;

  edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
  iEvent.getByToken(tok_EB_, barrelRecHitsHandle);
  edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
  iEvent.getByToken(tok_EE_, endcapRecHitsHandle);

  edm::Handle<HBHERecHitCollection> hbhe;
  iEvent.getByToken(tok_HBHE_, hbhe);

  edm::Handle<reco::MuonCollection> _Muon;
  iEvent.getByToken(tok_Muon_, _Muon);

  //edm::Handle<reco::GenParticleCollection> _GenPart;
  //iEvent.getByToken(tok_GenPart_, _GenPart);
  // require a good vertex
  math::XYZPoint pvx;
  goodVertex_ = 0;
  LumiScaler_ = 0;
  BunchLumi_ = -5;
  pileup_ = 0;
  pileupRMS_ = 0;

  if (!vtx.isValid()) {
#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("HBHEMuon") << "No Good Vertex found == Reject\n";
#endif
    return;
  }
  reco::VertexCollection::const_iterator firstGoodVertex = vtx->end();
  for (reco::VertexCollection::const_iterator it = vtx->begin(); 
       it != vtx->end(); it++) {
    if (isGoodVertex(*it)) {
      if (firstGoodVertex == vtx->end()) firstGoodVertex = it;
      ++goodVertex_;
    }
  }
  if (firstGoodVertex != vtx->end()) pvx = firstGoodVertex->position();

  if(!lumiScalers.empty()) { 
    LumiScaler_ = lumiScalers.front().instantLumi();
    BunchLumi_ = lumiScalers.front().bunchLumi();
    pileup_ = lumiScalers.front().pileup();
    pileupRMS_ = lumiScalers.front().pileupRMS();
  }

  if (_Muon.isValid() && barrelRecHitsHandle.isValid() && 
      endcapRecHitsHandle.isValid() && hbhe.isValid()) { 
    for (reco::MuonCollection::const_iterator RecMuon = _Muon->begin(); RecMuon!= _Muon->end(); ++RecMuon)  {

      if (RecMuon->innerTrack().isNonnull()) {
        const reco::Track* pTrack = (RecMuon->innerTrack()).get();
        spr::propagatedTrackID trackID = spr::propagateCALO(pTrack, geo, bField, (((verbosity_/100)%10>0)));

        HcalDetId check;
        std::pair<bool,HcalDetId> info = spr::propagateHCALBack(pTrack,  geo, bField, (((verbosity_/100)%10>0)));
        if (info.first) {
          check = info.second;
        }

        DetId closestCell(trackID.detIdHCAL);
        HcalDetId hcidt(closestCell.rawId());
        bool         tmpmatch(false);
        int ieta(-1000), iphi(-1000), o_ieta(-1000), o_iphi(-1000);
        if ((hcidt.ieta() == check.ieta()) && (hcidt.iphi() == check.iphi())) tmpmatch = true;
#ifdef EDM_ML_DEBUG
        edm::LogVerbatim("HBHEMuon") << "Front " << hcidt << " Back "
                                       << info.first << ":" << check
                                       << " Match " << tmpmatch;
#endif

        HcalSubdetector subdet = hcidt.subdet();
        ieta   = hcidt.ieta();
        iphi   = hcidt.iphi();
        bool hborhe = (std::abs(ieta) == 16);

        bool muon_is_tight = muon::isTightMuon(*RecMuon,*firstGoodVertex);
        double isoR04 = ((RecMuon->pfIsolationR04().sumChargedHadronPt + std::max(0.,RecMuon->pfIsolationR04().sumNeutralHadronEt + RecMuon->pfIsolationR04().sumPhotonEt - (0.5 *RecMuon->pfIsolationR04().sumPUPt))) / RecMuon->pt()) ;
        if ((RecMuon->pt()<20.0) || (!trackID.okHCAL) || !tmpmatch || fabs(ieta)<21 || (!muon_is_tight) || isoR04>0.15) continue;

	h_cutflow->Fill(0); //initial muon selections
        TLorentzVector lvMuon;
        lvMuon.SetPtEtaPhiM(RecMuon->pt(), RecMuon->eta(), RecMuon->phi(), 0.10566);
        bool isMuonRec=false;
        bool Zmm=false;
        for (reco::MuonCollection::const_iterator iMuon = _Muon->begin(); iMuon!= _Muon->end(); ++iMuon)  {
         if(iMuon!=RecMuon && (iMuon->charge()!=RecMuon->charge())){
           bool muon1 = iMuon->pt()>20.0 && fabs(iMuon->eta())<2.5;
           bool muon2 = muon::isMediumMuon(*iMuon);
           bool muon3 = (iMuon->pfIsolationR04().sumChargedHadronPt + std::max(0.,iMuon->pfIsolationR04().sumNeutralHadronEt + iMuon->pfIsolationR04().sumPhotonEt - (0.5 *iMuon->pfIsolationR04().sumPUPt))) / iMuon->pt()<0.15;
           if(muon1 && muon2 && muon3){
            const reco::Track* ipTrack = (iMuon->innerTrack()).get();
            spr::propagatedTrackID itrackID = spr::propagateCALO(ipTrack, geo, bField, (((verbosity_/100)%10>0)));
            DetId iclosestCell(itrackID.detIdHCAL);
            HcalDetId ihcidt(iclosestCell.rawId());
            if((hcidt.ieta()==26 && ieta>0) || (hcidt.ieta()==-26 && ieta<0)){ isMuonRec=true;}

            if(fabs(hcidt.ieta())<26){
              TLorentzVector Muon_tmp;
              Muon_tmp.SetPtEtaPhiM(iMuon->pt(), iMuon->eta(), iMuon->phi(), 0.10566);
              double dimuon_M = (Muon_tmp + lvMuon).M();
                if(fabs(dimuon_M-91.188)< 2. ){
                  Zmm=true;
                  break;
                }
              }
            } 
          }
        }
        if(Zmm) continue;
	h_cutflow->Fill(1); // Zmm cut
        muon_is_good_.push_back(RecMuon->isPFMuon());
        muon_global_.push_back(RecMuon->isGlobalMuon());
        muon_tracker_.push_back(RecMuon->isTrackerMuon());
        cGlob_.push_back(RecMuon->charge());
        ptGlob_.push_back((RecMuon)->pt());
        etaGlob_.push_back(RecMuon->eta());
        phiGlob_.push_back(RecMuon->phi());
        energyMuon_.push_back(RecMuon->energy());
        pMuon_.push_back(RecMuon->p());
        muon_is_tight_.push_back(muon_is_tight);
        isolationR04_.push_back(isoR04);
        ecalEnergy_.push_back(RecMuon->calEnergy().emS9);		 
        hcalEnergy_.push_back(RecMuon->calEnergy().hadS9);
        hoEnergy_.push_back(RecMuon->calEnergy().hoS9);
        /*
	for (reco::GenParticleCollection::const_iterator iGenPart = _GenPart->begin(); iGenPart!= _GenPart->end(); ++iGenPart)  {
	  if(fabs(iGenPart->pdgId())==13){
	    genMuon_pt_.push_back(iGenPart->pt());
	    genMuon_eta_.push_back(iGenPart->eta());
	    genMuon_phi_.push_back(iGenPart->phi());
	    genMuon_energy_.push_back(iGenPart->energy());
	  }
	}
        */
        double eEcal1x1(0), eEcal3x3(0),eEcal5x5(0), eEcal15x15(0), eEcal25x25(0), eHcal(0), activeLengthTot(0); //activeLengthHotTot(0);
        double eHcalDepth[depthMax_], eHcalDepthHot[depthMax_];
        double eHcalDepth1[depthMax_], eHcalDepth2[depthMax_], eHcalDepth3[depthMax_], eHcalDepth4[depthMax_], eHcalDepth5[depthMax_], eHcalDepth6[depthMax_], eHcalDepth7[depthMax_], eHcalDepth8[depthMax_];
        double o_eHcalDepth[depthMax_], o_eHcalDepth1[depthMax_], o_eHcalDepth2[depthMax_], o_eHcalDepth3[depthMax_], o_eHcalDepth4[depthMax_], o_eHcalDepth5[depthMax_], o_eHcalDepth6[depthMax_], o_eHcalDepth7[depthMax_], o_eHcalDepth8[depthMax_];
        double eHcalDepthC[depthMax_], eHcalDepthHotC[depthMax_];
        //double cHcalDepthHot[depthMax_], cHcalDepthHotBG[depthMax_];
        double activeL[depthMax_], activeHotL[depthMax_];
        bool   matchDepth[depthMax_], matchDepthHot[depthMax_];
        HcalDetId eHcalDetId[depthMax_];
        unsigned int isHot(0), o_isHot(0);
        for (int i=0; i<depthMax_; ++i) {
       	   eHcalDepth[i]    = eHcalDepthHot[i]  = 0;
           eHcalDepth1[i] = eHcalDepth2[i] = eHcalDepth3[i] = eHcalDepth4[i] = eHcalDepth5[i] = eHcalDepth6[i] = eHcalDepth7[i] = eHcalDepth8[i] = 0;
           o_eHcalDepth[i] = o_eHcalDepth1[i] = o_eHcalDepth2[i] = o_eHcalDepth3[i] = o_eHcalDepth4[i] = o_eHcalDepth5[i] = o_eHcalDepth6[i] = o_eHcalDepth7[i] = o_eHcalDepth8[i] = 0.0;

           eHcalDepthC[i]   = eHcalDepthHotC[i] = 0;
	   //cHcalDepthHot[i] = cHcalDepthHotBG[i]= 0;
	   activeL[i]       = activeHotL[i]     = 0;
	   matchDepth[i]    = matchDepthHot[i]  = true;
        }

	
	bool okE = trackID.okECAL;
	if (okE) {
	  const DetId isoCell(trackID.detIdECAL);
	  std::pair<double,bool> e3x3 = spr::eECALmatrix(isoCell,barrelRecHitsHandle,endcapRecHitsHandle,*theEcalChStatus,geo,caloTopology,sevlv.product(),1,1,-100.0,-100.0,-500.0,500.0,false);
	  eEcal3x3 = e3x3.first;
	  okE   = e3x3.second;

          std::pair<double,bool> e1x1 = spr::eECALmatrix(isoCell,barrelRecHitsHandle,endcapRecHitsHandle,*theEcalChStatus,geo,caloTopology,sevlv.product(),0,0,-100.0,-100.0,-500.0,500.0,false);
          eEcal1x1 = e1x1.first;
          std::pair<double,bool> e5x5 = spr::eECALmatrix(isoCell,barrelRecHitsHandle,endcapRecHitsHandle,*theEcalChStatus,geo,caloTopology,sevlv.product(),2,2,-100.0,-100.0,-500.0,500.0,false);
          eEcal5x5 = e5x5.first;
          std::pair<double,bool> e15x15 = spr::eECALmatrix(isoCell,barrelRecHitsHandle,endcapRecHitsHandle,*theEcalChStatus,geo,caloTopology,sevlv.product(),7,7,-100.0,-100.0,-500.0,500.0,false);
          eEcal15x15 = e15x15.first;
          std::pair<double,bool> e25x25 = spr::eECALmatrix(isoCell,barrelRecHitsHandle,endcapRecHitsHandle,*theEcalChStatus,geo,caloTopology,sevlv.product(),12,12,-100.0,-100.0,-500.0,500.0,false);
          eEcal25x25 = e25x25.first;

	}
#ifdef EDM_ML_DEBUG
	edm::LogVerbatim("HBHEMuon") << "Propagate Track to ECAL: " << okE 
				     << ":" << trackID.okECAL << " E " <<eEcal3x3;
#endif

        //check the opposite side of the detector
        o_ieta = 26;
        o_iphi = iphi - 36;
        if(iphi<=36) o_iphi = iphi + 36;
        if(ieta<0) o_ieta = -26;
          
        HcalSubdetector o_subdet = HcalEndcap;
        HcalDetId o_hcid0(o_subdet,o_ieta,o_iphi,1);
        HcalDetId o_hcid1(o_subdet,o_ieta-1,o_iphi-2,1);
        HcalDetId o_hcid2(o_subdet,o_ieta-1,o_iphi,1);
        HcalDetId o_hcid3(o_subdet,o_ieta-1,o_iphi+2,1);

        int o_ieta_temp = o_ieta+1;

        HcalDetId o_hcid4(o_subdet,o_ieta_temp,o_iphi-2,1);
        HcalDetId o_hcid5(o_subdet,o_ieta_temp,o_iphi,1);
        HcalDetId o_hcid6(o_subdet,o_ieta_temp,o_iphi+2,1);

        HcalDetId o_hcid7(o_subdet,o_ieta,o_iphi-2,1);
        HcalDetId o_hcid8(o_subdet,o_ieta,o_iphi+2,1);

        const DetId o_HcalCell(o_hcid0);
        HcalDetId     o_hotCell;
        spr::eHCALmatrix(geo, theHBHETopology_, o_HcalCell, hbhe, 1,1, o_hotCell, false, useRaw_, false);
        o_isHot = matchId(o_HcalCell,o_hotCell);

        std::vector<std::pair<double,int> > o_ehdepth0, o_ehdepth1, o_ehdepth2, o_ehdepth3, o_ehdepth4, o_ehdepth5, o_ehdepth6, o_ehdepth7, o_ehdepth8;
        spr::energyHCALCell((HcalDetId)o_hcid0, hbhe, o_ehdepth0, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)o_hcid1, hbhe, o_ehdepth1, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)o_hcid2, hbhe, o_ehdepth2, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)o_hcid3, hbhe, o_ehdepth3, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)o_hcid4, hbhe, o_ehdepth4, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)o_hcid5, hbhe, o_ehdepth5, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)o_hcid6, hbhe, o_ehdepth6, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)o_hcid7, hbhe, o_ehdepth7, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)o_hcid8, hbhe, o_ehdepth8, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(o_ieta,o_iphi), (((verbosity_/1000)%10)>0));
        std::pair<double,bool> ecal = spr::eECALmatrix(o_HcalCell,barrelRecHitsHandle,endcapRecHitsHandle,geo,ctmap,sevlv.product(),-100.0,-100.0,-500.0,500.0,false);
        eEcal_.push_back(ecal.first);
        Ecal_label_.push_back(ecal.second);
        for (unsigned int i=0; i<o_ehdepth1.size(); ++i) {
               o_eHcalDepth1[o_ehdepth1[i].second-1] = o_ehdepth1[i].first;
        }

        for (unsigned int i=0; i<o_ehdepth2.size(); ++i) {
               o_eHcalDepth2[o_ehdepth2[i].second-1] = o_ehdepth2[i].first;
        }

        for (unsigned int i=0; i<o_ehdepth3.size(); ++i) {
               o_eHcalDepth3[o_ehdepth3[i].second-1] = o_ehdepth3[i].first;
        }

        for (unsigned int i=0; i<o_ehdepth4.size(); ++i) {
               o_eHcalDepth4[o_ehdepth4[i].second-1] = o_ehdepth4[i].first;
        }

        for (unsigned int i=0; i<o_ehdepth5.size(); ++i) {
               o_eHcalDepth5[o_ehdepth5[i].second-1] = o_ehdepth5[i].first;
        }

        for (unsigned int i=0; i<o_ehdepth6.size(); ++i) {
               o_eHcalDepth6[o_ehdepth6[i].second-1] = o_ehdepth6[i].first;
        }

        for (unsigned int i=0; i<o_ehdepth7.size(); ++i) {
               o_eHcalDepth7[o_ehdepth7[i].second-1] = o_ehdepth7[i].first;
        }

        for (unsigned int i=0; i<o_ehdepth8.size(); ++i) {
               o_eHcalDepth8[o_ehdepth8[i].second-1] = o_ehdepth8[i].first;
        }

        for (unsigned int i=0; i<o_ehdepth0.size(); ++i) {
               o_eHcalDepth[o_ehdepth0[i].second-1] = o_ehdepth0[i].first;
               if(o_ehdepth0[i].second==2) edm::LogVerbatim("HBHEMuon") <<"depth "<<o_ehdepth0[i].second-1<<"energy: "<<o_ehdepth0[i].first;
        }


        eHcal = spr::eHCALmatrix(theHBHETopology_, closestCell, hbhe,0,0, false, true, -100.0, -100.0, -100.0, -100.0, -500.,500.,useRaw_);
	std::vector<std::pair<double,int> > ehdepth;
	spr::energyHCALCell((HcalDetId)closestCell, hbhe, ehdepth, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));
        int iphi_temp1 = iphi-1;
        int iphi_temp2 = iphi+1;
        if(ieta>20){ 
             iphi_temp1 = iphi-2; 
             iphi_temp2 = iphi+2;
        }
          
        HcalDetId hcid1(subdet,ieta-1,iphi_temp1,1);
        HcalDetId hcid2(subdet,ieta-1,iphi,1);
        HcalDetId hcid3(subdet,ieta-1,iphi_temp2,1);
        int ieta_temp = ieta+1;
        HcalDetId hcid4(subdet,ieta_temp,iphi_temp1,1);
        HcalDetId hcid5(subdet,ieta_temp,iphi,1);
        HcalDetId hcid6(subdet,ieta_temp,iphi_temp2,1);
        HcalDetId hcid7(subdet,ieta,iphi_temp1,1);
        HcalDetId hcid8(subdet,ieta,iphi_temp2,1);

        std::vector<std::pair<double,int> > ehdepth1, ehdepth2, ehdepth3, ehdepth4, ehdepth5, ehdepth6, ehdepth7, ehdepth8;
          //std::vector<std::pair<int,int> > numhits1, numhits2, numhits3, numhits4, numhits5, numhits6, numhits7, numhits8;
        spr::energyHCALCell((HcalDetId)hcid1, hbhe, ehdepth1, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)hcid2, hbhe, ehdepth2, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)hcid3, hbhe, ehdepth3, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)hcid4, hbhe, ehdepth4, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)hcid5, hbhe, ehdepth5, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)hcid6, hbhe, ehdepth6, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)hcid7, hbhe, ehdepth7, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));
        spr::energyHCALCell((HcalDetId)hcid8, hbhe, ehdepth8, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depth16HE(ieta,iphi), (((verbosity_/1000)%10)>0));

        for (unsigned int i=0; i<ehdepth1.size(); ++i) {
            eHcalDepth1[ehdepth1[i].second-1] = ehdepth1[i].first;
            //eHcalNumHits1[ehdepth1[i].second-1] = numhits1[i].first;
        }

        for (unsigned int i=0; i<ehdepth2.size(); ++i) {
            eHcalDepth2[ehdepth2[i].second-1] = ehdepth2[i].first;
            //eHcalNumHits2[ehdepth2[i].second-1] = numhits2[i].first;
        }

        for (unsigned int i=0; i<ehdepth3.size(); ++i) {
            eHcalDepth3[ehdepth3[i].second-1] = ehdepth3[i].first;
            //eHcalNumHits3[ehdepth3[i].second-1] = numhits3[i].first;
        }

        for (unsigned int i=0; i<ehdepth4.size(); ++i) {
            eHcalDepth4[ehdepth4[i].second-1] = ehdepth4[i].first;
            //eHcalNumHits4[ehdepth4[i].second-1] = numhits4[i].first;
        }

        for (unsigned int i=0; i<ehdepth5.size(); ++i) {
            eHcalDepth5[ehdepth5[i].second-1] = ehdepth5[i].first;
            //eHcalNumHits5[ehdepth5[i].second-1] = numhits5[i].first;
        }

        for (unsigned int i=0; i<ehdepth6.size(); ++i) {
            eHcalDepth6[ehdepth6[i].second-1] = ehdepth6[i].first;
            //eHcalNumHits6[ehdepth6[i].second-1] = numhits6[i].first;
        }

        for (unsigned int i=0; i<ehdepth7.size(); ++i) {
            eHcalDepth7[ehdepth7[i].second-1] = ehdepth7[i].first;
            //eHcalNumHits7[ehdepth7[i].second-1] = numhits7[i].first;
        }

        for (unsigned int i=0; i<ehdepth8.size(); ++i) {
            eHcalDepth8[ehdepth8[i].second-1] = ehdepth8[i].first;
            //eHcalNumHits8[ehdepth8[i].second-1] = numhits8[i].first;
        }

	for (int i=0; i<depthMax_; ++i) eHcalDetId[i] = HcalDetId();
	for (unsigned int i=0; i<ehdepth.size(); ++i) {
	    HcalSubdetector subdet0 = (hborhe) ? ((ehdepth[i].second >= depth16HE(ieta,iphi)) ? HcalEndcap : HcalBarrel) : subdet;
	    HcalDetId hcid0(subdet0,ieta,iphi,ehdepth[i].second);
	    double actL = activeLength(DetId(hcid0));
	    double ene  = ehdepth[i].first;
            //int nhits = numhits[i].first;
	    bool   tmpC(false);
	    if (ene > 0.0) {
	      if (!(theHBHETopology_->validHcal(hcid0))) {
		edm::LogWarning("HBHEMuon") << "(1) Invalid ID " << hcid0 
					    << " with E = " << ene;
		edm::LogWarning("HBHEMuon") << HcalDetId(closestCell) 
					    << " with " << ehdepth.size() 
					    << " depths:";
		for (const auto& ehd : ehdepth) 
		  edm::LogWarning("HBHEMuon") << " " << ehd.second << ":" 
					      << ehd.first;
	      } else {
		tmpC = goodCell(hcid0, pTrack, geo, bField);
		double enec(ene);
		if (unCorrect_) {
		  double corr = (ignoreHECorr_ && (subdet0==HcalEndcap)) ? 1.0 : respCorr(DetId(hcid0));
		  if (corr != 0) ene /= corr;
#ifdef EDM_ML_DEBUG
		  HcalDetId id = (isItPlan1_ && isItPreRecHit_) ? hdc_->mergedDepthDetId(hcid0) : hcid0;
		  edm::LogVerbatim("HBHEMuon") << hcid0 << ":" << id << " Corr "
					       << corr;
#endif
		}
		int depth = ehdepth[i].second - 1;
		if (collapseDepth_) {
		  HcalDetId id = hdc_->mergedDepthDetId(hcid0);
		  depth        = id.depth() - 1;
		}
		eHcalDepth[depth] += ene;
                //eHcalNumHits[depth] += nhits;
		eHcalDepthC[depth]+= enec;
		activeL[depth]    += actL;
		activeLengthTot   += actL;
		matchDepth[depth]  = (matchDepth[depth] && tmpC);
#ifdef EDM_ML_DEBUG
		if ((verbosity_%10) > 0)
		  edm::LogVerbatim("HBHEMuon") << hcid0 << " E " << ene << ":"
					       << enec << " L " << actL 
					       << " Match " << tmpC;
#endif
	      }
	    }
	}
#ifdef EDM_ML_DEBUG
	if ((verbosity_%10) > 0) {
	    edm::LogVerbatim("HBHEMuon") << hcidt << " Match " << tmpmatch 
					 << " Depths " << ehdepth.size();
	    for (unsigned int k=0; k<ehdepth.size(); ++k)
	      edm::LogVerbatim("HBHEMuon") << " [" << k << ":"
					   << ehdepth[k].second << "] " 
					   << matchDepth[k];
	}
#endif
	HcalDetId           hotCell;
	spr::eHCALmatrix(geo, theHBHETopology_, closestCell, hbhe, 1,1, hotCell, false, useRaw_, false);
	isHot = matchId(closestCell,hotCell);

        isMuonRec_.push_back(isMuonRec);
        ecalDetId_.push_back((trackID.detIdECAL)());
        hcalDetId_.push_back((trackID.detIdHCAL)());
        ehcalDetId_.push_back((trackID.detIdEHCAL)());
        matchedId_.push_back(tmpmatch); 
        ecal1x1Energy_.push_back(eEcal1x1);
        ecal3x3Energy_.push_back(eEcal3x3);
        ecal5x5Energy_.push_back(eEcal5x5);
        ecal15x15Energy_.push_back(eEcal15x15);
        ecal25x25Energy_.push_back(eEcal25x25);
        hcal1x1Energy_.push_back(eHcal);
        hcal_ieta_.push_back(ieta);
        hcal_iphi_.push_back(iphi);
        o_hcal_iphi_.push_back(o_iphi);
        o_hcal_ieta_.push_back(o_ieta); 
        o_hcalDepth1Energy_.push_back(o_eHcalDepth[0]);
        o_hcalDepth2Energy_.push_back(o_eHcalDepth[1]); edm::LogVerbatim("HBHEMuon") <<"depth 2: "<<o_eHcalDepth[1];
        o_hcalDepth3Energy_.push_back(o_eHcalDepth[2]);
        o_hcalDepth4Energy_.push_back(o_eHcalDepth[3]);
        o_hcalDepth5Energy_.push_back(o_eHcalDepth[4]);
        o_hcalDepth6Energy_.push_back(o_eHcalDepth[5]);
        o_hcalDepth7Energy_.push_back(o_eHcalDepth[6]);

        o_hcalDepth1Energy1_.push_back(o_eHcalDepth1[0]);
        o_hcalDepth2Energy1_.push_back(o_eHcalDepth1[1]);
        o_hcalDepth3Energy1_.push_back(o_eHcalDepth1[2]);
        o_hcalDepth4Energy1_.push_back(o_eHcalDepth1[3]);
        o_hcalDepth5Energy1_.push_back(o_eHcalDepth1[4]);
        o_hcalDepth6Energy1_.push_back(o_eHcalDepth1[5]);
        o_hcalDepth7Energy1_.push_back(o_eHcalDepth1[6]);

             o_hcalDepth1Energy2_.push_back(o_eHcalDepth2[0]);
             o_hcalDepth2Energy2_.push_back(o_eHcalDepth2[1]);
             o_hcalDepth3Energy2_.push_back(o_eHcalDepth2[2]);
             o_hcalDepth4Energy2_.push_back(o_eHcalDepth2[3]);
             o_hcalDepth5Energy2_.push_back(o_eHcalDepth2[4]);
             o_hcalDepth6Energy2_.push_back(o_eHcalDepth2[5]);
             o_hcalDepth7Energy2_.push_back(o_eHcalDepth2[6]);

             o_hcalDepth1Energy3_.push_back(o_eHcalDepth3[0]);
             o_hcalDepth2Energy3_.push_back(o_eHcalDepth3[1]);
             o_hcalDepth3Energy3_.push_back(o_eHcalDepth3[2]);
             o_hcalDepth4Energy3_.push_back(o_eHcalDepth3[3]);
             o_hcalDepth5Energy3_.push_back(o_eHcalDepth3[4]);
             o_hcalDepth6Energy3_.push_back(o_eHcalDepth3[5]);
             o_hcalDepth7Energy3_.push_back(o_eHcalDepth3[6]);

             o_hcalDepth1Energy4_.push_back(o_eHcalDepth4[0]);
             o_hcalDepth2Energy4_.push_back(o_eHcalDepth4[1]);
             o_hcalDepth3Energy4_.push_back(o_eHcalDepth4[2]);
             o_hcalDepth4Energy4_.push_back(o_eHcalDepth4[3]);
             o_hcalDepth5Energy4_.push_back(o_eHcalDepth4[4]);
             o_hcalDepth6Energy4_.push_back(o_eHcalDepth4[5]);
             o_hcalDepth7Energy4_.push_back(o_eHcalDepth4[6]);

             o_hcalDepth1Energy5_.push_back(o_eHcalDepth5[0]);
             o_hcalDepth2Energy5_.push_back(o_eHcalDepth5[1]);
             o_hcalDepth3Energy5_.push_back(o_eHcalDepth5[2]);
             o_hcalDepth4Energy5_.push_back(o_eHcalDepth5[3]);
             o_hcalDepth5Energy5_.push_back(o_eHcalDepth5[4]);
             o_hcalDepth6Energy5_.push_back(o_eHcalDepth5[5]);
             o_hcalDepth7Energy5_.push_back(o_eHcalDepth5[6]);

             o_hcalDepth1Energy6_.push_back(o_eHcalDepth6[0]);
             o_hcalDepth2Energy6_.push_back(o_eHcalDepth6[1]);
             o_hcalDepth3Energy6_.push_back(o_eHcalDepth6[2]);
             o_hcalDepth4Energy6_.push_back(o_eHcalDepth6[3]);
             o_hcalDepth5Energy6_.push_back(o_eHcalDepth6[4]);
             o_hcalDepth6Energy6_.push_back(o_eHcalDepth6[5]);
             o_hcalDepth7Energy6_.push_back(o_eHcalDepth6[6]);

             o_hcalDepth1Energy7_.push_back(o_eHcalDepth7[0]);
             o_hcalDepth2Energy7_.push_back(o_eHcalDepth7[1]);
             o_hcalDepth3Energy7_.push_back(o_eHcalDepth7[2]);
             o_hcalDepth4Energy7_.push_back(o_eHcalDepth7[3]);
             o_hcalDepth5Energy7_.push_back(o_eHcalDepth7[4]);
             o_hcalDepth6Energy7_.push_back(o_eHcalDepth7[5]);
             o_hcalDepth7Energy7_.push_back(o_eHcalDepth7[6]);

             o_hcalDepth1Energy8_.push_back(o_eHcalDepth8[0]);
             o_hcalDepth2Energy8_.push_back(o_eHcalDepth8[1]);
             o_hcalDepth3Energy8_.push_back(o_eHcalDepth8[2]);
             o_hcalDepth4Energy8_.push_back(o_eHcalDepth8[3]);
             o_hcalDepth5Energy8_.push_back(o_eHcalDepth8[4]);
             o_hcalDepth6Energy8_.push_back(o_eHcalDepth8[5]);
             o_hcalDepth7Energy8_.push_back(o_eHcalDepth8[6]);
      //edm::LogVerbatim("HBHEMuon")  <<"iphi push_back: "<<iphi;
      for (int i=0; i<depthMax_; ++i)  {

        //uncalibrated energy
	hcalDepthEnergy_[i].push_back(eHcalDepth[i]);
        //calibrated energy of cells surrounding
        hcalDepthEnergy1_[i].push_back(eHcalDepth1[i]);
        hcalDepthEnergy2_[i].push_back(eHcalDepth2[i]);
        hcalDepthEnergy3_[i].push_back(eHcalDepth3[i]);
        hcalDepthEnergy4_[i].push_back(eHcalDepth4[i]);
        hcalDepthEnergy5_[i].push_back(eHcalDepth5[i]);
        hcalDepthEnergy6_[i].push_back(eHcalDepth6[i]);
        hcalDepthEnergy7_[i].push_back(eHcalDepth7[i]);
        hcalDepthEnergy8_[i].push_back(eHcalDepth8[i]);
        /*
        hcalDepthNumHits_[i].push_back(eHcalNumHits[i]);

        hcalDepthNumHits1_[i].push_back(eHcalNumHits1[i]);
        hcalDepthNumHits2_[i].push_back(eHcalNumHits2[i]);
        hcalDepthNumHits3_[i].push_back(eHcalNumHits3[i]);
        hcalDepthNumHits4_[i].push_back(eHcalNumHits4[i]);
        hcalDepthNumHits5_[i].push_back(eHcalNumHits5[i]);
        hcalDepthNumHits6_[i].push_back(eHcalNumHits6[i]);
        hcalDepthNumHits7_[i].push_back(eHcalNumHits7[i]);
        hcalDepthNumHits8_[i].push_back(eHcalNumHits8[i]);
        */
        //active pathlength
	hcalDepthActiveLength_[i].push_back(activeL[i]);
        //uncalibrated energy of hot cell
	//hcalDepthEnergyHot_[i].push_back(eHcalDepthHot[i]);
        //active pathlength
	hcalDepthActiveLengthHot_[i].push_back(activeHotL[i]);
        //calibrated energy
	hcalDepthEnergyCorr_[i].push_back(eHcalDepthC[i]);
        //calibrated energy of hot cell
	//hcalDepthEnergyHotCorr_[i].push_back(eHcalDepthHotC[i]);
        //charge
	//hcalDepthChargeHot_[i].push_back(cHcalDepthHot[i]);
	//hcalDepthChargeHotBG_[i].push_back(cHcalDepthHotBG[i]);
        //make sure track matches the HCAL cell 
	hcalDepthMatch_[i].push_back(matchDepth[i]);
	//hcalDepthMatchHot_[i].push_back(matchDepthHot[i]);
        hcalHot_.push_back(isHot);
        o_hcalHot_.push_back(o_isHot);
        hcalActiveLength_.push_back(activeLengthTot);
        //hcalActiveLengthHot_.push_back(activeLengthHotTot);
        }//Hcal OK
      } //muon has inner track
      //hcalActiveLength_.push_back(activeLengthTot);
      //hcalHot_.push_back(isHot);
      //o_hcalHot_.push_back(o_isHot);
      //hcalActiveLengthHot_.push_back(activeLengthHotTot);
    } //loop muons
  }//Hcal Hits
  if (o_hcalHot_.size()>0) {
#ifdef EDM_ML_DEBUG
    for (unsigned int i=0; i<hcal_ieta_.size(); ++i)
      edm::LogVerbatim("HBHEMuon") << "[" << i << "] ieta/iphi for entry to "
				   << "HCAL has value of " << hcal_ieta_[i]
				   << ":" << hcal_iphi_[i];
#endif
    tree_->Fill();
  }
  
}

// ------------ method called once each job just before starting event loop  ------------
void HcalHBHEMuonAnalyzer::beginJob() {

  tree_ = fs->make<TTree>("TREE", "TREE");
  tree_->Branch("Event_No",                         &eventNumber_);
  tree_->Branch("Run_No",                           &runNumber_);
  tree_->Branch("LumiNumber",                       &lumiNumber_);
  tree_->Branch("BXNumber",                         &bxNumber_);
  tree_->Branch("GoodVertex",                       &goodVertex_);
  tree_->Branch("LumiScaler",                       &LumiScaler_);
  tree_->Branch("BunchLumi",                        &BunchLumi_);
  tree_->Branch("pileup",                           &pileup_);
  tree_->Branch("pileupRMS",                        &pileupRMS_);
  tree_->Branch("PF_Muon",                          &muon_is_good_);
  tree_->Branch("Global_Muon",                      &muon_global_);
  tree_->Branch("Tracker_muon",                     &muon_tracker_);
  tree_->Branch("MuonIsTight",                      &muon_is_tight_);
  tree_->Branch("charge_of_muon",                   &cGlob_);
  tree_->Branch("pt_of_muon",                       &ptGlob_);
  tree_->Branch("eta_of_muon",                      &etaGlob_);
  tree_->Branch("phi_of_muon",                      &phiGlob_);
  tree_->Branch("energy_of_muon",                   &energyMuon_);
  tree_->Branch("p_of_muon",                        &pMuon_);

  tree_->Branch("IsolationR04",                     &isolationR04_);
  tree_->Branch("ecal_3into3",                      &ecalEnergy_);
  tree_->Branch("hcal_3into3",                      &hcalEnergy_);
  tree_->Branch("tracker_3into3",                   &hoEnergy_);

  tree_->Branch("matchedId",                        &matchedId_);
  tree_->Branch("hcal_cellHot",                     &hcalHot_);
  tree_->Branch("hcal_cellHot_o",                   &o_hcalHot_);
  tree_->Branch("pt_of_genMuon",                    &genMuon_pt_);
  tree_->Branch("eta_of_genMuon",                   &genMuon_eta_);
  tree_->Branch("phi_of_genMuon",                   &genMuon_phi_);
  tree_->Branch("energy_of_genMuon",                &genMuon_energy_);
  tree_->Branch("ecal_1x1",                         &ecal1x1Energy_); 
  tree_->Branch("ecal_3x3",                         &ecal3x3Energy_);
  tree_->Branch("ecal_5x5",                         &ecal5x5Energy_);
  tree_->Branch("ecal_15x15",                       &ecal15x15Energy_);
  tree_->Branch("ecal_25x25",                       &ecal25x25Energy_);
  tree_->Branch("hcal_1x1",                         &hcal1x1Energy_);
  tree_->Branch("ecal_detID",                       &ecalDetId_);
  tree_->Branch("hcal_detID",                       &hcalDetId_);
  tree_->Branch("ehcal_detID",                      &ehcalDetId_);
  tree_->Branch("hcal_ieta",                        &hcal_ieta_);
  tree_->Branch("hcal_iphi",                        &hcal_iphi_);
  tree_->Branch("isMuonRec",                        &isMuonRec_);

  char name[100], namen[100];
  for (int k=0; k<maxDepth_; ++k) {
    sprintf (name, "hcal_edepthuncalibrated%d", (k+1));
    sprintf (namen, "hcal_ndepth%d", (k+1));
    tree_->Branch(name, &hcalDepthEnergy_[k]);
    tree_->Branch(namen, &hcalDepthNumHits_[k]);

    for (int jj=1; jj<9; ++jj) {
         sprintf (name, "hcal_edepth%d%d", jj, (k+1));
         sprintf (namen, "hcal_ndepth%d%d", jj, (k+1));
         if(jj==1){
            tree_->Branch(name, &hcalDepthEnergy1_[k]);
            tree_->Branch(namen, &hcalDepthNumHits1_[k]);
         }
         else if(jj==2){
            tree_->Branch(name, &hcalDepthEnergy2_[k]);
            tree_->Branch(namen, &hcalDepthNumHits2_[k]);
         }
         else if(jj==3){
            tree_->Branch(name, &hcalDepthEnergy3_[k]);
            tree_->Branch(namen, &hcalDepthNumHits3_[k]);
         }
         else if(jj==4){
            tree_->Branch(name, &hcalDepthEnergy4_[k]);
            tree_->Branch(namen, &hcalDepthNumHits4_[k]);
         }
         else if(jj==5){
            tree_->Branch(name, &hcalDepthEnergy5_[k]);
            tree_->Branch(namen, &hcalDepthNumHits5_[k]);
         }
         else if(jj==6){
            tree_->Branch(name, &hcalDepthEnergy6_[k]);
            tree_->Branch(namen, &hcalDepthNumHits6_[k]);
         }
         else if(jj==7){
            tree_->Branch(name, &hcalDepthEnergy7_[k]);
            tree_->Branch(namen, &hcalDepthNumHits7_[k]);
         }
         else{
            tree_->Branch(name, &hcalDepthEnergy8_[k]);
            tree_->Branch(namen, &hcalDepthNumHits8_[k]);
         }
    }

    sprintf (name, "hcal_activeL%d", (k+1));
    tree_->Branch(name,  &hcalDepthActiveLength_[k]);
    //sprintf (name, "hcal_edepthHot%d", (k+1));
    //tree_->Branch(name,  &hcalDepthEnergyHot_[k]);
    //sprintf (name, "hcal_activeHotL%d", (k+1));
    //tree_->Branch(name, &hcalDepthActiveLengthHot_[k]);
    //sprintf (name, "hcal_cdepthHot%d", (k+1));
    //tree_->Branch(name,  &hcalDepthChargeHot_[k]);
    //sprintf (name, "hcal_cdepthHotBG%d", (k+1));
    //tree_->Branch(name,  &hcalDepthChargeHotBG_[k]);
    sprintf (name, "hcal_edepth%d", (k+1));
    tree_->Branch(name, &hcalDepthEnergyCorr_[k]);
    //sprintf (name, "hcal_edepthHotCorrect%d", (k+1));
    //tree_->Branch(name,  &hcalDepthEnergyHotCorr_[k]);
    sprintf (name, "hcal_depthMatch%d", (k+1));
    tree_->Branch(name,  &hcalDepthMatch_[k]);
    //sprintf (name, "hcal_depthMatchHot%d", (k+1));
    //tree_->Branch(name,  &hcalDepthMatchHot_[k]);
  }
  
  tree_->Branch("activeLength",                     &hcalActiveLength_);
  //tree_->Branch("activeLengthHot",                  &hcalActiveLengthHot_);
 
  tree_->Branch("eEcal",             &eEcal_);
  tree_->Branch("Ecal_label",        &Ecal_label_);
  tree_->Branch("hcal_ieta_o",                        &o_hcal_ieta_);
  tree_->Branch("hcal_iphi_o",                        &o_hcal_iphi_);
  tree_->Branch("hcal_edepth1_o",     &o_hcalDepth1Energy_);
  tree_->Branch("hcal_edepth2_o",     &o_hcalDepth2Energy_);
  tree_->Branch("hcal_edepth3_o",     &o_hcalDepth3Energy_);
  tree_->Branch("hcal_edepth4_o",     &o_hcalDepth4Energy_);
  //1
  tree_->Branch("hcal_edepth11_o",     &o_hcalDepth1Energy1_);
  tree_->Branch("hcal_edepth12_o",     &o_hcalDepth2Energy1_);
  tree_->Branch("hcal_edepth13_o",     &o_hcalDepth3Energy1_);
  tree_->Branch("hcal_edepth14_o",     &o_hcalDepth4Energy1_);
  //2
  tree_->Branch("hcal_edepth21_o",     &o_hcalDepth1Energy2_);
  tree_->Branch("hcal_edepth22_o",     &o_hcalDepth2Energy2_);
  tree_->Branch("hcal_edepth23_o",     &o_hcalDepth3Energy2_);
  tree_->Branch("hcal_edepth24_o",     &o_hcalDepth4Energy2_);
  //3
  tree_->Branch("hcal_edepth31_o",     &o_hcalDepth1Energy3_);
  tree_->Branch("hcal_edepth32_o",     &o_hcalDepth2Energy3_);
  tree_->Branch("hcal_edepth33_o",     &o_hcalDepth3Energy3_);
  tree_->Branch("hcal_edepth34_o",     &o_hcalDepth4Energy3_);
  //4
  tree_->Branch("hcal_edepth41_o",     &o_hcalDepth1Energy4_);
  tree_->Branch("hcal_edepth42_o",     &o_hcalDepth2Energy4_);
  tree_->Branch("hcal_edepth43_o",     &o_hcalDepth3Energy4_);
  tree_->Branch("hcal_edepth44_o",     &o_hcalDepth4Energy4_);
  //5
  tree_->Branch("hcal_edepth51_o",     &o_hcalDepth1Energy5_);
  tree_->Branch("hcal_edepth52_o",     &o_hcalDepth2Energy5_);
  tree_->Branch("hcal_edepth53_o",     &o_hcalDepth3Energy5_);
  tree_->Branch("hcal_edepth54_o",     &o_hcalDepth4Energy5_);
  //6
  tree_->Branch("hcal_edepth61_o",     &o_hcalDepth1Energy6_);
  tree_->Branch("hcal_edepth62_o",     &o_hcalDepth2Energy6_);
  tree_->Branch("hcal_edepth63_o",     &o_hcalDepth3Energy6_);
  tree_->Branch("hcal_edepth64_o",     &o_hcalDepth4Energy6_);
  //7
  tree_->Branch("hcal_edepth71_o",     &o_hcalDepth1Energy7_);
  tree_->Branch("hcal_edepth72_o",     &o_hcalDepth2Energy7_);
  tree_->Branch("hcal_edepth73_o",     &o_hcalDepth3Energy7_);
  tree_->Branch("hcal_edepth74_o",     &o_hcalDepth4Energy7_);
  //8
  tree_->Branch("hcal_edepth81_o",     &o_hcalDepth1Energy8_);
  tree_->Branch("hcal_edepth82_o",     &o_hcalDepth2Energy8_);
  tree_->Branch("hcal_edepth83_o",     &o_hcalDepth3Energy8_);
  tree_->Branch("hcal_edepth84_o",     &o_hcalDepth4Energy8_);
  if (maxDepth_ > 4) {
    tree_->Branch("hcal_edepth5_o",     &o_hcalDepth5Energy_);
    tree_->Branch("hcal_edepth15_o",    &o_hcalDepth5Energy1_);
    tree_->Branch("hcal_edepth25_o",    &o_hcalDepth5Energy2_);
    tree_->Branch("hcal_edepth35_o",    &o_hcalDepth5Energy3_);
    tree_->Branch("hcal_edepth45_o",    &o_hcalDepth5Energy4_);
    tree_->Branch("hcal_edepth55_o",    &o_hcalDepth5Energy5_);
    tree_->Branch("hcal_edepth65_o",    &o_hcalDepth5Energy6_);
    tree_->Branch("hcal_edepth75_o",    &o_hcalDepth5Energy7_);
    tree_->Branch("hcal_edepth85_o",    &o_hcalDepth5Energy8_);

    if (maxDepth_ > 5) {
      tree_->Branch("hcal_edepth6_o",     &o_hcalDepth6Energy_);
      tree_->Branch("hcal_edepth16_o",    &o_hcalDepth6Energy1_);
      tree_->Branch("hcal_edepth26_o",    &o_hcalDepth6Energy2_);
      tree_->Branch("hcal_edepth36_o",    &o_hcalDepth6Energy3_);
      tree_->Branch("hcal_edepth46_o",    &o_hcalDepth6Energy4_);
      tree_->Branch("hcal_edepth56_o",    &o_hcalDepth6Energy5_);
      tree_->Branch("hcal_edepth66_o",    &o_hcalDepth6Energy6_);
      tree_->Branch("hcal_edepth76_o",    &o_hcalDepth6Energy7_);
      tree_->Branch("hcal_edepth86_o",    &o_hcalDepth6Energy8_);

      if (maxDepth_ > 6) {
        tree_->Branch("hcal_edepth7_o",     &o_hcalDepth7Energy_);
        tree_->Branch("hcal_edepth17_o",    &o_hcalDepth7Energy1_);
        tree_->Branch("hcal_edepth27_o",    &o_hcalDepth7Energy2_);
        tree_->Branch("hcal_edepth37_o",    &o_hcalDepth7Energy3_);
        tree_->Branch("hcal_edepth47_o",    &o_hcalDepth7Energy4_);
        tree_->Branch("hcal_edepth57_o",    &o_hcalDepth7Energy5_);
        tree_->Branch("hcal_edepth67_o",    &o_hcalDepth7Energy6_);
        tree_->Branch("hcal_edepth77_o",    &o_hcalDepth7Energy7_);
        tree_->Branch("hcal_edepth87_o",    &o_hcalDepth7Energy8_);

      }
    }
  }

  tree_->Branch("hltresults",                       &hltresults_);
  tree_->Branch("all_triggers",                     &all_triggers_);
  tree_->Branch("h_cutflow","TH1F",&h_cutflow,32000,0);
  
}

// ------------ method called when starting to processes a run  ------------
void HcalHBHEMuonAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  edm::ESHandle<HcalDDDRecConstants> pHRNDC;
  iSetup.get<HcalRecNumberingRecord>().get(pHRNDC);
  hdc_ = pHRNDC.product();
  actHB.clear();
  actHE.clear();
  actHB = hdc_->getThickActive(0);
  actHE = hdc_->getThickActive(1);
#ifdef EDM_ML_DEBUG
  unsigned int k1(0), k2(0);
  edm::LogVerbatim("HBHEMuon") << actHB.size() << " Active Length for HB";
  for (const auto& act : actHB) {
    edm::LogVerbatim("HBHEMuon") << "[" << k1 << "] ieta " << act.ieta
				 << " depth " << act.depth << " zside "
				 << act.zside << " type " << act.stype
				 << " phi " << act.iphis.size() << ":"
				 << act.iphis[0] << " L " << act.thick;
    HcalDetId hcid1(HcalBarrel,(act.ieta)*(act.zside),act.iphis[0],act.depth);
    HcalDetId hcid2 = mergedDepth_ ? hdc_->mergedDepthDetId(hcid1) : hcid1;
    edm::LogVerbatim("HBHEMuon") << hcid1 << " | " << hcid2 << " L "
				 << activeLength(DetId(hcid2));
    ++k1;
  }
  edm::LogVerbatim("HBHEMuon") << actHE.size() << " Active Length for HE";
  for (const auto& act : actHE) {
    edm::LogVerbatim("HBHEMuon") << "[" << k2 << "] ieta " << act.ieta
				 << " depth " << act.depth << " zside "
				 << act.zside << " type " << act.stype
				 << " phi " << act.iphis.size() << ":"
				 << act.iphis[0] << " L " << act.thick;
    HcalDetId hcid1(HcalEndcap,(act.ieta)*(act.zside),act.iphis[0],act.depth);
    HcalDetId hcid2 = mergedDepth_ ? hdc_->mergedDepthDetId(hcid1) : hcid1;
    edm::LogVerbatim("HBHEMuon") << hcid1 << " | " << hcid2 << " L "
				 << activeLength(DetId(hcid2));
    ++k2;
  }
#endif

  bool changed = true;
  all_triggers_.clear();
  if (hltConfig_.init(iRun, iSetup, "HLT" , changed)) {
    // if init returns TRUE, initialisation has succeeded!
#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("HBHEMuon") << "HLT config with process name " 
				 << "HLT" << " successfully extracted"
				 << std::endl;
#endif
    unsigned int ntriggers = hltConfig_.size();
    for (unsigned int t=0;t<ntriggers;++t) {
      std::string hltname(hltConfig_.triggerName(t));
      for (unsigned int ik=0; ik<6; ++ik) {
	if (hltname.find(triggers_[ik])!=std::string::npos ){
	  all_triggers_.push_back(hltname);
	  break;
	}
      }
    }//loop over ntriggers
    edm::LogVerbatim("HBHEMuon") << "All triggers size in begin run " 
				 << all_triggers_.size() << std::endl;
  } else {
    edm::LogError("HBHEMuon") << "Error! HLT config extraction with process "
			      << "name HLT failed";
  }

  edm::ESHandle<HcalTopology> htopo;
  iSetup.get<HcalRecNumberingRecord>().get(htopo);
  theHBHETopology_ = htopo.product();

  edm::ESHandle<HcalRespCorrs> resp;
  iSetup.get<HcalRespCorrsRcd>().get(resp);
  respCorrs_ = new HcalRespCorrs(*resp.product());
  respCorrs_->setTopo(theHBHETopology_);

  // Write correction factors for all HB/HE events
  if (writeRespCorr_) {
    edm::ESHandle<CaloGeometry> pG;
    iSetup.get<CaloGeometryRecord>().get(pG);
    const CaloGeometry* geo = pG.product();
    const HcalGeometry* gHcal = (const HcalGeometry*)(geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel));
    const std::vector<DetId>& ids = gHcal->getValidDetIds(DetId::Hcal,0);
    edm::LogVerbatim("HBHEMuon") << "\nTable of Correction Factors for Run "
				 << iRun.run() << "\n";
    for (auto const& id: ids) {
      if ((id.det() == DetId::Hcal) && 
	  ((id.subdetId() == HcalBarrel) || (id.subdetId() == HcalEndcap))) {
	edm::LogVerbatim("HBHEMuon") << HcalDetId(id) << " " << id.rawId() <<" "
				     << (respCorrs_->getValues(id))->getValue();
      }
    }
  }
  h_cutflow->Write();
}
void HcalHBHEMuonAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  h_cutflow->Write();
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HcalHBHEMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("hlTriggerResults",edm::InputTag("TriggerResults","","HLT"));
  desc.add<edm::InputTag>("labelEBRecHit",edm::InputTag("ecalRecHit","EcalRecHitsEB"));
  desc.add<edm::InputTag>("labelEERecHit",edm::InputTag("ecalRecHit","EcalRecHitsEE"));
  desc.add<edm::InputTag>("labelHBHERecHit",edm::InputTag("hbhereco"));
  desc.add<edm::InputTag>("labelLumiScalers", edm::InputTag("scalersRawToDigi"));
  desc.add<std::string>("labelVertex","offlinePrimaryVertices");
  desc.add<std::string>("labelMuon","muons");
  //desc.add<std::string>("labelGenPart","genParticles");
  std::vector<std::string> trig = {"HLT_IsoMu17","HLT_IsoMu20",
				   "HLT_IsoMu24","HLT_IsoMu27",
				   "HLT_Mu45","HLT_Mu50"};
  desc.add<std::vector<std::string>>("triggers",trig);
  desc.addUntracked<int>("verbosity",0);
  desc.add<int>("useRaw",0);
  desc.add<bool>("unCorrect",false);
  desc.add<bool>("getCharge",false);
  desc.add<bool>("collapseDepth",false);
  desc.add<bool>("isItPlan1",false);
  desc.addUntracked<bool>("ignoreHECorr",false);
  desc.addUntracked<bool>("isItPreRecHit",false);
  desc.addUntracked<std::string>("moduleName","");
  desc.addUntracked<std::string>("processName","");
  desc.addUntracked<int>("maxDepth",4);
  desc.addUntracked<std::string>("fileInCorr","");
  desc.addUntracked<bool>("writeRespCorr",false);
  descriptions.add("hcalHBHEMuon",desc);
}

void HcalHBHEMuonAnalyzer::clearVectors() {
  ///clearing vectots
  eventNumber_ = -99999;
  runNumber_   = -99999;
  lumiNumber_  = -99999;
  bxNumber_    = -99999;
  goodVertex_  = -99999;
  LumiScaler_ = -99999;
  BunchLumi_ = -99999;

  muon_is_good_.clear();
  muon_global_.clear();
  muon_tracker_.clear();
  cGlob_.clear();
  ptGlob_.clear();
  etaGlob_.clear(); 
  phiGlob_.clear(); 
  energyMuon_.clear();
  pMuon_.clear();
  muon_is_tight_.clear();

  isolationR04_.clear();
  ecalEnergy_.clear();
  hcalEnergy_.clear();
  hoEnergy_.clear();
  matchedId_.clear();
  hcalHot_.clear();
  o_hcalHot_.clear();
  
  genMuon_pt_.clear();
  genMuon_eta_.clear();
  genMuon_phi_.clear();
  genMuon_energy_.clear();

  ecal1x1Energy_.clear();
  ecal3x3Energy_.clear();
  ecal5x5Energy_.clear();
  ecal15x15Energy_.clear();
  ecal25x25Energy_.clear();
  hcal1x1Energy_.clear();
  ecalDetId_.clear();
  hcalDetId_.clear();
  ehcalDetId_.clear();
  hcal_ieta_.clear();
  hcal_iphi_.clear();
  o_hcal_ieta_.clear();
  o_hcal_iphi_.clear();
  eEcal_.clear();
  Ecal_label_.clear();
  isMuonRec_.clear();
 
  o_hcalDepth1Energy_.clear();
  o_hcalDepth1Energy1_.clear();
  o_hcalDepth1Energy2_.clear();
  o_hcalDepth1Energy3_.clear();
  o_hcalDepth1Energy4_.clear();
  o_hcalDepth1Energy5_.clear();
  o_hcalDepth1Energy6_.clear();
  o_hcalDepth1Energy7_.clear();
  o_hcalDepth1Energy8_.clear();

  o_hcalDepth2Energy_.clear(); 
  o_hcalDepth2Energy1_.clear();
  o_hcalDepth2Energy2_.clear();
  o_hcalDepth2Energy3_.clear();
  o_hcalDepth2Energy4_.clear();
  o_hcalDepth2Energy5_.clear();
  o_hcalDepth2Energy6_.clear();
  o_hcalDepth2Energy7_.clear();
  o_hcalDepth2Energy8_.clear();

  o_hcalDepth3Energy_.clear();
  o_hcalDepth3Energy1_.clear();
  o_hcalDepth3Energy2_.clear();
  o_hcalDepth3Energy3_.clear();
  o_hcalDepth3Energy4_.clear();
  o_hcalDepth3Energy5_.clear();
  o_hcalDepth3Energy6_.clear();
  o_hcalDepth3Energy7_.clear();
  o_hcalDepth3Energy8_.clear();

  o_hcalDepth4Energy_.clear();
  o_hcalDepth4Energy1_.clear();
  o_hcalDepth4Energy2_.clear();
  o_hcalDepth4Energy3_.clear();
  o_hcalDepth4Energy4_.clear();
  o_hcalDepth4Energy5_.clear();
  o_hcalDepth4Energy6_.clear();
  o_hcalDepth4Energy7_.clear();
  o_hcalDepth4Energy8_.clear();

  o_hcalDepth5Energy_.clear();
  o_hcalDepth5Energy1_.clear();
  o_hcalDepth5Energy2_.clear();
  o_hcalDepth5Energy3_.clear();
  o_hcalDepth5Energy4_.clear();
  o_hcalDepth5Energy5_.clear();
  o_hcalDepth5Energy6_.clear();
  o_hcalDepth5Energy7_.clear();
  o_hcalDepth5Energy8_.clear();

  o_hcalDepth6Energy_.clear();
  o_hcalDepth6Energy1_.clear();
  o_hcalDepth6Energy2_.clear();
  o_hcalDepth6Energy3_.clear();
  o_hcalDepth6Energy4_.clear();
  o_hcalDepth6Energy5_.clear();
  o_hcalDepth6Energy6_.clear();
  o_hcalDepth6Energy7_.clear();
  o_hcalDepth6Energy8_.clear();

  o_hcalDepth7Energy_.clear();
  o_hcalDepth7Energy1_.clear();
  o_hcalDepth7Energy2_.clear();
  o_hcalDepth7Energy3_.clear();
  o_hcalDepth7Energy4_.clear();
  o_hcalDepth7Energy5_.clear();
  o_hcalDepth7Energy6_.clear();
  o_hcalDepth7Energy7_.clear();
  o_hcalDepth7Energy8_.clear();
  
  for (int i=0; i<maxDepth_; ++i) {
    hcalDepthEnergy_[i].clear();
    hcalDepthEnergy1_[i].clear();
    hcalDepthEnergy2_[i].clear();
    hcalDepthEnergy3_[i].clear();
    hcalDepthEnergy4_[i].clear();
    hcalDepthEnergy5_[i].clear();
    hcalDepthEnergy6_[i].clear();
    hcalDepthEnergy7_[i].clear();
    hcalDepthEnergy8_[i].clear();
    hcalDepthNumHits_[i].clear();
    hcalDepthNumHits1_[i].clear();
    hcalDepthNumHits2_[i].clear();
    hcalDepthNumHits3_[i].clear();
    hcalDepthNumHits4_[i].clear();
    hcalDepthNumHits5_[i].clear();
    hcalDepthNumHits6_[i].clear();
    hcalDepthNumHits7_[i].clear();
    hcalDepthNumHits8_[i].clear();
    hcalDepthActiveLength_[i].clear();
    hcalDepthEnergyHot_[i].clear();
    hcalDepthActiveLengthHot_[i].clear();
    hcalDepthChargeHot_[i].clear();
    hcalDepthChargeHotBG_[i].clear();
    hcalDepthEnergyCorr_[i].clear();
    hcalDepthEnergyHotCorr_[i].clear();
    hcalDepthMatch_[i].clear();
    hcalDepthMatchHot_[i].clear();
  }
  hcalActiveLength_.clear();
  //hcalActiveLengthHot_.clear();
  hltresults_.clear();
}

int HcalHBHEMuonAnalyzer::matchId(const HcalDetId& id1, const HcalDetId& id2) {

  HcalDetId kd1(id1.subdet(),id1.ieta(),id1.iphi(),1);
  HcalDetId kd2(id1.subdet(),id2.ieta(),id2.iphi(),1);
  int match = ((kd1 == kd2) ? 1 : 0);
  return match;
}

double HcalHBHEMuonAnalyzer::activeLength(const DetId& id_) {
  HcalDetId id(id_);
  int ieta = id.ietaAbs();
  int zside= id.zside();
  int iphi = id.iphi();
  std::vector<int> dpths;
  if (mergedDepth_) {
    std::vector<HcalDetId> ids;
    hdc_->unmergeDepthDetId(id,ids);
    for (auto idh : ids) 
      dpths.emplace_back(idh.depth());
  } else {
    dpths.emplace_back(id.depth());
  }
  double lx(0);
  if (id.subdet() == HcalBarrel) {
    for (unsigned int i=0; i<actHB.size(); ++i) {
      if ((ieta == actHB[i].ieta) && (zside == actHB[i].zside) && 
	  (std::find(dpths.begin(),dpths.end(),actHB[i].depth) != dpths.end())&&
	  (std::find(actHB[i].iphis.begin(),actHB[i].iphis.end(),iphi) !=
	   actHB[i].iphis.end())) {
	lx += actHB[i].thick;
      }
    }
  } else {
    for (unsigned int i=0; i<actHE.size(); ++i) {
      if ((ieta == actHE[i].ieta) && (zside == actHE[i].zside) && 
	  (std::find(dpths.begin(),dpths.end(),actHE[i].depth) != dpths.end())&&
	  (std::find(actHE[i].iphis.begin(),actHE[i].iphis.end(),iphi) !=
	   actHE[i].iphis.end())) {
	lx += actHE[i].thick;
      }
    }
  }
  return lx;
}

bool HcalHBHEMuonAnalyzer::isGoodVertex(const reco::Vertex& vtx) {
  if (vtx.isFake())                   return false;
  if (vtx.ndof() < 4)                 return false;
  if (vtx.position().Rho() > 2.)      return false;
  if (fabs(vtx.position().Z()) > 24.) return false;
  return true;
}

double HcalHBHEMuonAnalyzer::respCorr(const DetId& id) {
  double cfac(1.0);
  if (useMyCorr_) {
    auto itr = corrValue_.find(id);
    if (itr != corrValue_.end()) cfac = itr->second;
  } else if (respCorrs_ != nullptr) {
    cfac = (respCorrs_->getValues(id))->getValue();
  }
  return cfac;
}

double HcalHBHEMuonAnalyzer::gainFactor(const edm::ESHandle<HcalDbService>& conditions, const HcalDetId& id) {
  double gain(0.0);
  const HcalCalibrations& calibs=conditions->getHcalCalibrations(id);
  for (int capid=0; capid<4; ++capid) 
    gain += (0.25*calibs.respcorrgain(capid));
  return gain;
}

int HcalHBHEMuonAnalyzer::depth16HE(int ieta, int iphi) {
  // Transition between HB/HE is special 
  // For Run 1 or for Plan1 standard reconstruction it is 3
  // For runs beyond 2018 or in Plan1 for HEP17 it is 4
  int zside = (ieta > 0) ? 1 : -1;
  int depth = theHBHETopology_->dddConstants()->getMinDepth(1,16,iphi,zside);
  if (isItPlan1_ && (!isItPreRecHit_)) depth = 3;
#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("HBHEMuon") << "Plan1 " << isItPlan1_ << " PreRecHit " 
			       << isItPreRecHit_ << " phi " << iphi
			       << " depth " << depth;
#endif
  return depth;
}

bool HcalHBHEMuonAnalyzer::goodCell(const HcalDetId& hcid, 
				    const reco::Track* pTrack, 
				    const CaloGeometry* geo, 
				    const MagneticField* bField) {

  std::pair<double,double> rz = hdc_->getRZ(hcid);
  bool typeRZ = (hcid.subdet() == HcalEndcap) ? false : true;
  bool match = spr::propagateHCAL(pTrack, geo, bField, typeRZ, rz, (((verbosity_/10000)%10)>0));
  return match;
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(HcalHBHEMuonAnalyzer);
