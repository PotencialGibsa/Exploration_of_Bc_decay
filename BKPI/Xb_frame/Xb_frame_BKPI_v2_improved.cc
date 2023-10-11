#include <memory>

/// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/FWLite/interface/EventBase.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
/// triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
/// tracks
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

/// muons
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

/// vertex fits
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//// gen ?
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/Error.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TLorentzVector.h"
#include <utility>
#include <string>
#include <map>
#include "Particles_parameters.h"

ParticleMass PM_PDG_MUON_MASS = PDG_MUON_MASS;
ParticleMass PM_PDG_JPSI_MASS = PDG_JPSI_MASS;
ParticleMass PM_PDG_KAON_MASS = PDG_KAON_MASS;
ParticleMass PM_PDG_PION_MASS = PDG_PION_MASS;
ParticleMass PM_PDG_PROTON_MASS = PDG_PROTON_MASS;
ParticleMass PM_PDG_PSI2S_MASS = PDG_PSI2S_MASS;

using namespace edm;
using namespace std;
using namespace reco;

class Xb_frame : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

    public:
        explicit Xb_frame(const edm::ParameterSet&);
        ~Xb_frame();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

	// ----------member data ---------------------------
	//-----------tokens---------------------------------

	edm::EDGetTokenT<edm::TriggerResults> hlTriggerResults_;
	edm::EDGetTokenT<trigger::TriggerEvent> TriggerEventToken_;
	std::vector<std::string> TriggerFilters_;
	edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
	edm::EDGetTokenT<reco::VertexCollection> vtxSample;
	edm::EDGetTokenT<vector < pat::GenericParticle > > tracks_;
	edm::EDGetTokenT< vector < pat::Muon >> muons_;
	edm::EDGetTokenT<reco::TrackCollection> revtxtrks_;

	std::vector<float> *BK_mass;
	std::vector<float> *BK_px            , *BK_py           , *BK_pz;
	std::vector<float> *BK_DecayVtxX     , *BK_DecayVtxY    , *BK_DecayVtxZ;
	std::vector<float> *BK_DecayVtxXE    , *BK_DecayVtxYE   , *BK_DecayVtxZE;
	std::vector<float> *BK_Prob;

	std::vector<float> *B_mass          , *B_mass_c0;
	std::vector<float> *B_px            , *B_py           , *B_pz;
	std::vector<float> *B_px_CV , *B_py_CV  , *B_pz_CV;
	std::vector<float> *B_DecayVtxX     , *B_DecayVtxY    , *B_DecayVtxZ;
	std::vector<float> *B_DecayVtxXE    , *B_DecayVtxYE   , *B_DecayVtxZE;
	std::vector<float> *B_Prob          , *B_J_Prob;

	std::vector<float> *B_J_mass;
	std::vector<float> *B_J_px          , *B_J_py         , *B_J_pz;
	std::vector<float> *B_J_DecayVtxX   , *B_J_DecayVtxY  , *B_J_DecayVtxZ;
	std::vector<float> *B_J_DecayVtxXE  , *B_J_DecayVtxYE , *B_J_DecayVtxZE;

	std::vector<float> *B_mu_px1_cjp    , *B_mu_py1_cjp   , *B_mu_pz1_cjp;
	std::vector<float> *B_mu_px1        , *B_mu_py1       , *B_mu_pz1;
	std::vector<float> *B_mu_px2_cjp    , *B_mu_py2_cjp   , *B_mu_pz2_cjp;
	std::vector<float> *B_mu_px2        , *B_mu_py2       , *B_mu_pz2;

	std::vector<float> *K1_px    , *K1_py     , *K1_pz;
	std::vector<float> *K1_px_CV , *K1_py_CV  , *K1_pz_CV;
	std::vector<int> *K1_ch;
	std::vector<float> *K1_vx    , *K1_vy     , *K1_vz;
	std::vector<float> *K1_ips;
	std::vector<float> *KA1_drTRG , *KA1_dptTRG;

	std::vector<float> *K2_px    , *K2_py     , *K2_pz;
	std::vector<float> *K2_px_CV , *K2_py_CV  , *K2_pz_CV;
	std::vector<int> *K2_ch;
	std::vector<float> *K2_vx    , *K2_vy     , *K2_vz;
	std::vector<float> *K2_ips;

	std::vector<float> *PI3_px    , *PI3_py     , *PI3_pz;
	std::vector<float> *PI3_px_CV , *PI3_py_CV  , *PI3_pz_CV;
	std::vector<int> *PI3_ch;
	std::vector<float> *PI3_vx    , *PI3_vy     , *PI3_vz;
	std::vector<float> *PI3_ips;
	//
	std::vector<float> *PV_becos_XX , *PV_becos_YY  , *PV_becos_ZZ;
	std::vector<float> *PV_becos_EX , *PV_becos_EY  , *PV_becos_EZ;
	std::vector<float> *PV_becos_CL;
	std::vector<int>   *PV_becos_dN;

	std::vector<bool>  *trig0_fi;
	std::vector<bool>  *trig1_fi;
	std::vector<bool>  *trig2_fi;
	std::vector<bool>  *trig3_fi;
	std::vector<bool>  *trig4_fi;
	std::vector<bool>  *trig5_fi;
	std::vector<bool>  *trig6_fi;
	std::vector<bool>  *trig7_fi;
	std::vector<bool>  *trig8_fi;
	std::vector<bool>  *trig9_fi;
	std::vector<bool>  *trig10_fi;
	std::vector<bool>  *trig11_fi;
	std::vector<bool>  *trig12_fi;
	std::vector<bool>  *trig13_fi;

	std::vector<bool>  *trig0_ma;
	std::vector<bool>  *trig1_ma;
	std::vector<bool>  *trig2_ma;
	std::vector<bool>  *trig3_ma;
	std::vector<bool>  *trig4_ma;
	std::vector<bool>  *trig5_ma;
	std::vector<bool>  *trig6_ma;
	std::vector<bool>  *trig7_ma;
	std::vector<bool>  *trig8_ma;
	std::vector<bool>  *trig9_ma;
	std::vector<bool>  *trig10_ma;
	std::vector<bool>  *trig11_ma;
	std::vector<bool>  *trig12_ma;
	std::vector<bool>  *trig13_ma;


	Int_t       nCand;

	Int_t       run;
	Int_t       event;
	Float_t     lumi;

	Int_t       numPV;
	Int_t       numTrack;
	Int_t       numV0;

	TTree *wwtree;
	TFile *f;
	std::string fileName;

};

//
Xb_frame::Xb_frame(const edm::ParameterSet& iConfig) :
	hlTriggerResults_( consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTriggerResults"))),
	TriggerEventToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerEvent"))),
	TriggerFilters_(   iConfig.getParameter<std::vector<std::string>>("TriggerFilters")),
	thebeamspot_(      consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
	vtxSample(         consumes< reco::VertexCollection > (iConfig.getParameter<edm::InputTag>("VtxSample"))),
	tracks_(           consumes< vector < pat::GenericParticle > > (iConfig.getParameter<edm::InputTag>("Trak"))),
	muons_(            consumes< vector <pat::Muon >> (iConfig.getParameter<edm::InputTag>("muons"))),
	revtxtrks_(        consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("revtxtrks"))),

	BK_mass(0),
	BK_px(0),            BK_py(0),            BK_pz(0),
	BK_DecayVtxX(0),     BK_DecayVtxY(0),     BK_DecayVtxZ(0),
	BK_DecayVtxXE(0),    BK_DecayVtxYE(0),    BK_DecayVtxZE(0),
	BK_Prob(0),

	B_mass(0),          B_mass_c0(0),
	B_px(0),            B_py(0),            B_pz(0),
	B_px_CV(0)  ,       B_py_CV(0),         B_pz_CV(0),
	B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
	B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
	B_Prob(0),

	B_J_Prob(0),        B_J_mass(0),
	B_J_px(0),          B_J_py(0),          B_J_pz(0),
	B_J_DecayVtxX(0),   B_J_DecayVtxY(0),   B_J_DecayVtxZ(0),
	B_J_DecayVtxXE(0),  B_J_DecayVtxYE(0),  B_J_DecayVtxZE(0),

	B_mu_px1_cjp(0),    B_mu_py1_cjp(0),    B_mu_pz1_cjp(0),
	B_mu_px1(0),        B_mu_py1(0),        B_mu_pz1(0),

	B_mu_px2_cjp(0),    B_mu_py2_cjp(0),    B_mu_pz2_cjp(0),
	B_mu_px2(0),        B_mu_py2(0),        B_mu_pz2(0),

	K1_px(0)     , K1_py(0)   , K1_pz(0),
	K1_px_CV(0)  , K1_py_CV(0), K1_pz_CV(0),
	K1_ch(0),
	K1_vx(0)     , K1_vy(0)   , K1_vz(0),
	K1_ips(0),
	KA1_drTRG(0) , KA1_dptTRG(0),

	K2_px(0)     , K2_py(0)   , K2_pz(0),
	K2_px_CV(0)  , K2_py_CV(0), K2_pz_CV(0),
	K2_ch(0),
	K2_vx(0)     , K2_vy(0)   , K2_vz(0),
	K2_ips(0),

	PI3_px(0)     , PI3_py(0)   , PI3_pz(0),
	PI3_px_CV(0)  , PI3_py_CV(0), PI3_pz_CV(0),
	PI3_ch(0),
	PI3_vx(0)     , PI3_vy(0)   , PI3_vz(0),
	PI3_ips(0),

	PV_becos_XX(0)  , PV_becos_YY(0), PV_becos_ZZ(0),
	PV_becos_EX(0)  , PV_becos_EY(0), PV_becos_EZ(0),
	PV_becos_CL(0)  , PV_becos_dN(0),

	trig0_fi(0),    trig1_fi(0),    trig2_fi(0),
	trig3_fi(0),    trig4_fi(0),    trig5_fi(0),
	trig6_fi(0),    trig7_fi(0),    trig8_fi(0),
	trig9_fi(0),    trig10_fi(0),   trig11_fi(0),
	trig12_fi(0),   trig13_fi(0),

	trig0_ma(0),    trig1_ma(0),    trig2_ma(0),
	trig3_ma(0),    trig4_ma(0),    trig5_ma(0),
	trig6_ma(0),    trig7_ma(0),    trig8_ma(0),
	trig9_ma(0),    trig10_ma(0),   trig11_ma(0),
	trig12_ma(0),   trig13_ma(0),


	//
	nCand(0),

	run(0),
	event(0),
	lumi(0),

	numPV(0),
	numTrack(0),
	numV0(0)
{
	fileName = iConfig.getUntrackedParameter<std::string>("fileName","BFinder.root");
   	//now do what ever initialization is needed
   	usesResource("TFileService");
}

Xb_frame::~Xb_frame()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

// ------------ method called for each event  ------------
void
Xb_frame::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;
	using namespace std;
	using reco::MuonCollection;

	run   = iEvent.id().run();
	event = iEvent.id().event();

	lumi = 1;
	lumi = iEvent.luminosityBlock();

	ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

	// Get HLT results HLT   HLT  HLT   HLT  HLT   HLT
	edm::Handle<edm::TriggerResults> triggerResults_handle;
	iEvent.getByToken(hlTriggerResults_, triggerResults_handle);

	edm::Handle<trigger::TriggerEvent> triggerEvent;
	iEvent.getByToken(TriggerEventToken_, triggerEvent);
	auto triggerObjects = (triggerEvent->getObjects());
	// PAT trigger helper for trigger matching information
	// const pat::helper::TriggerMatchHelper matchHelper;

	//// PV PV  PV  PV  PV  PV  PV  PV  PV  PV  PV  PV  PV
	Handle < VertexCollection > recVtxs ;
	iEvent.getByToken(vtxSample, recVtxs);
	///
	Vertex thePrimaryV;
	thePrimaryV = Vertex(*(recVtxs->begin()));
	const reco::VertexCollection & vertices = *recVtxs.product();



	///  MUONS + TRACKS    MUONS + TRACKS    MUONS + TRACKS
	Handle <  vector <pat::Muon  >>thePATMuonHandle;
	iEvent.getByToken(muons_, thePATMuonHandle);
	Handle < vector < pat::GenericParticle > >thePATTrackHandle;
	iEvent.getByToken(tracks_, thePATTrackHandle);

	////  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0

	numPV       = vertices.size();
	numTrack    = thePATTrackHandle->size();
	numV0       = 0;

	if (numTrack< 3 ) return;


	int trig = 0;

	unsigned int NTRIGGERS = 14;

	///// 2017 , 2018. not 2016 !
	std::string TriggersToTest[NTRIGGERS] = {
	"HLT_Dimuon25_Jpsi","HLT_Dimuon20_Jpsi_Barrel_Seagulls", // 0, 1: inclusive dimuon jpsi
	"HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_JpsiTrkTrk_Displaced", // 2, 3: displaced jpistrk or jpsi trktrk
	"HLT_DoubleMu4_3_Jpsi_Displaced", "HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_Jpsi_Displaced",  // 4, 5, 6: prescaled, 4 for 2017, 5 for 2018
	"HLT_Dimuon18_PsiPrime", "HLT_Dimuon10_PsiPrime_Barrel_Seagulls", // 7, 8: inclusive dimuon psi2s
	"HLT_DoubleMu4_PsiPrimeTrk_Displaced", // 9: displaced psi2s trk
	"HLT_Dimuon0_Jpsi3p5_Muon2", // 10: triple-mu (jpsi + muon)
	"HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", // 11: jpsi + 2 trkmu (phi->mumu)
	"HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05" // 12, 13: jpsi+2 trk(phi->KK), 12 for 2017, 13 for 2018
	};

	std::string TriggerFilters[NTRIGGERS] = {
	"hltDisplacedmumuFilterDimuon25Jpsis","hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow", // 0, 1: inclusive dimuon jpsi
	"hltJpsiTkVertexFilter","hltJpsiTkTkVertexFilterPhiKstar", // 2, 3: displaced jpistrk or jpsi trktrk
	"hltDisplacedmumuFilterDoubleMu43Jpsi", "hltmumuFilterDoubleMu43Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi",  // 4, 5, 6: prescaled, 4 for 2017, 5 for 2018
	"hltDisplacedmumuFilterDimuon18PsiPrimes", "hltDisplacedmumuFilterDimuon10PsiPrimeBarrelnoCow", // 7, 8: inclusive dimuon psi2s
	"hltPsiPrimeTkVertexFilter", // 9: displaced psi2s trk
	"hltVertexmumuFilterJpsiMuon3p5", // 10: triple-mu (jpsi + muon)
	"hltMumuFilterDoubleMu2Jpsi", // 11: jpsi + 2 trkmu (phi->mumu)
	"hltJpsiTkTkVertexFilterPhiDoubleTrk1v2", "hltJpsiTkTkVertexFilterPhiDoubleTrk1v4" // 12, 13: jpsi+2 trk(phi->KK), 12 for 2017, 13 for 2018
	};

	if (run < 290000 )  ///  2016 !!! FOR MC DO MANUALLY
	{
	    /////  2016 !!!
	    std::string TriggersToTest2016[NTRIGGERS] = {
	    "HLT_Dimuon16_Jpsi","HLT_Dimuon10_Jpsi_Barrel", // 0, 1: inclusive dimuon jpsi, 2016 version
	    "HLT_DoubleMu4_JpsiTrk_Displaced","HLT_Dimuon20_Jpsi", // 2, 3: displaced jpsitrk, 3 REPLACED with 20_Jpsi !!!
	    "HLT_DoubleMu4_3_Jpsi_Displaced", "HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_Jpsi_Displaced",  // 4, 5, 6: prescaled, only 4 for 2016
	    "HLT_Dimuon13_PsiPrime", "HLT_Dimuon8_PsiPrime_Barrel", // 7, 8: inclusive dimuon psi2s, 2016 version
	    "HLT_DoubleMu4_PsiPrimeTrk_Displaced", // 9: displaced psi2s trk
	    "HLT_Dimuon0_Jpsi_Muon", // 10: triple-mu (jpsi + muon)
	    // "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", // 11: jpsi + 2 trkmu , NOT in 2016 !
	    "HLT_Dimuon6_Jpsi_NoVertexing", // 11: in 2016 jpsi no vertex very prescaled !
	    // "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05" // 12, 13: jpsi+2 trk(phi->KK), not in 2016 !
	    "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing", "HLT_Dimuon0er16_Jpsi_NoVertexing" // 12, 13: in 2016 jpsi no vertexing paths !
	    };
	    //
	    std::string TriggerFilters2016[NTRIGGERS] = {
	    "hltDisplacedmumuFilterDimuon16Jpsi","hltDisplacedmumuFilterDimuon10JpsiBarrel", // 0, 1: inclusive dimuon jpsi, 2016 version
	    "hltJpsiTkVertexFilter","hltDisplacedmumuFilterDimuon20Jpsi", // 2, 3: displaced jpsitrk, 3 REPLACED with 20_Jpsi !!!
	    "hltDisplacedmumuFilterDoubleMu43Jpsi", "hltmumuFilterDoubleMu43Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi",  // 4, 5, 6: prescaled, only 4 for 2016
	    "hltDisplacedmumuFilterDimuon13PsiPrime", "hltDisplacedmumuFilterDimuon8PsiPrimeBarrel", // 7, 8: inclusive dimuon psi2s, 2016 version
	    "hltPsiPrimeTkVertexFilter", // 9: displaced psi2s trk
	    "hltVertexmumuFilterJpsiMuon", // 10: triple-mu (jpsi + muon), 2016 version
	    // "hltDiMuonGlbOrTrk0zFiltered0p2v2", // 11: jpsi + 2 trkmu , not in 2016 !
	    "hltDimuon6JpsiL3Filtered", // 11: in 2016 jpsi no vertex very prescaled !
	    // "hltJpsiTkTkVertexFilterPhiDoubleTrk1v2", "hltJpsiTkTkVertexFilterPhiDoubleTrk1v4" // 12, 13: jpsi+2 trk(phi->KK), 12 for 2017, 13 for 2018
	    "hltDimuon0JpsiNoOSL3Filtered", "hltDimuon0JpsiOSL3Filtered" // 12, 13: in 2016 jpsi no vertexing paths !
	    };
	    //
	    for (unsigned int i = 0; i < NTRIGGERS; i++)
	    {
	        TriggersToTest[i] = TriggersToTest2016[i];
	        TriggerFilters[i] = TriggerFilters2016[i];
	    }
	}

	unsigned short int TriggersFired[NTRIGGERS] = {};

	std::string strfired;

	if (triggerResults_handle.isValid())
	{
	    const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
	    for (unsigned int i = 0; i < NTRIGGERS; i++)
	    {
	        for (int version = 1; version < 30; version++)
	        {
	            std::stringstream ss; // full trigger name
	            ss << TriggersToTest[i] << "_v" << version;
	            //
	            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
	            if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit))
	            {
	                trig += (1<<i);
	                TriggersFired[i] = 1;
	                strfired.append( ss.str());
	                strfired.append( " " );
	            }
	        }
	    }
	}
	else
	{

	}

	// ACCESS TRIGGER FILTERS INFO FOR TRK MATCHING in AOD {{{
	std::vector<TLorentzVector> p4_TRIG_TRKS_JpsiTrk;
	TLorentzVector p4TrgTrk;

	auto processName = triggerEvent->usedProcessName();
	// std::cout << "----> Input Filters" << std::endl;
	for (size_t i = 0; i < TriggerFilters_.size(); i++)
	{
	    // std::cout << "> " << TriggerFilters_[i] << std::endl;
	    const unsigned int filterIndex(triggerEvent->filterIndex(edm::InputTag(TriggerFilters_[i], "", processName)));
	    //
	    if (filterIndex < triggerEvent->sizeFilters())
	    { // this tells us if the filter has been fired
	        const trigger::Keys& keys(triggerEvent->filterKeys(filterIndex));
	        // const trigger::Vids& ids(triggerEvent->filterIds(filterIndex));
	        const size_t nObjects(keys.size());

	        for (size_t i = 0; i < nObjects; i++)
	        {
	            const trigger::TriggerObject& tObject(triggerObjects[keys[i]]);
	            // std::cout << "   " << i << " " << ids[i] << "/" << keys[i] << ": " << tObject.id() << " " << tObject.pt() << " " << tObject.eta() << " " << tObject.phi() << " " << tObject.mass() << endl;
	            if (abs(tObject.id())==321) //// trigger code assigns PID, in our case 321 = kaon
	            {
	                p4TrgTrk.SetPtEtaPhiM(tObject.pt(), tObject.eta(), tObject.phi(), tObject.mass());
	                p4_TRIG_TRKS_JpsiTrk.push_back(p4TrgTrk);
	            }
	        }
	    }
	}

	nCand = 0;

	float PM_sigma = 1.e-7;
	KinematicParticleFactoryFromTransientTrack pFactory;

	// get J_psi -> mu mu: first mu+
	for (vector <pat::Muon >::const_iterator iMuonP = thePATMuonHandle->begin(); iMuonP != thePATMuonHandle->end(); ++iMuonP)
	{
	    const reco::Muon * rmu1 = dynamic_cast < const reco::Muon * >(iMuonP->originalObject()); // first muo
	    if (iMuonP->pt()<3.0) continue;
	    if (! iMuonP->isSoftMuon(thePrimaryV)) continue;

	    // next search for mu-
	    for ( vector <pat::Muon >::const_iterator iMuonM = thePATMuonHandle->begin(); iMuonM != thePATMuonHandle->end(); ++iMuonM)
	    {
	        if (iMuonM->charge() * iMuonP->charge() > -0.5) continue;
	        const reco::Muon * rmu2 = dynamic_cast < const reco::Muon * >(iMuonM->originalObject());
	        if (muon::overlap(*rmu1, *rmu2)) continue;

	        TrackRef muTrackM = iMuonM->track();
	        if (muTrackM.isNull()) continue;
	        if (! iMuonM->isSoftMuon(thePrimaryV)) continue;

	        TrackRef glbTrackM, glbTrackP;
	        glbTrackP = iMuonP->track();
	        glbTrackM = iMuonM->track();
	        if( glbTrackP.isNull() || glbTrackM.isNull() )          continue;
	        //
	        TransientTrack muonPTT(glbTrackP, &(*bFieldHandle));
	        TransientTrack muonMTT(glbTrackM, &(*bFieldHandle));
	        //
	        if(!muonPTT.isValid()) continue;
	        if(!muonMTT.isValid()) continue;
	        if(iMuonM->pt()<3.0) continue;
	        if(fabs(iMuonP->eta())>2.4 || fabs(iMuonM->eta())>2.4)  continue;
	        if(!(glbTrackM->quality(reco::TrackBase::highPurity)))  continue;
	        if(!(glbTrackP->quality(reco::TrackBase::highPurity)))  continue;

	        // Get The J/Psi information from the tracks
	        TLorentzVector p4mup_0c,p4mum_0c;
	        p4mup_0c.SetPtEtaPhiM(iMuonP->pt(), iMuonP->eta(), iMuonP->phi(), PDG_MUON_MASS);
	        p4mum_0c.SetPtEtaPhiM(iMuonM->pt(), iMuonM->eta(), iMuonM->phi(), PDG_MUON_MASS);
	        //
	        float muon_sigma = PDG_MUON_MASS * 1.e-6;
	        float chi = 0.;
	        float ndf = 0.;
	        //
	        vector < RefCountedKinematicParticle > muonParticles;
	        muonParticles.push_back(pFactory.particle(muonPTT, PM_PDG_MUON_MASS, chi, ndf, muon_sigma));
	        muonParticles.push_back(pFactory.particle(muonMTT, PM_PDG_MUON_MASS, chi, ndf, muon_sigma));
	        KinematicParticleVertexFitter fitter;
	        RefCountedKinematicTree psiVFT_noC;
	        int fitgood = 1;
	        try
	        {
	            psiVFT_noC = fitter.fit(muonParticles);
	        }
	        catch (const VertexException &)
	        {
	            fitgood = 0;
	        }
	        if (fitgood == 0) continue;

	        if (!psiVFT_noC->isValid()) continue;
	        //
	        psiVFT_noC->movePointerToTheTop();
	        RefCountedKinematicParticle MUMUparticle    = psiVFT_noC->currentParticle();
	        RefCountedKinematicVertex   MUMUvtx         = psiVFT_noC->currentDecayVertex();
	        //
	        double MUMU_mass_c0 = MUMUparticle->currentState().mass();
	        if ( MUMU_mass_c0 < 2.90) continue;
	        if ( MUMU_mass_c0 > 3.95) continue;
	        //
	        int psi2 = 0;
	        if (MUMU_mass_c0 > 3.4)  {psi2=1;}
	        if (psi2==1)    continue;
	        ParticleMass    PM_MYPSI_MASS = PM_PDG_JPSI_MASS;
	        //
	        double MUMU_vtxprob = TMath::Prob(MUMUvtx->chiSquared(), MUMUvtx->degreesOfFreedom());
	        if(MUMU_vtxprob < 0.01) continue;
	        //
	        psiVFT_noC->movePointerToTheFirstChild();
	        RefCountedKinematicParticle MUP_cMUMU       = psiVFT_noC->currentParticle();
	        psiVFT_noC->movePointerToTheNextChild();
	        RefCountedKinematicParticle MUM_cMUMU       = psiVFT_noC->currentParticle();
	        //
	        ///// {{{ MUON TRIGGER MATCHING  ( WORKING !!!! )
	        unsigned short int TriggersMathed[NTRIGGERS] = {};

	        for (unsigned int i = 0; i < NTRIGGERS; i++)
	        {
	            const pat::TriggerObjectStandAloneCollection muHLTMatches1 = iMuonP->triggerObjectMatchesByFilter(TriggerFilters[i]);
	            const pat::TriggerObjectStandAloneCollection muHLTMatches2 = iMuonM->triggerObjectMatchesByFilter(TriggerFilters[i]);
	            //
	            if ( muHLTMatches1.size() > 0 && muHLTMatches2.size() > 0 )
	            {
	                TriggersMathed[i] = 1;
	            }
	        }

	        for (vector < pat::GenericParticle >::const_iterator iTrack1 = thePATTrackHandle->begin(); iTrack1 != thePATTrackHandle->end(); ++iTrack1)
	        {
	            if(iTrack1->pt() < 1.2) continue;
	            pat::GenericParticle patTrack1 = *iTrack1;
	            if (!(patTrack1.track()->quality(reco::TrackBase::highPurity))) continue;


	            TLorentzVector p4k1;
	            p4k1.SetPtEtaPhiM(iTrack1->pt(),iTrack1->eta(),iTrack1->phi(), PDG_KAON_MASS);
	            if (p4k1.DeltaR(p4mum_0c) < 0.01) continue; // корень из разности углов этта и фи
	            if (p4k1.DeltaR(p4mup_0c) < 0.01) continue;
	            //
	            if(iTrack1->track().key() == rmu1->track().key() || iTrack1->track().key() == rmu2->track().key()) continue;
	            if( fabs(iTrack1->eta()) > 3.0 ) continue;

	            TransientTrack kaon1TT(iTrack1->track(), &(*bFieldHandle) );
	            if (!kaon1TT.isValid()) continue;
	            //
	            if ((p4mup_0c + p4k1 + p4mum_0c).M() < 5.00) continue;
	            if ((p4mup_0c + p4k1 + p4mum_0c).M() > 5.60) continue;
	            //
	            ///////////////////////////////  B VERTEX FIT /////////////////////////////////////////
	            //
	            std::vector<RefCountedKinematicParticle> W_candidate_init;

	            W_candidate_init.push_back(pFactory.particle(muonPTT, PM_PDG_MUON_MASS, chi,ndf, muon_sigma));
	            W_candidate_init.push_back(pFactory.particle(muonMTT, PM_PDG_MUON_MASS, chi,ndf, muon_sigma));
	            W_candidate_init.push_back(pFactory.particle(kaon1TT, PM_PDG_KAON_MASS, chi,ndf, PM_sigma));
	            RefCountedKinematicTree ebVFT, eertexFitTree;

	            std::vector<RefCountedKinematicParticle> W_candidate = W_candidate_init;
	            KinematicParticleVertexFitter eFitter;
	            fitgood = 1;
	            try
	            {
	                ebVFT = eFitter.fit(W_candidate);
	            }
	            catch (const VertexException &)
	            {
	                fitgood = 0;
	            }
	            if (fitgood == 0) continue;
	            if (!ebVFT->isValid()) continue;
	            ebVFT->movePointerToTheTop();

	            W_candidate = W_candidate_init;
	            MultiTrackKinematicConstraint *ConstraintJpsiMass_1 = new TwoTrackMassKinematicConstraint(PM_MYPSI_MASS);
	            KinematicConstrainedVertexFitter ekcvFitter;
	            eertexFitTree = ekcvFitter.fit(W_candidate, ConstraintJpsiMass_1);
	            if (!eertexFitTree->isValid()) continue;

	            eertexFitTree->movePointerToTheTop();
	            RefCountedKinematicParticle wCandMC         = eertexFitTree->currentParticle();
	            RefCountedKinematicVertex   wDecayVertexMC  = eertexFitTree->currentDecayVertex();
	            if (!wDecayVertexMC->vertexIsValid())  continue;
	            //
	            if(wDecayVertexMC->chiSquared() < 0) continue;
	            double W_Prob_tmp   = TMath::Prob(wDecayVertexMC->chiSquared(), (int) wDecayVertexMC->degreesOfFreedom());
	            if(W_Prob_tmp < 0.02) continue;

	            double W_mass_cjp_tmp = wCandMC->currentState().mass()  ;
	            if(W_mass_cjp_tmp < 5.10) continue;
	            if(W_mass_cjp_tmp > 5.50) continue;
	            //
	            // get children

	            eertexFitTree->movePointerToTheFirstChild();
	            RefCountedKinematicParticle mu1CandMCW    = eertexFitTree->currentParticle();
	            eertexFitTree->movePointerToTheNextChild();
	            RefCountedKinematicParticle mu2CandMCW    = eertexFitTree->currentParticle();
	            eertexFitTree->movePointerToTheNextChild();
	            RefCountedKinematicParticle K1CandMCW    = eertexFitTree->currentParticle();
	            eertexFitTree->movePointerToTheTop();
	            //
	            GlobalVector mu1CandMCW_p = mu1CandMCW->currentState().kinematicParameters().momentum();
	            GlobalVector mu2CandMCW_p = mu2CandMCW->currentState().kinematicParameters().momentum();
	            GlobalVector K1CandMCW_p = K1CandMCW->currentState().kinematicParameters().momentum();
	            //
	            TLorentzVector p4Bp;
	            p4Bp.SetXYZM(wCandMC->currentState().globalMomentum().x(),wCandMC->currentState().globalMomentum().y(),wCandMC->currentState().globalMomentum().z(), W_mass_cjp_tmp );

	            TLorentzVector p4fitW_mu1, p4fitW_mu2, p4fitW_K1;
	            p4fitW_mu1.SetXYZM(mu1CandMCW_p.x(), mu1CandMCW_p.y(), mu1CandMCW_p.z(), PDG_MUON_MASS);
	            p4fitW_mu2.SetXYZM(mu2CandMCW_p.x(), mu2CandMCW_p.y(), mu2CandMCW_p.z(), PDG_MUON_MASS);
	            p4fitW_K1.SetXYZM(K1CandMCW_p.x(), K1CandMCW_p.y(), K1CandMCW_p.z(), PDG_KAON_MASS);

	            ///
	            //K2 for B+K-
	            for (vector < pat::GenericParticle >::const_iterator iTrack2 = thePATTrackHandle->begin(); iTrack2 != thePATTrackHandle->end(); ++iTrack2)
	            {
	                if(iTrack2 == iTrack1) continue;
	                if(iTrack2->track().key() == rmu1->track().key() || iTrack2->track().key() == rmu2->track().key()) continue;
	                if(iTrack2->pt() < 1.0) continue;
	                pat::GenericParticle patTrack2 = *iTrack2;
	                if (!(patTrack2.track()->quality(reco::TrackBase::highPurity))) continue;

	                if (iTrack2->charge()*iTrack1->charge() > -0.5 ) continue;

	                TLorentzVector p4k2;
	                p4k2.SetPtEtaPhiM(iTrack2->pt(),iTrack2->eta(),iTrack2->phi(), PDG_KAON_MASS);

	                if( fabs(iTrack2->eta()) > 3.0 ) continue;

	                TransientTrack kaon2TT(iTrack2->track(), &(*bFieldHandle) );
	                if (!kaon2TT.isValid()) continue;
	                //cuts for speed

	                if ((p4k2 + p4Bp).M() - p4Bp.M() + PDG_B_MASS > 6.9) continue;

	                //
	                ////////Pi meson///////////////
	                for (vector < pat::GenericParticle >::const_iterator iTrack3 = thePATTrackHandle->begin(); iTrack3 != thePATTrackHandle->end(); ++iTrack3)
	                {
	                    if(iTrack3 == iTrack2) continue;
	                    if(iTrack3 == iTrack1) continue;
	                    if(iTrack3->track().key() == rmu1->track().key() || iTrack3->track().key() == rmu2->track().key()) continue;
	                    if(iTrack3->pt() < 0.4) continue;
	                    pat::GenericParticle patTrack3 = *iTrack3;
	                    if (!(patTrack3.track()->quality(reco::TrackBase::highPurity))) continue;

	                    if (iTrack3->charge()*iTrack2->charge() > -0.5) continue;

	                    TLorentzVector p4pi3;
	                    p4pi3.SetPtEtaPhiM(iTrack3->pt(),iTrack3->eta(),iTrack3->phi(), PDG_PION_MASS);

	                    if( fabs(iTrack3->eta()) > 3.0 ) continue;

	                    TransientTrack pion3TT(iTrack3->track(), &(*bFieldHandle) );
	                    if (!pion3TT.isValid()) continue;
	                    //
	                    //Cuts for speed

	                    if ((p4pi3 + p4k2 +p4Bp).M() - p4Bp.M() + PDG_B_MASS > 6.7) continue;
	                    ///////////B+SECOND FIT//////////////
	                    //
	                    std::vector<RefCountedKinematicParticle> B_candidate_init;

	                    B_candidate_init.push_back(pFactory.particle(muonPTT, PM_PDG_MUON_MASS, chi,ndf, muon_sigma));
	                    B_candidate_init.push_back(pFactory.particle(muonMTT, PM_PDG_MUON_MASS, chi,ndf, muon_sigma));
	                    B_candidate_init.push_back(pFactory.particle(kaon1TT, PM_PDG_KAON_MASS, chi,ndf, PM_sigma));
	                    RefCountedKinematicTree xbVFT, vertexFitTree;

	                    std::vector<RefCountedKinematicParticle> B_candidate = B_candidate_init;
	                    KinematicParticleVertexFitter pFitter;
	                    fitgood = 1;
	                    try
	                    {
	                      xbVFT = pFitter.fit(B_candidate);
	                    }
	                    catch (const VertexException &)
	                    {
	                      fitgood = 0;
	                    }
	                    if (fitgood == 0) continue;
	                    if (!xbVFT->isValid()) continue;
	                    xbVFT->movePointerToTheTop();
	                    double B_mass_c0_tmp = xbVFT->currentParticle()->currentState().mass();
	                    //J/Psi mass constraint
	                    B_candidate = B_candidate_init;
	                    MultiTrackKinematicConstraint *ConstraintJpsiMass_2 = new TwoTrackMassKinematicConstraint(PM_MYPSI_MASS);
	                    KinematicConstrainedVertexFitter kcvFitter;
	                    vertexFitTree = kcvFitter.fit(B_candidate, ConstraintJpsiMass_2);
	                    if (!vertexFitTree->isValid()) continue;

	                    vertexFitTree->movePointerToTheTop();
	                    RefCountedKinematicParticle bCandMC         = vertexFitTree->currentParticle();
	                    RefCountedKinematicVertex   bDecayVertexMC  = vertexFitTree->currentDecayVertex();
	                    if (!bDecayVertexMC->vertexIsValid())  continue;
	                    //
	                    if(bDecayVertexMC->chiSquared() < 0) continue;
	                    double B_Prob_tmp   = TMath::Prob(bDecayVertexMC->chiSquared(), (int) bDecayVertexMC->degreesOfFreedom());

	                    double B_mass_cjp_tmp = bCandMC->currentState().mass();
	                    //
	                    // get children from final B fit

	                    vertexFitTree->movePointerToTheFirstChild();
	                    RefCountedKinematicParticle mu1CandMC    = vertexFitTree->currentParticle();
	                    vertexFitTree->movePointerToTheNextChild();
	                    RefCountedKinematicParticle mu2CandMC    = vertexFitTree->currentParticle();
	                    vertexFitTree->movePointerToTheNextChild();
	                    RefCountedKinematicParticle K1CandMC    = vertexFitTree->currentParticle();
	                    vertexFitTree->movePointerToTheTop();
	                    //
	                    GlobalVector mu1CandMC_p = mu1CandMC->currentState().kinematicParameters().momentum();
	                    GlobalVector mu2CandMC_p = mu2CandMC->currentState().kinematicParameters().momentum();
	                    GlobalVector K1CandMC_p = K1CandMC->currentState().kinematicParameters().momentum();
	                    //
	                    TLorentzVector p4fit_mu1, p4fit_mu2, p4fit_K1;
	                    p4fit_mu1.SetXYZM(mu1CandMC_p.x(), mu1CandMC_p.y(), mu1CandMC_p.z(), PDG_MUON_MASS);
	                    p4fit_mu2.SetXYZM(mu2CandMC_p.x(), mu2CandMC_p.y(), mu2CandMC_p.z(), PDG_MUON_MASS);
	                    p4fit_K1.SetXYZM(K1CandMC_p.x(), K1CandMC_p.y(), K1CandMC_p.z(), PDG_KAON_MASS);

	                    ///////////B+K-//////////////////////
	                    std::vector<RefCountedKinematicParticle> BK_candidate_init;

	                    BK_candidate_init.push_back(bCandMC);
	                    BK_candidate_init.push_back(pFactory.particle(kaon2TT, PM_PDG_KAON_MASS, chi,ndf, PM_sigma));
	                    BK_candidate_init.push_back(pFactory.particle(pion3TT, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
	                    RefCountedKinematicTree dertexFitTree;

	                    std::vector<RefCountedKinematicParticle> BK_candidate = BK_candidate_init;
	                    KinematicParticleVertexFitter dFitter; //KinematicParticleVertexFitter
	                    fitgood = 1;
	                    try
	                    {
	                        dertexFitTree = dFitter.fit(BK_candidate);
	                    }
	                    catch (const VertexException &)
	                    {
	                        fitgood = 0;
	                    }
	                    if (fitgood == 0) continue;
	                    if (!dertexFitTree->isValid()) continue;
	                    dertexFitTree->movePointerToTheTop();

	                    RefCountedKinematicParticle bkCandMC         = dertexFitTree->currentParticle();
	                    RefCountedKinematicVertex   bkDecayVertexMC  = dertexFitTree->currentDecayVertex();
	                    if (!bkDecayVertexMC->vertexIsValid())  continue;
	                    //
	                    if(bkDecayVertexMC->chiSquared() < 0) continue;
	                    double BK_Prob_tmp   = TMath::Prob(bkDecayVertexMC->chiSquared(), (int) bkDecayVertexMC->degreesOfFreedom());
	                    if(BK_Prob_tmp < 0.05) continue;


	                    //
	                    // get children from final B fit

	                    dertexFitTree->movePointerToTheFirstChild();
	                    RefCountedKinematicParticle BCandMC    = dertexFitTree->currentParticle();
	                    dertexFitTree->movePointerToTheNextChild();
	                    RefCountedKinematicParticle K2CandMC    = dertexFitTree->currentParticle();
	                    dertexFitTree->movePointerToTheTop();
	                    RefCountedKinematicParticle PI3CandMC    = dertexFitTree->currentParticle();
	                    dertexFitTree->movePointerToTheTop();
	                    //
	                    GlobalVector BCandMC_p = BCandMC->currentState().kinematicParameters().momentum();
	                    GlobalVector K2CandMC_p = K2CandMC->currentState().kinematicParameters().momentum();
	                    GlobalVector PI3CandMC_p =PI3CandMC->currentState().kinematicParameters().momentum();
	                    //
	                    TLorentzVector p4fit_B, p4fit_K2, p4fit_PI3;
	                    p4fit_B.SetXYZM(BCandMC_p.x(), BCandMC_p.y(), BCandMC_p.z(), BCandMC->currentState().mass());
	                    p4fit_K2.SetXYZM(K2CandMC_p.x(), K2CandMC_p.y(), K2CandMC_p.z(), PDG_KAON_MASS);
	                    p4fit_PI3.SetXYZM(PI3CandMC_p.x(), PI3CandMC_p.y(), PI3CandMC_p.z(), PDG_PION_MASS);

	                    double BK_mass_tmp = bkCandMC->currentState().mass() - BCandMC->currentState().mass() + PDG_B_MASS ;

	                    if(BK_mass_tmp > 6.60) continue;
	                    if(BK_mass_tmp < 5.00) continue;
	                    ////
	                    reco::Vertex bestVtxBSIP;  /// using XB fit and XB vtx HERE !!
	                     // {{{ GET THE BEST PV BY CHOSING THE BEST POINTING ANGLE AND REMOVE B TRACKS FROM ITS FIT
	                      // ********************* todos los vertices primarios con constrain del Beam-Spot y escogemos el de mejor pointing angle ****************

	                    reco::Vertex vtxBSrf ;

	                    Double_t pVtxBSIPX_temp = -10000.0;
	                    Double_t pVtxBSIPY_temp = -10000.0;
	                    Double_t pVtxBSIPZ_temp = -10000.0;
	                    Double_t pVtxBSIPXE_temp = -10000.0;
	                    Double_t pVtxBSIPYE_temp = -10000.0;
	                    Double_t pVtxBSIPZE_temp = -10000.0;
	                    Double_t pVtxBSIPCL_temp = -10000.0;
	                    Double_t pVtxBSIPdN_temp = 0;
	                    Double_t lip = -100000.0;
	                    for(size_t i = 0; i < recVtxs->size(); ++i)
	                    {
	                        const Vertex &vtxBS = (*recVtxs)[i];
	                        vector<reco::TransientTrack> vertexTracks;
	                        for ( std::vector<TrackBaseRef >::const_iterator iTrack = vtxBS.tracks_begin(); iTrack != vtxBS.tracks_end(); ++iTrack)
	                        {
	                            TrackRef trackRef = iTrack->castTo<TrackRef>();
	                            // the  tracks in the BK cand are  patTrack1  glbTrackP glbTrackM
	                            if (  !(      (glbTrackM            ==trackRef)
	                                        ||(glbTrackP             ==trackRef)
	                                        ||(patTrack1.track()     ==trackRef)
	                                        ||(patTrack2.track()     ==trackRef)
	                                        ||(patTrack3.track()     ==trackRef) ) )
	                            {
	                                TransientTrack tt(trackRef, &(*bFieldHandle) );
	                                vertexTracks.push_back(tt);
	                            }
	                        }

	                        // if no tracks in primary or no reco track included in primary then don't do anything
	                        vtxBSrf = vtxBS;
	                        GlobalPoint PVRfP = GlobalPoint( vtxBS.x(), vtxBS.y(), vtxBS.z() );
	                        if ( vertexTracks.size()>0 && (vtxBS.tracksSize()!=vertexTracks.size()) )
	                        {
	                            AdaptiveVertexFitter theFitter;
	                            TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);
	                            if ( v.isValid() )
	                            {
	                                vtxBSrf = reco::Vertex(v);
	                            }
	                        }
	                        Double_t dx = (*bkDecayVertexMC).position().x() - vtxBSrf.x();
	                        Double_t dy = (*bkDecayVertexMC).position().y() - vtxBSrf.y();
	                        Double_t dz = (*bkDecayVertexMC).position().z() - vtxBSrf.z();
	                        Double_t cosAlphaXYb = ( bkCandMC->currentState().globalMomentum().x() * dx + bkCandMC->currentState().globalMomentum().y()*dy + bkCandMC->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* bkCandMC->currentState().globalMomentum().mag() );
	                        if(cosAlphaXYb>lip)
	                        {
	                            lip = cosAlphaXYb ;
	                            pVtxBSIPX_temp     = vtxBSrf.x();
	                            pVtxBSIPY_temp     = vtxBSrf.y();
	                            pVtxBSIPZ_temp     = vtxBSrf.z();
	                            pVtxBSIPXE_temp    = vtxBSrf.covariance(0, 0);
	                            pVtxBSIPYE_temp    = vtxBSrf.covariance(1, 1);
	                            pVtxBSIPZE_temp    = vtxBSrf.covariance(2, 2);
	                            pVtxBSIPCL_temp    = (TMath::Prob(vtxBSrf.chi2(), (int)vtxBSrf.ndof()) );
	                            pVtxBSIPdN_temp    = vtxBS.tracksSize() - vertexTracks.size();
	                            bestVtxBSIP = vtxBSrf;
	                        }
	                    }
	                     // }}}
	                     //otlet B (imenno B a ne Bc)
	                    Double_t SIG_B_pv = -100000.0;
	                    SIG_B_pv = sqrt (    (bDecayVertexMC->position().x() - pVtxBSIPX_temp)*(bDecayVertexMC->position().x() - pVtxBSIPX_temp) / (pVtxBSIPXE_temp*pVtxBSIPXE_temp + bDecayVertexMC->error().cxx())    +
	                                        (bDecayVertexMC->position().y() - pVtxBSIPY_temp)*(bDecayVertexMC->position().y() - pVtxBSIPY_temp) / (pVtxBSIPYE_temp*pVtxBSIPYE_temp + bDecayVertexMC->error().cyy())    );

	                    if ( SIG_B_pv < 3.0 ) continue;
	                    //otlet Bc
	                    Double_t SIG_Bc_pv = -100000.0;
	                    SIG_Bc_pv = sqrt (    (bkDecayVertexMC->position().x() - pVtxBSIPX_temp)*(bkDecayVertexMC->position().x() - pVtxBSIPX_temp) / (pVtxBSIPXE_temp*pVtxBSIPXE_temp + bkDecayVertexMC->error().cxx())    +
	                                        (bkDecayVertexMC->position().y() - pVtxBSIPY_temp)*(bkDecayVertexMC->position().y() - pVtxBSIPY_temp) / (pVtxBSIPYE_temp*pVtxBSIPYE_temp + bkDecayVertexMC->error().cyy())    );
	                    if ( SIG_Bc_pv < 2.0 ) continue;

	                    ///
	                    ///// FOR TRIGGER MATCHING OF THE TRACK
	                    //
	                    Double_t tmp_drtrig = 99;
	                    Double_t tmp_dpttrig = 99;
	                    for (TLorentzVector p4TrgTrk : p4_TRIG_TRKS_JpsiTrk)
	                    {
	                        if (p4TrgTrk.DeltaR(p4k1) < tmp_drtrig)
	                        {
	                            tmp_drtrig = p4TrgTrk.DeltaR(p4k1);
	                            tmp_dpttrig = p4TrgTrk.Pt() / p4k1.Pt() - 1.0;
	                        }
	                    }
	                    ///

	                    //////////////////   SAVE
	                    BK_mass          ->push_back(BK_mass_tmp);
	                    BK_px            ->push_back(bkCandMC->currentState().globalMomentum().x());
	                    BK_py            ->push_back(bkCandMC->currentState().globalMomentum().y());
	                    BK_pz            ->push_back(bkCandMC->currentState().globalMomentum().z());
	                    BK_DecayVtxX     ->push_back(bkDecayVertexMC->position().x());
	                    BK_DecayVtxY     ->push_back(bkDecayVertexMC->position().y());
	                    BK_DecayVtxZ     ->push_back(bkDecayVertexMC->position().z());
	                    BK_DecayVtxXE    ->push_back(bkDecayVertexMC->error().cxx());
	                    BK_DecayVtxYE    ->push_back(bkDecayVertexMC->error().cyy());
	                    BK_DecayVtxZE    ->push_back(bkDecayVertexMC->error().czz());
	                    BK_Prob          ->push_back(BK_Prob_tmp);
	                    //
	                    B_mass          ->push_back(B_mass_cjp_tmp);
	                    B_mass_c0       ->push_back(B_mass_c0_tmp);
	                    B_px            ->push_back(bCandMC->currentState().globalMomentum().x());
	                    B_py            ->push_back(bCandMC->currentState().globalMomentum().y());
	                    B_pz            ->push_back(bCandMC->currentState().globalMomentum().z());
	                    B_px_CV         ->push_back( BCandMC_p.x());
	                    B_py_CV         ->push_back( BCandMC_p.y());
	                    B_pz_CV         ->push_back( BCandMC_p.z());
	                    B_DecayVtxX     ->push_back(bDecayVertexMC->position().x());
	                    B_DecayVtxY     ->push_back(bDecayVertexMC->position().y());
	                    B_DecayVtxZ     ->push_back(bDecayVertexMC->position().z());
	                    B_DecayVtxXE    ->push_back(bDecayVertexMC->error().cxx());
	                    B_DecayVtxYE    ->push_back(bDecayVertexMC->error().cyy());
	                    B_DecayVtxZE    ->push_back(bDecayVertexMC->error().czz());
	                    B_Prob          ->push_back(B_Prob_tmp);

	                    B_J_Prob        ->push_back(MUMU_vtxprob);
	                    B_J_mass        ->push_back( MUMU_mass_c0 );
	                    B_J_px          ->push_back( MUMUparticle->currentState().globalMomentum().x() );
	                    B_J_py          ->push_back( MUMUparticle->currentState().globalMomentum().y() );
	                    B_J_pz          ->push_back( MUMUparticle->currentState().globalMomentum().z() );
	                    B_J_DecayVtxX   ->push_back( MUMUvtx->position().x() );
	                    B_J_DecayVtxY   ->push_back( MUMUvtx->position().y() );
	                    B_J_DecayVtxZ   ->push_back( MUMUvtx->position().z() );
	                    B_J_DecayVtxXE  ->push_back( MUMUvtx->error().cxx() );
	                    B_J_DecayVtxYE  ->push_back( MUMUvtx->error().cyy() );
	                    B_J_DecayVtxZE  ->push_back( MUMUvtx->error().czz() );

	                    B_mu_px1_cjp    ->push_back(mu1CandMC_p.x());
	                    B_mu_py1_cjp    ->push_back(mu1CandMC_p.y());
	                    B_mu_pz1_cjp    ->push_back(mu1CandMC_p.z());
	                    B_mu_px1        ->push_back(p4mup_0c.Px());
	                    B_mu_py1        ->push_back(p4mup_0c.Py());
	                    B_mu_pz1        ->push_back(p4mup_0c.Pz());

	                    B_mu_px2_cjp    ->push_back(mu2CandMC_p.x());
	                    B_mu_py2_cjp    ->push_back(mu2CandMC_p.y());
	                    B_mu_pz2_cjp    ->push_back(mu2CandMC_p.z());
	                    B_mu_px2        ->push_back(p4mum_0c.Px());
	                    B_mu_py2        ->push_back(p4mum_0c.Py());
	                    B_mu_pz2        ->push_back(p4mum_0c.Pz());

	                    K1_px          ->push_back( p4k1.Px()                                             );
	                    K1_py          ->push_back( p4k1.Py()                                             );
	                    K1_pz          ->push_back( p4k1.Pz()                                             );
	                    K1_px_CV       ->push_back( K1CandMC_p.x());
	                    K1_py_CV       ->push_back( K1CandMC_p.y());
	                    K1_pz_CV       ->push_back( K1CandMC_p.z());
	                    K1_ch          ->push_back( iTrack1->charge()                                      );
	                    K1_vx          ->push_back( iTrack1->vx()                                          );
	                    K1_vy          ->push_back( iTrack1->vy()                                          );
	                    K1_vz          ->push_back( iTrack1->vz()                                          );
	                    K1_ips         ->push_back( fabs(patTrack1.track()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(patTrack1.track()->dxyError())) );
	                    KA1_drTRG       ->push_back( tmp_drtrig                                             );
	                    KA1_dptTRG      ->push_back( tmp_dpttrig                                            );

	                    K2_px          ->push_back( p4k2.Px()                                             );
	                    K2_py          ->push_back( p4k2.Py()                                             );
	                    K2_pz          ->push_back( p4k2.Pz()                                             );
	                    K2_px_CV       ->push_back( K2CandMC_p.x());
	                    K2_py_CV       ->push_back( K2CandMC_p.y());
	                    K2_pz_CV       ->push_back( K2CandMC_p.z());
	                    K2_ch          ->push_back( iTrack2->charge()                                      );
	                    K2_vx          ->push_back( iTrack2->vx()                                          );
	                    K2_vy          ->push_back( iTrack2->vy()                                          );
	                    K2_vz          ->push_back( iTrack2->vz()                                          );
	                    K2_ips         ->push_back( fabs(patTrack2.track()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(patTrack2.track()->dxyError())) );
	                    //
	                    PI3_px          ->push_back( p4pi3.Px()                                             );
	                    PI3_py          ->push_back( p4pi3.Py()                                             );
	                    PI3_pz          ->push_back( p4pi3.Pz()                                             );
	                    PI3_px_CV       ->push_back( PI3CandMC_p.x());
	                    PI3_py_CV       ->push_back( PI3CandMC_p.y());
	                    PI3_pz_CV       ->push_back( PI3CandMC_p.z());
	                    PI3_ch          ->push_back( iTrack3->charge()                                      );
	                    PI3_vx          ->push_back( iTrack3->vx()                                          );
	                    PI3_vy          ->push_back( iTrack3->vy()                                          );
	                    PI3_vz          ->push_back( iTrack3->vz()                                          );
	                    PI3_ips         ->push_back( fabs(patTrack3.track()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(patTrack3.track()->dxyError())) );
	                    //
	                    PV_becos_XX     ->push_back( pVtxBSIPX_temp                                         );
	                    PV_becos_YY     ->push_back( pVtxBSIPY_temp                                         );
	                    PV_becos_ZZ     ->push_back( pVtxBSIPZ_temp                                         );
	                    PV_becos_EX     ->push_back( pVtxBSIPXE_temp                                        );
	                    PV_becos_EY     ->push_back( pVtxBSIPYE_temp                                        );
	                    PV_becos_EZ     ->push_back( pVtxBSIPZE_temp                                        );
	                    PV_becos_CL     ->push_back( pVtxBSIPCL_temp                                        );
	                    PV_becos_dN     ->push_back( pVtxBSIPdN_temp                                        );

	                    trig0_fi        ->push_back( TriggersFired[0]                                       );
	                    trig1_fi        ->push_back( TriggersFired[1]                                       );
	                    trig2_fi        ->push_back( TriggersFired[2]                                       );
	                    trig3_fi        ->push_back( TriggersFired[3]                                       );
	                    trig4_fi        ->push_back( TriggersFired[4]                                       );
	                    trig5_fi        ->push_back( TriggersFired[5]                                       );
	                    trig6_fi        ->push_back( TriggersFired[6]                                       );
	                    trig7_fi        ->push_back( TriggersFired[7]                                       );
	                    trig8_fi        ->push_back( TriggersFired[8]                                       );
	                    trig9_fi        ->push_back( TriggersFired[9]                                       );
	                    trig10_fi       ->push_back( TriggersFired[10]                                      );
	                    trig11_fi       ->push_back( TriggersFired[11]                                      );
	                    trig12_fi       ->push_back( TriggersFired[12]                                      );
	                    trig13_fi       ->push_back( TriggersFired[13]                                      );

	                    trig0_ma        ->push_back( TriggersMathed[0]                                      );
	                    trig1_ma        ->push_back( TriggersMathed[1]                                      );
	                    trig2_ma        ->push_back( TriggersMathed[2]                                      );
	                    trig3_ma        ->push_back( TriggersMathed[3]                                      );
	                    trig4_ma        ->push_back( TriggersMathed[4]                                      );
	                    trig5_ma        ->push_back( TriggersMathed[5]                                      );
	                    trig6_ma        ->push_back( TriggersMathed[6]                                      );
	                    trig7_ma        ->push_back( TriggersMathed[7]                                      );
	                    trig8_ma        ->push_back( TriggersMathed[8]                                      );
	                    trig9_ma        ->push_back( TriggersMathed[9]                                      );
	                    trig10_ma       ->push_back( TriggersMathed[10]                                     );
	                    trig11_ma       ->push_back( TriggersMathed[11]                                     );
	                    trig12_ma       ->push_back( TriggersMathed[12]                                     );
	                    trig13_ma       ->push_back( TriggersMathed[13]                                     );

	                    nCand ++;
	                    BK_candidate_init.clear();
	                    BK_candidate.clear();
	                    B_candidate_init.clear();
	                    B_candidate.clear();
	                }//PI3
	            } //K2
	            W_candidate_init.clear();
	            W_candidate.clear();
	        }// tr1
	        muonParticles.clear();
	    } // muonN
	} // muonP

	// ===================== END OF EVENT : WRITE ETC ++++++++++++++++++++++

	if (nCand > 0)
	{
		wwtree->Fill();
	} // nCand > 0

	//
	BK_mass->clear();
	BK_px->clear();             BK_py->clear();              BK_pz->clear();
	BK_DecayVtxX->clear();      BK_DecayVtxY->clear();       BK_DecayVtxZ->clear();
	BK_DecayVtxXE->clear();     BK_DecayVtxYE->clear();      BK_DecayVtxZE->clear();
	BK_Prob->clear();
	//
	B_mass->clear();            B_mass_c0->clear();
	B_px->clear();              B_py->clear();              B_pz->clear();
	B_px_CV->clear();   B_py_CV->clear();   B_pz_CV->clear();
	B_DecayVtxX->clear();       B_DecayVtxY->clear();       B_DecayVtxZ->clear();
	B_DecayVtxXE->clear();      B_DecayVtxYE->clear();      B_DecayVtxZE->clear();
	B_Prob->clear();            B_J_Prob->clear();

	B_J_mass->clear();
	B_J_px->clear();            B_J_py->clear();            B_J_pz->clear();
	B_J_DecayVtxX->clear();     B_J_DecayVtxY->clear();     B_J_DecayVtxZ->clear();
	B_J_DecayVtxXE->clear();    B_J_DecayVtxYE->clear();    B_J_DecayVtxZE->clear();
	B_J_Prob->clear();

	B_mu_px1_cjp->clear();      B_mu_py1_cjp->clear();      B_mu_pz1_cjp->clear();
	B_mu_px1->clear();          B_mu_py1->clear();          B_mu_pz1->clear();

	B_mu_px2_cjp->clear();      B_mu_py2_cjp->clear();      B_mu_pz2_cjp->clear();
	B_mu_px2->clear();          B_mu_py2->clear();          B_mu_pz2->clear();


	K1_px->clear();      K1_py->clear();      K1_pz->clear();
	K1_px_CV->clear();   K1_py_CV->clear();   K1_pz_CV->clear();
	K1_ch->clear();
	K1_vx->clear();      K1_vy->clear();      K1_vz->clear();
	K1_ips->clear();
	KA1_drTRG->clear();   KA1_dptTRG->clear();

	K2_px->clear();      K2_py->clear();      K2_pz->clear();
	K2_px_CV->clear();   K2_py_CV->clear();   K2_pz_CV->clear();
	K2_ch->clear();
	K2_vx->clear();      K2_vy->clear();      K2_vz->clear();
	K2_ips->clear();

	PI3_px->clear();      PI3_py->clear();      PI3_pz->clear();
	PI3_px_CV->clear();   PI3_py_CV->clear();   PI3_pz_CV->clear();
	PI3_ch->clear();
	PI3_vx->clear();      PI3_vy->clear();      PI3_vz->clear();
	PI3_ips->clear();
	//
	PV_becos_XX->clear();   PV_becos_YY->clear();   PV_becos_ZZ->clear();
	PV_becos_EX->clear();   PV_becos_EY->clear();   PV_becos_EZ->clear();
	PV_becos_CL->clear();   PV_becos_dN->clear();

	trig0_fi->clear();      trig1_fi->clear();      trig2_fi->clear();
	trig3_fi->clear();      trig4_fi->clear();      trig5_fi->clear();
	trig6_fi->clear();      trig7_fi->clear();      trig8_fi->clear();
	trig9_fi->clear();      trig10_fi->clear();     trig11_fi->clear();
	trig12_fi->clear();     trig13_fi->clear();

	trig0_ma->clear();      trig1_ma->clear();      trig2_ma->clear();
	trig3_ma->clear();      trig4_ma->clear();      trig5_ma->clear();
	trig6_ma->clear();      trig7_ma->clear();      trig8_ma->clear();
	trig9_ma->clear();      trig10_ma->clear();     trig11_ma->clear();
	trig12_ma->clear();     trig13_ma->clear();

}


// ------------ method called once each job just before starting event loop  ------------
void Xb_frame::beginJob()
{
	using namespace std;
	using namespace reco;
	//
	cout << "------------------------------->>>>> Begin Job" << endl;

	f = new TFile(fileName.c_str(), "RECREATE");
	wwtree  = new TTree("wztree", "muons tree");

	wwtree->Branch("nCand"              , &nCand            , "nCand/I"     );

	wwtree->Branch("run"                , &run              , "run/I"       );
	wwtree->Branch("event"              , &event            , "event/I"     );
	wwtree->Branch("lumi"               , &lumi             , "lumi/F"      );

	wwtree->Branch("numPV"              , &numPV            , "numPV/I"     );
	wwtree->Branch("numTrack"           , &numTrack         , "numTrack/I"  );
	wwtree->Branch("numV0"              , &numV0            , "numV0/I"     );

	wwtree->Branch("BK_mass"             , &BK_mass           );
	wwtree->Branch("BK_px"               , &BK_px             );
	wwtree->Branch("BK_py"               , &BK_py             );
	wwtree->Branch("BK_pz"               , &BK_pz             );
	wwtree->Branch("BK_DecayVtxX"        , &BK_DecayVtxX      );
	wwtree->Branch("BK_DecayVtxY"        , &BK_DecayVtxY      );
	wwtree->Branch("BK_DecayVtxZ"        , &BK_DecayVtxZ      );
	wwtree->Branch("BK_DecayVtxXE"       , &BK_DecayVtxXE     );
	wwtree->Branch("BK_DecayVtxYE"       , &BK_DecayVtxYE     );
	wwtree->Branch("BK_DecayVtxZE"       , &BK_DecayVtxZE     );
	wwtree->Branch("BK_Prob"             , &BK_Prob           );

	wwtree->Branch("B_mass"             , &B_mass           );
	wwtree->Branch("B_mass_c0"          , &B_mass_c0        );
	wwtree->Branch("B_px"               , &B_px             );
	wwtree->Branch("B_py"               , &B_py             );
	wwtree->Branch("B_pz"               , &B_pz             );
	wwtree->Branch("B_px_CV"            , &B_px_CV        );
	wwtree->Branch("B_py_CV"            , &B_py_CV        );
	wwtree->Branch("B_pz_CV"            , &B_pz_CV        );
	wwtree->Branch("B_DecayVtxX"        , &B_DecayVtxX      );
	wwtree->Branch("B_DecayVtxY"        , &B_DecayVtxY      );
	wwtree->Branch("B_DecayVtxZ"        , &B_DecayVtxZ      );
	wwtree->Branch("B_DecayVtxXE"       , &B_DecayVtxXE     );
	wwtree->Branch("B_DecayVtxYE"       , &B_DecayVtxYE     );
	wwtree->Branch("B_DecayVtxZE"       , &B_DecayVtxZE     );
	wwtree->Branch("B_Prob"             , &B_Prob           );

	wwtree->Branch("B_J_Prob"           , &B_J_Prob         );
	wwtree->Branch("B_J_mass"           , &B_J_mass         );
	wwtree->Branch("B_J_px"             , &B_J_px           );
	wwtree->Branch("B_J_py"             , &B_J_py           );
	wwtree->Branch("B_J_pz"             , &B_J_pz           );
	wwtree->Branch("B_J_DecayVtxX"      , &B_J_DecayVtxX    );
	wwtree->Branch("B_J_DecayVtxY"      , &B_J_DecayVtxY    );
	wwtree->Branch("B_J_DecayVtxZ"      , &B_J_DecayVtxZ    );
	wwtree->Branch("B_J_DecayVtxXE"     , &B_J_DecayVtxXE   );
	wwtree->Branch("B_J_DecayVtxYE"     , &B_J_DecayVtxYE   );
	wwtree->Branch("B_J_DecayVtxZE"     , &B_J_DecayVtxZE   );

	wwtree->Branch("B_mu_px1_cjp"       , &B_mu_px1_cjp     );
	wwtree->Branch("B_mu_py1_cjp"       , &B_mu_py1_cjp     );
	wwtree->Branch("B_mu_pz1_cjp"       , &B_mu_pz1_cjp     );
	wwtree->Branch("B_mu_px1"           , &B_mu_px1         );
	wwtree->Branch("B_mu_py1"           , &B_mu_py1         );
	wwtree->Branch("B_mu_pz1"           , &B_mu_pz1         );
	wwtree->Branch("B_mu_px2_cjp"       , &B_mu_px2_cjp     );
	wwtree->Branch("B_mu_py2_cjp"       , &B_mu_py2_cjp     );
	wwtree->Branch("B_mu_pz2_cjp"       , &B_mu_pz2_cjp     );
	wwtree->Branch("B_mu_px2"           , &B_mu_px2         );
	wwtree->Branch("B_mu_py2"           , &B_mu_py2         );
	wwtree->Branch("B_mu_pz2"           , &B_mu_pz2         );

	//
	wwtree->Branch("K1_px"             , &K1_px           );
	wwtree->Branch("K1_py"             , &K1_py           );
	wwtree->Branch("K1_pz"             , &K1_pz           );
	wwtree->Branch("K1_px_CV"          , &K1_px_CV        );
	wwtree->Branch("K1_py_CV"          , &K1_py_CV        );
	wwtree->Branch("K1_pz_CV"          , &K1_pz_CV        );
	wwtree->Branch("K1_ch"             , &K1_ch           );
	wwtree->Branch("K1_vx"             , &K1_vx           );
	wwtree->Branch("K1_vy"             , &K1_vy           );
	wwtree->Branch("K1_vz"             , &K1_vz           );
	wwtree->Branch("K1_ips"            , &K1_ips          );
	wwtree->Branch("KA1_drTRG"          , &KA1_drTRG        );
	wwtree->Branch("KA1_dptTRG"         , &KA1_dptTRG       );
	//
	wwtree->Branch("K2_px"             , &K2_px           );
	wwtree->Branch("K2_py"             , &K2_py           );
	wwtree->Branch("K2_pz"             , &K2_pz           );
	wwtree->Branch("K2_px_CV"          , &K2_px_CV        );
	wwtree->Branch("K2_py_CV"          , &K2_py_CV        );
	wwtree->Branch("K2_pz_CV"          , &K2_pz_CV        );
	wwtree->Branch("K2_ch"             , &K2_ch           );
	wwtree->Branch("K2_vx"             , &K2_vx           );
	wwtree->Branch("K2_vy"             , &K2_vy           );
	wwtree->Branch("K2_vz"             , &K2_vz           );
	wwtree->Branch("K2_ips"            , &K2_ips          );
	//
	wwtree->Branch("PI3_px"             , &PI3_px           );
	wwtree->Branch("PI3_py"             , &PI3_py           );
	wwtree->Branch("PI3_pz"             , &PI3_pz           );
	wwtree->Branch("PI3_px_CV"          , &PI3_px_CV        );
	wwtree->Branch("PI3_py_CV"          , &PI3_py_CV        );
	wwtree->Branch("PI3_pz_CV"          , &PI3_pz_CV        );
	wwtree->Branch("PI3_ch"             , &PI3_ch           );
	wwtree->Branch("PI3_vx"             , &PI3_vx           );
	wwtree->Branch("PI3_vy"             , &PI3_vy           );
	wwtree->Branch("PI3_vz"             , &PI3_vz           );
	wwtree->Branch("PI3_ips"            , &PI3_ips          );
	//
	wwtree->Branch("PV_becos_XX"        , &PV_becos_XX      );
	wwtree->Branch("PV_becos_YY"        , &PV_becos_YY      );
	wwtree->Branch("PV_becos_ZZ"        , &PV_becos_ZZ      );
	wwtree->Branch("PV_becos_EX"        , &PV_becos_EX      );
	wwtree->Branch("PV_becos_EY"        , &PV_becos_EY      );
	wwtree->Branch("PV_becos_EZ"        , &PV_becos_EZ      );
	wwtree->Branch("PV_becos_CL"        , &PV_becos_CL      );
	wwtree->Branch("PV_becos_dN"        , &PV_becos_dN      );

	wwtree->Branch("trig0_fi"           , &trig0_fi         );
	wwtree->Branch("trig1_fi"           , &trig1_fi         );
	wwtree->Branch("trig2_fi"           , &trig2_fi         );
	wwtree->Branch("trig3_fi"           , &trig3_fi         );
	wwtree->Branch("trig4_fi"           , &trig4_fi         );
	wwtree->Branch("trig5_fi"           , &trig5_fi         );
	wwtree->Branch("trig6_fi"           , &trig6_fi         );
	wwtree->Branch("trig7_fi"           , &trig7_fi         );
	wwtree->Branch("trig8_fi"           , &trig8_fi         );
	wwtree->Branch("trig9_fi"           , &trig9_fi         );
	wwtree->Branch("trig10_fi"          , &trig10_fi         );
	wwtree->Branch("trig11_fi"          , &trig11_fi         );
	wwtree->Branch("trig12_fi"          , &trig12_fi         );
	wwtree->Branch("trig13_fi"          , &trig13_fi         );

	wwtree->Branch("trig0_ma"           , &trig0_ma         );
	wwtree->Branch("trig1_ma"           , &trig1_ma         );
	wwtree->Branch("trig2_ma"           , &trig2_ma         );
	wwtree->Branch("trig3_ma"           , &trig3_ma         );
	wwtree->Branch("trig4_ma"           , &trig4_ma         );
	wwtree->Branch("trig5_ma"           , &trig5_ma         );
	wwtree->Branch("trig6_ma"           , &trig6_ma         );
	wwtree->Branch("trig7_ma"           , &trig7_ma         );
	wwtree->Branch("trig8_ma"           , &trig8_ma         );
	wwtree->Branch("trig9_ma"           , &trig9_ma         );
	wwtree->Branch("trig10_ma"          , &trig10_ma         );
	wwtree->Branch("trig11_ma"          , &trig11_ma         );
	wwtree->Branch("trig12_ma"          , &trig12_ma         );
	wwtree->Branch("trig13_ma"          , &trig13_ma         );

}

// ------------ method called once each job just after ending the event loop  ------------
void
Xb_frame::endJob()
{

	using namespace std;
	cout << "------------------------------->>>>> End Job" << endl;
	f->WriteTObject(wwtree);
	delete wwtree;
	f->Close();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Xb_frame::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Xb_frame);
