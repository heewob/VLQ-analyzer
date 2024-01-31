// system include files
#include <fastjet/JetDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/PtEtaPhiMass.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "PhysicsTools/CandUtils/interface/Thrust.h"
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include  "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <algorithm>
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <string>
using namespace reco;


class jetAnalyzer : public edm::EDAnalyzer
{
public:
   explicit jetAnalyzer(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);

   //need to have one "token" for each collection (=classes of reconstructed particles/objects) you want to include
   edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_;
   edm::EDGetTokenT<std::vector<pat::Jet>> fatJetToken_;
   edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;

   TTree * tree;
  
   double event_HT = 0;
   int nAK8 = 0;
   int nAK8particle = 0;
   int nCOMparticle = 0;

   double AK8_pt[100], AK8_eta[100], AK8_mass[100];
   double PUPPI_px[1000], PUPPI_py[1000], PUPPI_pz[1000], PUPPI_E[1000];
   double COM_px[1000], COM_py[1000], COM_pz[1000], COM_E[1000]; 
  
   //create a vector for 4 vectors
   std::vector<TLorentzVector> particle4vectors;
   std::vector<TLorentzVector> COM4vectors;
   
};


jetAnalyzer::jetAnalyzer(const edm::ParameterSet& iConfig)
{

   // these tokens are needed to grab the colletions of physics objects you're interested in. Each requires an "Input Tag" that is specified in the cfg file. 
   genPartToken_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPartCollection"));
   fatJetToken_ =    consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("fatJetCollection"));
   jetToken_    = consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetCollection"));

   edm::Service<TFileService> fs;


   //make a tree and set all the branches that you'll want to save 
   // a tree creates an entry for each event that holds all the variables (called "branches") you're interested in. After you are done running, you can loop over the branches and get these values, 
   // or simply take a look at them in a TBrowser.
   //the way this works for arrays is that you have to tell the tree how many values you want to save per event (for example, AK8_pt saves nAK8 values each event, one for each AK8 jet, and nAK8 must also be a tree variable)
   tree = fs->make<TTree>("tree", "tree");
   tree->Branch("nAK8", &nAK8, "nAK8/I");
   tree->Branch("nAK8particle", &nAK8particle, "nAK8particle/I");
   tree->Branch("nCOMparticle", &nCOMparticle, "nCOMparticle/I");

   tree->Branch("event_HT", &event_HT, "event_HT/D"); //the "D" is because this is a double type.

   tree->Branch("AK8_pt", AK8_pt, "AK8_pt[nAK8]/D");
   tree->Branch("AK8_eta", AK8_eta, "AK8_eta[nAK8]/D");
   tree->Branch("AK8_mass", AK8_mass, "AK8_mass[nAK8]/D");
  
   tree->Branch("PUPPI_px", PUPPI_px, "PUPPI_px[nAK8particle]/D");
   tree->Branch("PUPPI_py", PUPPI_py, "PUPPI_py[nAK8particle]/D");
   tree->Branch("PUPPI_pz", PUPPI_pz, "PUPPI_pz[nAK8particle]/D");
   tree->Branch("PUPPI_E", PUPPI_E, "PUPPI_E[nAK8particle]/D");
   tree->Branch("COM_px", COM_px, "COM_px[nCOMparticle]/D");
   tree->Branch("COM_py", COM_py, "COM_py[nCOMparticle]/D");
   tree->Branch("COM_pz", COM_pz, "COM_pz[nCOMparticle]/D");
   tree->Branch("COM_E", COM_E, "COM_E[nCOMparticle]/D");

}


void jetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   // This is where the actual analysis is done - you take the token that was defined earlier and use this to obtain the physics object collections (basically vectors with your physics objects).
   // You can see more physics object collections that can be imported here - https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017
   // You can then loop over all the objects in these collections and see the information for each.
   // Note: collections of physics objects are generally sorted in decreasing order of pt.

   ///////////////////////Gen Particles//////////////////////////////////

   // These are the two lines that take the collection token and return you your collection vector ( in this case, a collection vector generated (=truth level) particles ).
   edm::Handle<std::vector<reco::GenParticle>> genParticles;
   iEvent.getByToken(genPartToken_, genParticles);

   /////////////////////////////AK8 Jets///////////////////////////////
   // loop over the reconstructed large (R = 0.8) jets from the event 
   edm::Handle<std::vector<pat::Jet> > fatJets;
   iEvent.getByToken(fatJetToken_, fatJets);
   TLorentzVector LeadingAK8Jets(0,0,0,0); // create a TLorentz vector (https://root.cern.ch/doc/master/classTLorentzVector.html)  to find the combined mass of the leading two AK8 jets
   nAK8 = 0;
   nAK8particle = 0;
   nCOMparticle = 0;

   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   {
      pat::Jet Jet(*iJet);
      std::cout <<"is iJet PF Jet: " <<iJet->isPFJet()<<",   is Jet PF Jet:  " << Jet.isPFJet() << std::endl;
      if (Jet.isPFJet() == 0) {std::cout << Jet.eta() <<std::endl; } 
      if (!Jet.isPFJet()) continue;
      if( (Jet.pt() > 300.) && (abs(Jet.eta())< 2.5) ) continue;   //make a cut to AK8 jets
      std::cout << "passed AK8 Jet cut" << std::endl;

      if( (Jet.neutralHadronEnergyFraction() < 0.90) && (Jet.neutralEmEnergyFraction() < 0.90) && (Jet.muonEnergyFraction() < 0.80) && (Jet.chargedHadronEnergyFraction() > 0) && (Jet.chargedEmEnergyFraction() < 0.80) && (Jet.chargedMultiplicity() > 0)  && (Jet.numberOfDaughters() > 1) ) continue;    //jet ID
      std::cout << "passed Jet ID cut" << std::endl;

      AK8_mass[nAK8] = Jet.mass(); // save some quantities of the large jets
      AK8_pt[nAK8] = Jet.pt();
      AK8_eta[nAK8] = Jet.eta();

      if(nAK8<2)
      {
         LeadingAK8Jets+=TLorentzVector(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      }

      nAK8++;

     //look at the particles in the jet and find COM
     double pxSum = 0;
     double pySum = 0;
     double pzSum = 0;
     double massSum = 0;
     
     //get the number of daughters
     unsigned int numDaughters = iJet->numberOfDaughters();

     for (unsigned int i = 0; i < numDaughters; ++i)
     {
        const reco::Candidate* parti = iJet->daughter(i);
        const pat::PackedCandidate* candParti = (pat::PackedCandidate*) parti; //type casting to use puppi weight method
        double puppiweight = candParti->puppiWeight();
        TLorentzVector weighted4vec((parti->px())*puppiweight, (parti->py())*puppiweight, (parti->pz())*puppiweight, (parti->energy())*puppiweight);
        particle4vectors.push_back(weighted4vec);
        pxSum += (parti->px())*puppiweight;
        pySum += (parti->py())*puppiweight;
        pzSum += (parti->pz())*puppiweight;
        massSum += (parti->energy())*puppiweight;

     
        //storing weighted PUPPI 4 vec in array
        PUPPI_px[nAK8particle] = weighted4vec.Px();
        PUPPI_py[nAK8particle] = weighted4vec.Py();
        PUPPI_pz[nAK8particle] = weighted4vec.Pz();
        PUPPI_E[nAK8particle] = weighted4vec.E();
   
        nAK8particle++;  
     }

     //calculate COM 4vector
     TLorentzVector COM; //velocity of COM
     COM.SetPxPyPzE(pxSum/massSum, pySum/massSum, pzSum/massSum, 1);
   
     for (unsigned int i = 0; i < numDaughters; ++i)
     {
        const reco::Candidate* parti = iJet->daughter(i);
        TLorentzVector Boosted;
        const pat::PackedCandidate* candParti = (pat::PackedCandidate*) parti; //type casting to use puppi weight method
        double puppiweight = candParti->puppiWeight();
        double newpx = (parti->px())*puppiweight-(parti->energy())*puppiweight*(COM.Px()); //px' = px-m*V_c
        double newpy = (parti->py())*puppiweight-(parti->energy())*puppiweight*(COM.Py()); //px' = px-m*V_c
        double newpz = (parti->pz())*puppiweight-(parti->energy())*puppiweight*(COM.Pz()); //px' = px-m*V_c
        double newE = (parti->energy())*puppiweight;
        Boosted.SetPxPyPzE(newpx, newpy, newpz, newE);
        COM4vectors.push_back(Boosted);

        //storing boosted COM 4 vec in array
        COM_px[nCOMparticle] = Boosted.Px();
        COM_py[nCOMparticle] = Boosted.Py();
        COM_pz[nCOMparticle] = Boosted.Pz();
        COM_E[nCOMparticle] = Boosted.E();
  
        nCOMparticle++;
     }

   //here you've reached the end of the event, so this will create a new entry in your tree (representing the event) and fills the branches you created above.
   tree->Fill();

  }
}
DEFINE_FWK_MODULE(jetAnalyzer);

