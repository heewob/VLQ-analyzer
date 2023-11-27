////////////////////////////HELP////////////////////////////////
//////////////Test EDAnalyzer to be filled in//////////////


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


// The EDAnalyzer framework is based off the edm::EDAnalyzer class
// You create your own, customized version of the class by inhereting from the edm::EDAnalyzer class and then overwriting constructor and analyzer methods. 
// The way that this works (at least as far as I know) is that when you run cmsRun on your cfg file, the settings you pass are used to create an instance of this class,
// and the analyze method is called for each event. 
// There are more edm::EDAnalyzer methods here that are not being used, such as beginRun, endRun, that some people use. 

class skeletonAnalyzer : public edm::EDAnalyzer 
{
public:
   explicit skeletonAnalyzer(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);



   //need to have one "token" for each collection (=classes of reconstructed particles/objects) you want to include
   edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_; 
   edm::EDGetTokenT<std::vector<pat::Jet>> fatJetToken_;
   edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;

   TTree * tree;

   //define all the =variables here that you will want to use or save later on here in the class object definition
   int nAK4 = 0;
   int nAK8 = 0;
   double event_HT = 0;

   double AK8_pt[100], AK8_eta[100], AK8_mass[100];
   double AK4_pt[100], AK4_eta[100], AK4_mass[100];

   int nTprime = 0;
   int Tprime_pdgid = 8000001;     // this is the numeric identifier for Tprime particles, called a PDG ID = Particle Data Group ID. You can see this for more particles here - https://pdg.lbl.gov/2023/mcdata/mc_particle_id_contents.html
};



skeletonAnalyzer::skeletonAnalyzer(const edm::ParameterSet& iConfig)
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

   tree->Branch("nAK4", &nAK4, "nAK4/I");  //the "/I" in the last part here says that this is an integer. You can get a lot of weird problems if you have confilcting types. Notice objects are passed by reference to the tree.
   tree->Branch("nAK8", &nAK8, "nAK8/I");
   tree->Branch("event_HT", &event_HT, "event_HT/D"); //the "D" is because this is a double type.

   tree->Branch("AK8_pt", AK8_pt, "AK8_pt[nAK8]/D");
   tree->Branch("AK8_eta", AK8_eta, "AK8_eta[nAK8]/D");

   tree->Branch("AK4_pt", AK4_pt, "AK4_pt[nAK8]/D");
   tree->Branch("AK4_eta", AK4_eta, "AK4_eta[nAK8]/D");
   
   tree->Branch("AK4_mass",AK4_mass, "AK4_mass[nAK8]/D");
   tree->Branch("AK8_mass",AK8_mass, "AK8_mass[nAK8]/D");

}



void skeletonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   // This is where the actual analysis is done - you take the token that was defined earlier and use this to obtain the physics object collections (basically vectors with your physics objects).
   // You can see more physics object collections that can be imported here - https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017
   // You can then loop over all the objects in these collections and see the information for each.
   // Note: collections of physics objects are generally sorted in decreasing order of pt.

   ///////////////////////Gen Particles//////////////////////////////////

   // These are the two lines that take the collection token and return you your collection vector ( in this case, a collection vector generated (=truth level) particles ).
   edm::Handle<std::vector<reco::GenParticle>> genParticles;
   iEvent.getByToken(genPartToken_, genParticles);

   nTprime = 0;
   //loop over the gen particles (the "true" simulated particles that were generated for this event. These are not the "reconstructed" objects that have detector and software reconstruction effects factored in)
   for (auto iG = genParticles->begin(); iG != genParticles->end(); iG++) 
   {
      if ((abs(iG->pdgId()) == Tprime_pdgid) && (iG->isLastCopy()))
      {
         std::cout << "Here is a TPrime, px/py/pz/E is " << iG->px() << "/" << iG->py() << "/" << iG->pz() << "/" << iG->energy() << std::endl;
         nTprime++;
      }
   
   }
   ///////////////////////////AK4 Jets////////////////////////////////
   // loop over the reconstructed small (R = 0.4) jets from the event 
   edm::Handle<std::vector<pat::Jet> > smallJets;
   iEvent.getByToken(jetToken_, smallJets);
   nAK4 = 0;
   event_HT = 0;
   for(auto iJet = smallJets->begin(); iJet != smallJets->end(); iJet++) 
   {
      if( (iJet->pt() < 30.)|| (abs(iJet->eta()) >2.5)) continue;   //make a cut to AK4 jets, usually jets are only reliable out to |eta| < 2.5 because of detector acceptance
      event_HT+=iJet->pt();  
      AK4_mass[nAK4] = iJet->mass();   // save some quantities of the small jets
      AK4_pt[nAK4]   = iJet->pt();
      AK4_eta[nAK4]  = iJet->eta();
      nAK4++;
   }
   
   /////////////////////////////AK8 Jets///////////////////////////////
   // loop over the reconstructed large (R = 0.8) jets from the event 
   edm::Handle<std::vector<pat::Jet> > fatJets;
   iEvent.getByToken(fatJetToken_, fatJets);
   TLorentzVector LeadingAK8Jets(0,0,0,0); // create a TLorentz vector (https://root.cern.ch/doc/master/classTLorentzVector.html)  to find the combined mass of the leading two AK8 jets
   nAK8 = 0;
   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   {
      if( (iJet->pt() < 100.) || (abs(iJet->eta())> 2.5) ) continue;   //make a cut to AK8 jets
      AK8_mass[nAK8] = iJet->mass(); // save some quantities of the large jets
      AK8_pt[nAK8] = iJet->pt();
      AK8_eta[nAK8] = iJet->eta();

      if(nAK8<2)
      {
         LeadingAK8Jets+=TLorentzVector(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      }

      nAK8++;
   }

   // Now how can you find the mass of the leading two AK8 jets and store this as a variable in your tree?
   // What other quantities might be useful to save in order to understand the kinematics of your event? 

   std::cout << "--------------End of event - the HT was " << event_HT << " and there were " << nTprime << " Tprimes. --------------" << std::endl;


   //here you've reached the end of the event, so this will create a new entry in your tree (representing the event) and fills the branches you created above.
   tree->Fill();

}   
DEFINE_FWK_MODULE(skeletonAnalyzer);
