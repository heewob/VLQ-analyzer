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
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include  "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <algorithm>
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <string>
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "/uscms/home/hwchung/nobackup/testAnalysis/CMSSW_10_6_29/src/sortJets.h"
#include <TCanvas.h>
#include "TH2.h"
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

   int recluster_multiplicity;
   double recluster_energy[1000], recluster_mass[1000];
    
   double sortMass[2];    
 
   double posSJ_mass;
   double negSJ_mass;
   double posSJ_nAK8;
   double negSJ_nAK8;
   double pos_presort_mass;
   double neg_presort_mass;
   double misc_presort_mass;


   //create a vector for 4 vectors
   std::vector<TLorentzVector> particle4vectors;
   std::vector<TLorentzVector> COM4vectors;
   std::vector<fastjet::PseudoJet> PseudoJetvectors;
   std::vector<TLorentzVector> superjet_pos;
   std::vector<TLorentzVector> superjet_neg;  
   std::vector<TLorentzVector> leftoverjet;

  //create some plots of jet mass
  TH2D* pos_vs_neg_mass;
  TH2D* pos_vs_neg_nAK8;
  TH1D* pos_presort;
  TH1D* neg_presort;
  TH1D* misc_presort;


};

//convert Pseudojet to Candidate (to be precise, a concrete subclass of reco::Candidate i.e. reco::LeafCandidate)
reco::LeafCandidate convertPseudoJetToCandidate(const fastjet::PseudoJet& psJet) {
double px = psJet.px();
double py = psJet.py();
double pz = psJet.pz();
double energy = psJet.E();
math::XYZTLorentzVector vec;
vec.SetPxPyPzE(px, py, pz, energy);
reco::LeafCandidate candidate;
candidate.setP4(vec);
return candidate;
}

//convert Candidate to Pseudojet
fastjet::PseudoJet convertCandidateToPseudoJet(const reco::LeafCandidate& candJet) {
double px = candJet.px();
double py = candJet.py();
double pz = candJet.pz();
double energy = candJet.energy();
fastjet::PseudoJet pseudoJet(px, py, pz, energy);
return pseudoJet;
}

//convert candidate to TLorentzVector
TLorentzVector convertCandidateToTLorentz(const reco::LeafCandidate& candJet) {
double px = candJet.px();
double py = candJet.py();
double pz = candJet.pz();
double energy = candJet.energy();
TLorentzVector TLV(px, py, pz, energy);
return TLV;
}


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

   tree->Branch("recluster_multiplicity", &recluster_multiplicity, "recluster_multiplicity/I");
   tree->Branch("recluster_energy", recluster_energy, "recluster_energy[recluster_multiplicity]/D");
   tree->Branch("recluster_mass", recluster_mass, "recluster_mass[recluster_multiplicity]/D");
 
   tree->Branch("sortMass", sortMass, "sortMass[2]/D");

   tree->Branch("posSJ_mass", &posSJ_mass, "posSJ_mass/D");
   tree->Branch("negSJ_mass", &negSJ_mass, "negSJ_mass/D");
   tree->Branch("posSJ_nAK8", &posSJ_nAK8, "posSJ_nAK8/D");
   tree->Branch("negSJ_nAK8", &negSJ_nAK8, "negSJ_nAK8/D");
   tree->Branch("pos_presort_mass", &pos_presort_mass, "pos_presort_mass/D");
   tree->Branch("neg_presort_mass", &neg_presort_mass, "neg_presort_mass/D");
   tree->Branch("misc_presort_mass", &misc_presort_mass, "misc_presort_mass/D");


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

//register TFileService
edm::Service<TFileService> fileService;

//create histograms using TFileService
pos_vs_neg_mass = fileService->make<TH2D>("pos_vs_neg_mass", "pos SJ mass vs negative SJ mass", 4000, 0, 4000, 4000, 0, 4000);
pos_vs_neg_nAK8 = fileService->make<TH2D>("pos_vs_neg_nAK8", "pos SJ nAK8 vs negative SJ nAK8", 4000, 0, 4000, 4000, 0, 4000);
pos_presort = fileService->make<TH1D>("pos_presort", "pos SJ mass before sortmass", 4000, 0, 4000);
neg_presort =fileService->make<TH1D>("neg_presort", "neg SJ mass before sortmass", 4000, 0, 4000);
misc_presort = fileService->make<TH1D>("misc_presort", "misc SJ mass before sortmass", 4000, 0, 4000);


   for(auto iJet = fatJets->begin(); iJet != fatJets->end(); iJet++)         ////////Over AK8 Jets
   {
      pat::Jet Jet(*iJet);
      //std::cout <<"is iJet PF Jet: " <<iJet->isPFJet()<<",   is Jet PF Jet:  " << Jet.isPFJet() << std::endl;
      //if (Jet.isPFJet() == 0) {std::cout << Jet.eta() <<std::endl; } 
      if (!Jet.isPFJet()) continue;
      if( (Jet.pt() < 300.) || (abs(Jet.eta()) > 2.5) ) continue;   //make a cut to AK8 jets
      //std::cout << "passed AK8 Jet cut" << std::endl;

      if( (Jet.neutralHadronEnergyFraction() >= 0.90) || (Jet.neutralEmEnergyFraction() >= 0.90) || (Jet.muonEnergyFraction() >= 0.80) || (Jet.chargedHadronEnergyFraction() <= 0) || (Jet.chargedEmEnergyFraction() >= 0.80) || (Jet.chargedMultiplicity() <= 0)  || (Jet.numberOfDaughters() <= 1) ) continue;    //jet ID
      //std::cout << "passed Jet ID cut" << std::endl;

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
     double ESum = 0;
     
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
        ESum += (parti->energy())*puppiweight;

     
        //storing weighted PUPPI 4 vec in array
        PUPPI_px[nAK8particle] = weighted4vec.Px();
        PUPPI_py[nAK8particle] = weighted4vec.Py();
        PUPPI_pz[nAK8particle] = weighted4vec.Pz();
        PUPPI_E[nAK8particle] = weighted4vec.E();
   
        nAK8particle++;  
        //std::cout << "nAK8particle++" << nAK8particle << std::endl;
     }

     //calculate COM 4vector
     TLorentzVector COM; //velocity of COM
     COM.SetPxPyPzE(pxSum/ESum, pySum/ESum, pzSum/ESum, 1);
   
     for (unsigned int i = 0; i < numDaughters; ++i)
     {
        const reco::Candidate* parti = iJet->daughter(i);
        //TLorentzVector Boosted;
        //const pat::PackedCandidate* candParti = (pat::PackedCandidate*) parti; //type casting to use puppi weight method
        //double puppiweight = candParti->puppiWeight();
        //double newpx = (parti->px())*puppiweight-(parti->energy())*puppiweight*(COM.Px()); //px' = px-m*V_c
        //double newpy = (parti->py())*puppiweight-(parti->energy())*puppiweight*(COM.Py()); //px' = px-m*V_c
        //double newpz = (parti->pz())*puppiweight-(parti->energy())*puppiweight*(COM.Pz()); //px' = px-m*V_c
        //double newE = (parti->energy())*puppiweight;
       
        //Boosted.SetPxPyPzE(newpx, newpy, newpz, newE);
        
	TLorentzVector Tparti(parti->px(), parti->py(), parti->pz(), parti->energy());
	Tparti.Boost(-COM.X(), -COM.Y(), -COM.Z());
        COM4vectors.push_back(Tparti);

        //storing boosted COM 4 vec in array
        COM_px[nCOMparticle] = Tparti.Px();
        COM_py[nCOMparticle] = Tparti.Py();
        COM_pz[nCOMparticle] = Tparti.Pz();
        COM_E[nCOMparticle] = Tparti.E();
  
        nCOMparticle++;
     } 
   //end of jet loops in an event >> i moved a baket that was here so pur it back if it doesn't work
 
   //reclustering COM boosted jets  
   int nom = 0;
   for (const auto& vec : COM4vectors)
   {
      fastjet::PseudoJet pseudojet(vec.Px(), vec.Py(), vec.Pz(), vec.E());
      PseudoJetvectors.push_back(pseudojet);
      nom++;
      //std::cout<< "pseudojet pushback"<< nom <<std::endl;
   }
   double R0 = 0.8;
   fastjet::JetDefinition jet_def0(fastjet::cambridge_algorithm, R0);
   fastjet::ClusterSequence cs_jet0(PseudoJetvectors, jet_def0);  // first parameter = fastjet::PseudoJet vector w/ COM frame particles
   std::vector<fastjet::PseudoJet> jetsFJ_jet0 = fastjet::sorted_by_E(cs_jet0.inclusive_jets(10.)); // get the new jets clustered this way and sort them by energy, setting a 10 GeV minimum value for jets     
   recluster_multiplicity = jetsFJ_jet0.size();  //counting number of reclustered jets in the event 
   for (unsigned int i =0 ; i < jetsFJ_jet0.size(); i++)
   {
      if (jetsFJ_jet0[i].m()>50) 
      {
      recluster_mass[i] = jetsFJ_jet0[i].m();
      recluster_energy[i] = jetsFJ_jet0[i].E();
      }
   }
  
   //convert the pseudojet to candidate
   std::vector<reco::LeafCandidate> CandVec;
   for (const auto& PSjet : jetsFJ_jet0)
      {reco::LeafCandidate jet = convertPseudoJetToCandidate(PSjet);
      CandVec.push_back(jet);
      }

   //calculate the thrust axis in the boosted frame.
   math::XYZVector thrustAxis = Thrust(CandVec.begin(), CandVec.end()).axis();

   //calculate the jet angle with the thrust axis and sort them with negative and positive cosine
   for (const auto& jet : CandVec)
      {
      math::XYZVector jetMom(jet.px(), jet.py(), jet.pz());
      double dotProd = jetMom.Dot(thrustAxis);
      double jetMag = jetMom.R();
      double thrustMag = thrustAxis.R();
      double cos = dotProd / (jetMag*thrustMag);
      TLorentzVector TLVjet = convertCandidateToTLorentz(jet);  
      if (cos > 0.8)
         {
         superjet_pos.push_back(TLVjet);
         }
      else if (cos < -0.8)
         {
         superjet_neg.push_back(TLVjet);
         }
      else
         {
         leftoverjet.push_back(TLVjet);
         }
      }

   //calculating mass of pos/neg/misc jets before sortMass for check up
   TLorentzVector posSum_preSort(0, 0, 0, 0);
   TLorentzVector negSum_preSort(0, 0, 0, 0);
   TLorentzVector miscSum_preSort(0, 0, 0, 0);
   for (unsigned int i = 0; i < superjet_pos.size(); i++)
   {
   posSum_preSort += superjet_pos[i];
   }
   for (unsigned int i = 0; i < superjet_neg.size(); i++)
   {
   negSum_preSort += superjet_neg[i];
   }
   for (unsigned int i = 0; i < leftoverjet.size(); i++)
   {
   miscSum_preSort += leftoverjet[i];
   }

   //fill in hist of mass pre sortMass
   pos_presort->Fill(posSum_preSort.M());
   neg_presort->Fill(negSum_preSort.M());
   misc_presort->Fill(miscSum_preSort.M());
   //fill in the branch
   pos_presort_mass = posSum_preSort.M();
   neg_presort_mass = negSum_preSort.M();
   misc_presort_mass = miscSum_preSort.M();


   //sorting the leftover jets into SuperJet pos/neg with mass
   std::vector<TLorentzVector> posSuperJet;
   std::vector<TLorentzVector> negSuperJet;
   if(leftoverjet.size() > 0)
   {
   sortJets testing(superjet_pos, superjet_neg, leftoverjet);
   posSuperJet = testing.finalSuperJet1;
   negSuperJet = testing.finalSuperJet2;
   }
   else
   {
   posSuperJet = superjet_pos;
   negSuperJet = superjet_neg;
   }
 
   //sum of the mass of each superjets
   TLorentzVector posSum(0, 0, 0, 0);
   TLorentzVector negSum(0, 0, 0, 0);
   for (unsigned int i = 0; i < posSuperJet.size(); i++)
   {
   //std::cout<<"positive super jet: "<<posSuperJet[i].E()<<std::endl;
   posSum += posSuperJet[i];
   }
   for (unsigned int i = 0; i < negSuperJet.size(); i++)
   {
   //std::cout<<"negative super jet: "<<negSuperJet[i].E()<<std::endl;
   negSum += negSuperJet[i];
   }

   if(posSum.M() > 0)
   {
   sortMass[0] = posSum.M();
   }
   if(negSum.M() > 0)
   {
   sortMass[1] = negSum.M();
   }

   
   //make TH2D of mass of pos/neg superjet
   pos_vs_neg_mass->Fill(posSum.M(), negSum.M());
   //make TH2D of number of AK8 jets in each pos/neg superjet
   pos_vs_neg_nAK8->Fill(posSuperJet.size(), negSuperJet.size());
   //Fill branches
   posSJ_mass = posSum.M();
   negSJ_mass = negSum.M();
   posSJ_nAK8 = posSuperJet.size();
   negSJ_nAK8 = negSuperJet.size();

   //here you've reached the end of the event, so this will create a new entry in your tree (representing the event) and fills the branches you created above.
   tree->Fill();
   //empty the vectors and arrays
   particle4vectors.clear();
   COM4vectors.clear();
   PseudoJetvectors.clear();
   posSuperJet.clear();
   negSuperJet.clear();
   superjet_pos.clear();
   superjet_neg.clear();
   leftoverjet.clear();
   sortMass[0] = 0;
   sortMass[1] = 0;

   for (int i = 0; i < 100; ++i)
   {AK8_pt[i] = 0;
    AK8_eta[i] = 0;
    AK8_mass[i] = 0;
   }
   for (int i = 0; i < 1000; ++i)
   {PUPPI_px[i] = 0;
    PUPPI_py[i] = 0;
    PUPPI_pz[i] = 0;
    PUPPI_E[i] = 0;
    COM_px[i] = 0;
    COM_py[i] = 0;
    COM_pz[i] = 0;
    COM_E[i] = 0;
    recluster_energy[i] = 0;
    recluster_mass[i] = 0;
   }

}



//clean up
delete pos_vs_neg_mass;
delete pos_vs_neg_nAK8;
delete pos_presort;
delete neg_presort;
delete misc_presort;


}


DEFINE_FWK_MODULE(jetAnalyzer);

