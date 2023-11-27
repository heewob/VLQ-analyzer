// -*- C++ -*-
//
// Package:    leptonVeto/leptonVeto
// Class:      leptonVeto
// 
/**\class leptonVeto leptonVeto.cc leptonVeto/leptonVeto/plugins/leptonVeto.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ethan Cannaert
//         Created:  Mon, 29 Mar 2021 15:55:02 GMT
//
//


// system include files
// user include files
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include <fastjet/JetDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/Filter.hh>
#include <fastjet/ClusterSequence.hh>
//#include <fastjet/ActiveAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// new includes
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
//#include "PhysicsTools/CandUtils/interface/Thrust.h"
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
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <string>
// class declaration
//

class leptonVeto : public edm::stream::EDFilter<> {
   public:
      explicit leptonVeto(const edm::ParameterSet&);
      ~leptonVeto();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
      edm::EDGetTokenT<std::vector<pat::Electron>> electronToken_;
      edm::EDGetTokenT<std::vector<pat::Tau>> tauToken_;
	  std::string runType;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
leptonVeto::leptonVeto(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
    muonToken_ =    consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"));
    electronToken_ =    consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electronCollection"));
    tauToken_ =    consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("tauCollection"));

   //runType = iConfig.getParameter<std::string>("runType");

}


leptonVeto::~leptonVeto()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
leptonVeto::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //using namespace edm;
   /*#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
   #endif */
    int nE =0, nTau = 0, nMuon = 0;
    edm::Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonToken_, muons);
    for(auto iM = muons->begin(); iM != muons->end();iM++)
    {
      if( (iM->passed(reco::Muon::CutBasedIdMedium)) && (iM->pt() > 8.) && (abs(iM->eta()) < 2.4)  && (iM->passed(reco::Muon::PFIsoMedium))   ) nMuon++;
    }


    edm::Handle<std::vector<pat::Electron>> electrons;
    iEvent.getByToken(electronToken_, electrons);
    for(auto iE = electrons->begin(); iE != electrons->end();iE++)
    {
      if((iE->electronID("mvaEleID-Fall17-iso-V2-wp90")) && (iE->pt() > 12.) && (abs(iE->eta())<2.5 )    )nE++;   //medium WP
    }


    /* No tau vetoeing because it's too much work 

    
    edm::Handle<std::vector<pat::Tau>> taus;
    iEvent.getByToken(tauToken_, taus);
    for(auto iT = taus->begin(); iT != taus->end();iT++)
    {
      if((iT->pt() > 20.) && (abs(iT->eta()) < 2.3)  )
      {
        if (iT->tauID("byMediumDeepTau2017v2p1VSjet"))   //was medium before
        {
            if (iT->tauID("byLooseDeepTau2017v2p1VSe") && iT->tauID("byMediumDeepTau2017v2p1VSmu"))
            {
               nTau++;
            }
        }
      }
            ///byVTightDeepTau2017v2p1VSjet doesn't work on 2018 MC
    }
    */ 
    //std::cout << "E/Mu/Tau: " << nE << "/" << nMuon << "/" << nTau << std::endl;

    if( (nTau > 0) || (nE > 0) || (nMuon > 0) )return false;
    else {return true;}
   
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
leptonVeto::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
leptonVeto::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
leptonVeto::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
leptonVeto::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
leptonVeto::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
leptonVeto::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
leptonVeto::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(leptonVeto);
