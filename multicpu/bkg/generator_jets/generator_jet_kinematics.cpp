#include "/expsoft/fastjet-3.3.0-cms/include/fastjet/PseudoJet.hh"
#include "/expsoft/fastjet-3.3.0-cms/include/fastjet/ClusterSequence.hh"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TH1.h"

using namespace std;

void generator_jet_kinematics() {

	TStopwatch timer;

	// Import the tree with the monte carlo data
	TFile* f = TFile::Open("/data/cpapage/jet-gen200K-pt120.root");
	TTree* t = dynamic_cast< TTree* >(f->Get("Pythia"));
	TClonesArray* particles = new TClonesArray("Event");
	t->SetBranchAddress("Event", &particles);

	//Define some histograms
	TCanvas* c = new TCanvas("c","analysis window", 1);
	TH1F* h1 = new TH1F("h1","p_{T} of Jets in #hat{p}_{T min} = 120 GeV sample; p_{T} [GeV]; Entries / 1.3 GeV",100,50,180);
	TH1F* h2 = new TH1F("h2","#eta of Jets in #hat{p}_{T min} = 120 GeV sample; #eta",100,-8,8);
  	TH1F* h3 = new TH1F("h3","#phi of Jets in #hat{p}_{T min} = 120 GeV sample; #phi",100,0,6.3);
	TH1F* h4 = new TH1F("h4","Energy of Jets in #hat{p}_{T min} = 120 GeV sample; E [GeV]; Entries / 6 GeV",100,0,600);

	// Set up the jet finder FastJet
	Double_t Rparam = 0.5;
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
	fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rparam, recombScheme, strategy);
 
	// This is the Event loop
	for (Int_t iEvent = 0; iEvent < t->GetEntries(); ++iEvent) {
		t->GetEntry(iEvent);
		vector<fastjet::PseudoJet> fjInputs;

		// This is the particle loop
		for (Int_t ip = 0; ip < particles->GetEntries(); ++ip) {
			TParticle* part = reinterpret_cast< TParticle* >(particles->At(ip));
			if (part->GetStatusCode() <= 0) continue; 

			// Feed FastJet input
			fjInputs.push_back(fastjet::PseudoJet( part->Px(), part->Py(), part->Pz(), part->Energy() ) );
			fjInputs[fjInputs.size()-1].set_user_index(ip);
		}

		// Check that FastJet input is not empty
		if (fjInputs.size() == 0) {
			cout << "Error: event with no final state particles" << endl;
		}

		// Initiate FastJet
		vector <fastjet::PseudoJet> inclusiveJets;
		fastjet::ClusterSequence clustSeq(fjInputs, jetDef);

		// Output some info for the first run of the algorithm
		if (!iEvent) {
			cout << "Ran " << jetDef.description() << endl;
			cout << "Strategy adopted by FastJet was "
			     << clustSeq.strategy_string() << endl << endl;
		}

		// Populate the output vector
		inclusiveJets = clustSeq.inclusive_jets(50.0);
		fastjet::Selector select_rapidity = fastjet::SelectorAbsEtaMax(100);
		// fastjet::Selector select_ptmax = fastjet::SelectorPtMax(70.0);
	        fastjet::Selector select_all = select_rapidity;// && select_ptmax;
		vector <fastjet::PseudoJet> selectedJets = select_all(inclusiveJets);
		for (Int_t i = 0; i < selectedJets.size(); ++i) {
			h1->Fill(selectedJets[i].pt());
			h2->Fill(selectedJets[i].eta());
			h3->Fill(selectedJets[i].phi());
			h4->Fill(selectedJets[i].e());
		}
	}

	c->Divide(2,2);
	c->cd(1); h1->Draw();
	c->cd(2); h2->Draw();
	c->cd(3); h3->Draw();
	c->cd(4); h4->Draw();
	c->Print("jets_kinematics_pt120.eps");

	TFile* fsave = new TFile("jets_kinematics_pt120.root","RECREATE");
	h1->Write();
	h2->Write();
	h3->Write();
	h4->Write();
	fsave->Close();
	timer.Stop();
	timer.Print();
}
