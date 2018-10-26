#include "Pythia8/Pythia.h"
#include "/expsoft/fastjet-3.3.0-cms/include/fastjet/PseudoJet.hh"
#include "/expsoft/fastjet-3.3.0-cms/include/fastjet/ClusterSequence.hh"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TH1.h"

using namespace std;
using namespace Pythia8;

Float_t DeltaR(TParticle* part, TParticle* candidate) {
	Float_t DeltaPhi;
	Float_t DeltaPhip = part->Phi()-candidate->Phi();
	Float_t DeltaPhim = abs(part->Phi()-candidate->Phi()-TMath::TwoPi());
	Float_t DeltaEta = part->Eta()-candidate->Eta();
	if (DeltaPhip > DeltaPhim) DeltaPhi = DeltaPhim;
	else DeltaPhi = DeltaPhip;
	return sqrt(DeltaPhi*DeltaPhi+DeltaEta*DeltaEta);
}

Float_t Isolation(Int_t l, TParticle* part, TClonesArray* particles, vector<fastjet::PseudoJet> constituents) {
	Float_t isol = 0;
	for (Long64_t k = 0; k < constituents.size(); ++k) {
		TParticle* candidate = reinterpret_cast< TParticle* >(particles->At(constituents[k].user_index()));
		Bool_t cut2 = candidate->Pt() > 0.1 && candidate->Eta() < 2.4 && candidate->Eta() > -2.4;
		if (cut2 && candidate->GetStatusCode() > 0 && k!=part->GetDaughter(0) && k!=part->GetDaughter(1)) {
			Float_t deltaR = DeltaR(part, candidate);
			TParticlePDG* candpdg = candidate->GetPDG(0);
			if (deltaR < 0.4 && candpdg) {
				Float_t pt = candidate->Pt();
				Float_t mass = candidate->GetMass();
				Float_t Et = sqrt(pt*pt+mass*mass);
				isol += Et;
				if (candpdg->Charge()!=0) isol += pt;
			}
		}
	}
	return isol;
}

void jet_finder(TH1F* h1, string &input, Int_t a, Int_t b) {
	// Import the tree with the monte carlo data
	TFile* f = TFile::Open(input.c_str());
	TTree* t = dynamic_cast< TTree* >(f->Get("Pythia"));
	TClonesArray* particles = new TClonesArray("Event");
	t->SetBranchAddress("Event", &particles);

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

		// Populate the output vector
		inclusiveJets = clustSeq.inclusive_jets(a);
		fastjet::Selector select_rapidity = fastjet::SelectorAbsEtaMax(2.5);
		fastjet::Selector select_pt_max = fastjet::SelectorPtMax(b);
		fastjet::Selector select_all = select_rapidity && select_pt_max;
		vector <fastjet::PseudoJet> selectedJets = select_all(inclusiveJets);
		for (Int_t i = 0; i < selectedJets.size(); ++i) {
			vector <fastjet::PseudoJet> constituents = selectedJets[i].constituents();
			if (constituents.size() > 3) h1->Fill(selectedJets[i].pt());
		}
	}
}


void generator_jet_pt_comparison2() {
	TStopwatch timer;

	TCanvas* c = new TCanvas("c","c",1);
	TH1F* h11 = new TH1F("h11","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 60, 40, 100);
	TH1F* h21 = new TH1F("h21","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 60, 40, 100);
	TH1F* h31 = new TH1F("h31","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 60, 40, 100);
	TH1F* h41 = new TH1F("h41","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 60, 40, 100);
	TH1F* h12 = new TH1F("h12","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 100, 100, 300);
	TH1F* h22 = new TH1F("h22","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 100, 100, 300);
	TH1F* h32 = new TH1F("h32","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 100, 100, 300);
	TH1F* h42 = new TH1F("h42","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 100, 100, 300);
	TH1F* h13 = new TH1F("h13","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 50, 300, 500);
	TH1F* h23 = new TH1F("h23","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 50, 300, 500);
	TH1F* h33 = new TH1F("h33","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 50, 300, 500);
	TH1F* h43 = new TH1F("h43","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]", 50, 300, 500);

	string input1 = "/data/cpapage/jet-gen2K-pt40.root";
	string input2 = "/data/cpapage/jet-gen2K-pt80.root";
	string input3 = "/data/cpapage/jet-gen2K-pt120.root";
	string input4 = "/data/cpapage/jet-gen2K-pt160.root";

	jet_finder(h11,input1, 40, 100);
	jet_finder(h21,input2, 40, 100); h21->SetLineColor(632);
	jet_finder(h31,input3, 40, 100); h31->SetLineColor(417);
	jet_finder(h41,input4, 40, 100); h41->SetLineColor(401);
	jet_finder(h12,input1, 100, 300);
	jet_finder(h22,input2, 100, 300); h22->SetLineColor(632);
	jet_finder(h32,input3, 100, 300); h32->SetLineColor(417);
	jet_finder(h42,input4, 100, 300); h42->SetLineColor(401);
	jet_finder(h13,input1, 300, 500);
	jet_finder(h23,input2, 300, 500); h23->SetLineColor(632);
	jet_finder(h33,input3, 300, 500); h33->SetLineColor(417);
	jet_finder(h43,input4, 300, 500); h43->SetLineColor(401);

	THStack* h1 = new THStack("h1","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]");
	THStack* h2 = new THStack("h2","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]");
	THStack* h3 = new THStack("h3","pTHatMin comparison for p_{T}^{Jet}; p_{T}^{Jet} [GeV/c]");

	h1->Add(h11); h1->Add(h21); h1->Add(h31); h1->Add(h41);
	h2->Add(h12); h2->Add(h22); h2->Add(h32); h2->Add(h42);
	h3->Add(h13); h3->Add(h23); h3->Add(h33); h3->Add(h43);

	auto l1 = new TLegend(0.7,0.7,0.9,0.9); l1->SetHeader("pTHatMin","c");
	l1->AddEntry(h11,"40 GeV","l"); l1->AddEntry(h21,"80 GeV","l"); l1->AddEntry(h31,"120 GeV","l"); l1->AddEntry(h41,"160 GeV","l");

	c->Divide(2,2);
	c->cd(1); h1->Draw("nostack"); l1->Draw();
	c->cd(2); h2->Draw("nostack"); l1->Draw();
	c->cd(3); h3->Draw("nostack"); l1->Draw();
	c->Print("generator_jet_pthatmin_comparison.pdf");

	timer.Stop();
	timer.Print();
}
