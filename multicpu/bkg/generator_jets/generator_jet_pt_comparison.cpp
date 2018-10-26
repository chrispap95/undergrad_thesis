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

Int_t* jet_finder(Int_t count_jets, Int_t count_pi0, string &input, Float_t a, Float_t b) {
	// Import the tree with the monte carlo data
	TFile* f = TFile::Open(input.c_str());
	TTree* t = dynamic_cast< TTree* >(f->Get("Pythia"));
	TClonesArray* particles = new TClonesArray("Event");
	t->SetBranchAddress("Event", &particles);

	count_jets = 0;
	count_pi0 = 0;

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
			if (constituents.size() > 3) ++count_jets;
			for (Int_t j = 0; j < constituents.size(); ++j) {
				TParticle* part = reinterpret_cast< TParticle* >(particles->At(constituents[j].user_index()));
				Int_t mother1_id = part->GetMother(0);
				Int_t mother2_id = part->GetMother(1);
				Bool_t control = false;
				if (mother1_id>0) {
					TParticle* mother = reinterpret_cast< TParticle* >(particles->At(mother1_id));
					if (mother->GetPdgCode() == 111 && mother->Pt()>5.0) {
						if (Isolation(mother1_id,mother,particles,constituents)/mother->Pt() < 0.5) {
							++count_pi0;
							control = true;
						}
					}
				}
				if (mother2_id>0 && !control) {
					TParticle* mother = reinterpret_cast< TParticle* >(particles->At(mother2_id));
					if (mother->GetPdgCode() == 111 && mother->Pt() > 10.0) {
						if ((Isolation(mother2_id,mother,particles,constituents)/mother->Pt()) < 0.5) {
							++count_pi0;
							control = true;
						}
					}
				}
				if (control) break;
			}
		}
	}
	cout << input << '\t' << '\t' << count_jets << '\t' << count_pi0 << '\t' << (double)count_pi0/(double)count_jets << endl;
	Int_t output[2];
	output[0] = count_jets;
	output[1] = count_pi0;
	return output;
}


void generator_jet_pt_comparison() {

	TStopwatch timer;

	TCanvas* c = new TCanvas("c","c",1);
	TH1F* h1 = new TH1F("h1","pTHatMin comparison for jets in p_{T} = [30,50] GeV region; pTHatMin [GeV/c]", 300, 0, 300);
	TH1F* h2 = new TH1F("h2","pTHatMin comparison for jets in p_{T} = [50,70] GeV region; pTHatMin [GeV/c]", 300, 0, 300);
	TH1F* h3 = new TH1F("h3","pTHatMin comparison for jets in p_{T} > 70 GeV region; pTHatMin [GeV/c]", 300, 0, 300);

	Int_t count_jets, count_pi0;
	Int_t *output;

	string input = "/data/cpapage/jet-gen-pt20.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(20,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(20,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(20,*output); 
	input = "/data/cpapage/jet-gen-pt30.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(30,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(30,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(30,*output); 
	input = "/data/cpapage/jet-gen-pt35.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(35,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(35,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(35,*output); 
	input = "/data/cpapage/jet-gen-pt40.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(40,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(40,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(40,*output); 
	input = "/data/cpapage/jet-gen-pt45.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(45,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(45,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(45,*output); 
	input = "/data/cpapage/jet-gen-pt50.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(50,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(50,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(50,*output); 
	input = "/data/cpapage/jet-gen-pt55.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(55,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(55,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(55,*output); 
	input = "/data/cpapage/jet-gen-pt60.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(60,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(60,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(60,*output); 
	input = "/data/cpapage/jet-gen-pt65.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(65,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(65,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(65,*output); 
	input = "/data/cpapage/jet-gen-pt70.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(70,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(70,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(70,*output); 
	input = "/data/cpapage/jet-gen-pt80.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(80,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(80,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(80,*output); 
	input = "/data/cpapage/jet-gen-pt100.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(100,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(100,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(100,*output); 
	input = "/data/cpapage/jet-gen-pt125.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(125,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(125,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(125,*output); 
	input = "/data/cpapage/jet-gen-pt150.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(150,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(150,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(150,*output); 
	input = "/data/cpapage/jet-gen-pt175.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(175,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(175,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(175,*output); 
	input = "/data/cpapage/jet-gen-pt200.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(200,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(200,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(200,*output); 
	input = "/data/cpapage/jet-gen-pt225.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(225,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(225,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(225,*output); 
	input = "/data/cpapage/jet-gen-pt250.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(250,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(250,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(250,*output); 
	input = "/data/cpapage/jet-gen-pt275.root";
	output = jet_finder(count_jets, count_pi0, input, 30, 50);
	h1->SetBinContent(275,*output); 
	output = jet_finder(count_jets, count_pi0, input, 50, 70);
	h2->SetBinContent(275,*output); 
	output = jet_finder(count_jets, count_pi0, input, 70, 1000);
	h3->SetBinContent(275,*output); 

	gStyle->SetOptStat(0);
	h1->Draw("e1");
	c->Print("generator_jet_pthatmin_comparison.pdf(");
	h2->Draw("e1");
	c->Print("generator_jet_pthatmin_comparison.pdf");
	h3->Draw("e1");
	c->Print("generator_jet_pthatmin_comparison.pdf)");

	timer.Stop();
	timer.Print();
}
 
