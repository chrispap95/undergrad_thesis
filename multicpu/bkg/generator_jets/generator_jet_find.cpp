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

void generator_jet_find() {
	TStopwatch timer;

	// Import the tree with the monte carlo data
	TFile* f = TFile::Open("/data/cpapage/jet-gen200K-pt120.root");
	TTree* t = dynamic_cast< TTree* >(f->Get("Pythia"));
	TClonesArray* particles = new TClonesArray("Event");
	t->SetBranchAddress("Event", &particles);

	//Define some histograms
	TCanvas* c = new TCanvas("c","analysis window", 1);
	// TH1F* h1 = new TH1F("h1","p_{T} of Jets and pi0 jets; p_{T} [GeV]",100,70,200);
	// TH1F* h2 = new TH1F("h2","p_{T} of Jets and pi0 jets; p_{T} [GeV]",100,70,200);
	TH2F* h2d = new TH2F("h2d","p_{T}^{isolated #pi^0}/p_{T}^{jet} vs p_{T}^{jet}",200,30,200,200,0,1);

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
		vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
		fastjet::ClusterSequence clustSeq(fjInputs, jetDef);

		// Output some info for the first run of the algorithm
		if (!iEvent) {
			cout << "Ran " << jetDef.description() << endl;
			cout << "Strategy adopted by FastJet was "
			     << clustSeq.strategy_string() << endl << endl;
		}

		// Populate the output vector
		inclusiveJets = clustSeq.inclusive_jets(30.0);
		fastjet::Selector select_rapidity = fastjet::SelectorAbsEtaMax(2.37);
		fastjet::Selector select_ptmax = fastjet::SelectorPtMax(200.0);
		fastjet::Selector select_all = select_rapidity && select_ptmax;
		vector <fastjet::PseudoJet> selectedJets = select_all(inclusiveJets);
		for (Int_t i = 0; i < selectedJets.size(); ++i) {
			vector <fastjet::PseudoJet> constituents = selectedJets[i].constituents();
			//if (constituents.size() > 3) h1->Fill(selectedJets[i].pt());
			for (Int_t j = 0; j < constituents.size(); ++j) {
				TParticle* part = reinterpret_cast< TParticle* >(particles->At(constituents[j].user_index()));
				Int_t mother1_id = part->GetMother(0);
				Int_t mother2_id = part->GetMother(1);
				Bool_t control = false; //Variable that checks if a pi0 has been found during the check for the first mother 
				if (mother1_id > 0) {
					TParticle* mother = reinterpret_cast< TParticle* >(particles->At(mother1_id));
					if (mother->GetPdgCode() == 111 && mother->Pt() > 10.0) {
						Float_t isolation = 0;
						for (Long64_t k = 0; k < particles->GetEntries(); ++k) {
							TParticle* candidate = reinterpret_cast< TParticle* >(particles->At(k));
							//Bool_t cut2 = candidate->Eta() < 2.4 && candidate->Eta() > -2.4;
							if (candidate->GetStatusCode() > 0) {
								Float_t DR = DeltaR(mother,candidate); 
								TParticlePDG* candpdg = candidate->GetPDG(0);
								if (DR < 0.2 && k!=mother1_id && candpdg && k!=mother->GetDaughter(0) && k!=mother->GetDaughter(1)) {
									Float_t pt = candidate->Pt();
									Float_t mass = candidate->GetMass();
									Float_t Et = sqrt(pt*pt+mass*mass);
									isolation += Et;
									if (candpdg->Charge()!=0) isolation += pt;
								}
							}
						}
						if (isolation/mother->Pt() < 0.08) {
							h2d->Fill(selectedJets[i].pt(),mother->Pt()/selectedJets[i].pt());
							control = true;
						}
					}
				}
				if (mother2_id > 0 && !control) {
					TParticle* mother = reinterpret_cast< TParticle* >(particles->At(mother2_id));
					if (mother->GetPdgCode() == 111 && mother->Pt() > 10.0) {
						Float_t isolation = 0;
						for (Long64_t k = 0; k < particles->GetEntries(); ++k) {
							TParticle* candidate = reinterpret_cast< TParticle* >(particles->At(k));
							//Bool_t cut2 = candidate->Eta() < 2.4 && candidate->Eta() > -2.4;
							if (candidate->GetStatusCode() > 0) {
								Float_t DR = DeltaR(mother,candidate); 
								TParticlePDG* candpdg = candidate->GetPDG(0);
								if (DR < 0.2 && k!=mother2_id && candpdg && k!=mother->GetDaughter(0) && k!=mother->GetDaughter(1)) {
									Float_t pt = candidate->Pt();
									Float_t mass = candidate->GetMass();
									Float_t Et = sqrt(pt*pt+mass*mass);
									isolation += Et;
									if (candpdg->Charge()!=0) isolation += pt;
								}
							}
						}
						if (isolation/mother->Pt() < 0.08) {
							h2d->Fill(selectedJets[i].pt(),mother->Pt()/selectedJets[i].pt());
							control = true;
						}
					}
				}
				if (control) break;
			}
		}
	}

	TFile* fsave = new TFile("ptratio_vs_ptjet_high.root","RECREATE");
	// h1->Write();
	h2d->Write();
	fsave->Close();
	timer.Stop();
	timer.Print();
}
