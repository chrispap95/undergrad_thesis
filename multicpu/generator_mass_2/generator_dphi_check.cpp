#include "/expsoft/fastjet-3.3.0-cms/include/fastjet/PseudoJet.hh"
#include "/expsoft/fastjet-3.3.0-cms/include/fastjet/ClusterSequence.hh"
#ifndef __CINT__
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TRandom3.h>
#endif

Int_t FindLastParticle(Int_t iPart, Int_t pdgCode, TClonesArray* particles) { //Find the last Jpsi
	while (1) {
		TParticle* part = reinterpret_cast< TParticle* >(particles->At(iPart));
		Int_t daughter0id = part->GetDaughter(0);
		Int_t daughter1id = part->GetDaughter(1);
		Bool_t daughter0validity = true;
		Bool_t daughter1validity = true;
		if (daughter0id >= 0) {
			TParticle* daughter0 = reinterpret_cast< TParticle* >(particles->At(daughter0id));
			if (daughter0->GetPdgCode() == pdgCode) {
				iPart = part->GetDaughter(0);
				daughter0validity = false;
				continue;
			}
		}
		if (daughter1id >= 0) {
			TParticle* daughter1 = reinterpret_cast< TParticle* >(particles->At(daughter1id));
			if (daughter1->GetPdgCode() == pdgCode) {
				iPart = part->GetDaughter(1);
				daughter1validity = false;
				continue;
			}
		}
		if (daughter0validity && daughter1validity) break;
	}
	return iPart;
}

Float_t bkg_mass(string &file, Float_t countb, TH1D* h01b, TH1D* h01bc, TH1F* hpy1, TH1F* hpy2, TH1F* hpy3, TH1F* hpy4) {
	// Set up the jet finder FastJet
	Double_t Rparam = 0.5;
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
	fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rparam, recombScheme, strategy);

	// Insert the .root file and get the tree and its Event branch
        TFile* f1 = TFile::Open(file.c_str());
        TTree* t1 = dynamic_cast< TTree* >(f1->Get("Pythia"));
        TClonesArray* particles1 = new TClonesArray("Event");
        t1->SetBranchAddress("Event", &particles1);

	// This is the Event loop
        for (Long64_t i = 0; i < t1->GetEntries(); ++i) {
                t1->GetEntry(i);
		vector<fastjet::PseudoJet> fjInputs;

		// This is the particle loop
		for (Int_t ip = 0; ip < particles1->GetEntries(); ++ip) {
			TParticle* part = reinterpret_cast< TParticle* >(particles1->At(ip));
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
		if (!i) {
			cout << "Ran " << jetDef.description() << endl;
			cout << "Strategy adopted by FastJet was "
			     << clustSeq.strategy_string() << endl << endl;
		}

		// Populate the output vector
		inclusiveJets = clustSeq.inclusive_jets(30.0);
		fastjet::Selector select_rapidity = fastjet::SelectorAbsEtaMax(2.37);
		vector <fastjet::PseudoJet> selectedJets = select_rapidity(inclusiveJets);

		// This is a particle loop
                for (Long64_t l = 0; l < selectedJets.size(); ++l) {
			Int_t iJpsi = 0;
			for (Int_t k = 0; k < particles1->GetEntries(); ++k) {
				TParticle* jpsi = reinterpret_cast< TParticle* >(particles1->At(k));
				if (jpsi->GetPdgCode() == 443) {
					iJpsi = k;
					break;
				}
			}
			iJpsi = FindLastParticle(iJpsi, 443, particles1);
			TParticle* jpsi = reinterpret_cast< TParticle* >(particles1->At(iJpsi));
			Int_t iMuon = FindLastParticle(iJpsi, 13, particles1);
			TParticle* muon = reinterpret_cast< TParticle* >(particles1->At(iMuon));
			Int_t iAntimuon = FindLastParticle(iJpsi, -13, particles1);
			TParticle* antimuon = reinterpret_cast< TParticle* >(particles1->At(iAntimuon));
			Bool_t cut1 = muon->Eta() < 2.5 && muon->Eta() > -2.5;
			Bool_t cut2 = antimuon->Eta() < 2.5 && antimuon->Eta() > -2.5;
			Bool_t cut3 = muon->Pt() > 10 && antimuon->Pt() > 10;
			Bool_t cut4 = muon->Pt() > 20 || antimuon->Pt() > 20;
			Float_t deltaPhi = jpsi->Phi()-selectedJets[l].phi();
			Bool_t cut5 = selectedJets[l].pt() > 36;
			Bool_t cut6 = jpsi->Pt() > 36;
			//Isolation cut for leading muon needed here

			// Momentum fraction of pi0
			Float_t fraction = 0;
			Float_t weight = 0;
			if (selectedJets[l].pt() < 40) {
				fraction = hpy1->GetRandom();
				weight = 0.000431215;
			}
			else if (selectedJets[l].pt() < 50) {
				fraction = hpy2->GetRandom();
				weight = 0.000431215;
			}
			else if (selectedJets[l].pt() < 60) {
				fraction = hpy3->GetRandom();
				weight = 0.000366135;
			}
			else if (selectedJets[l].pt() < 70) {
				fraction = hpy4->GetRandom();
				weight = 0.000366135;
			}
			else {
				fraction = hpy4->GetRandom();
				weight = 0.000155525;
			}
			if (cut1 && cut2) {
				h01b->Fill(deltaPhi,weight);
				if (cut2 && cut3 && cut4 && cut5 && cut6) {
					h01bc->Fill(deltaPhi,weight);
				}
			}
		}
	}
	return countb;
}

void generator_dphi_check() {
	TStopwatch t;
	TRandom3 r(0);

	Double_t countb1 = 0;
	Double_t countb2 = 0;
	Double_t countb3 = 0;
	Double_t countb4 = 0;
        Double_t counts = 0;

	// Define all the histograms
	TCanvas* c2 = new TCanvas("c2", "trial analysis", 1);
	
	TH1D* h01b1 = new TH1D("h01b1","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01b2 = new TH1D("h01b2","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01b3 = new TH1D("h01b3","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01b4 = new TH1D("h01b4","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01s = new TH1D("h01s","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01b1c = new TH1D("h01b1c","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01b2c = new TH1D("h01b2c","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01b3c = new TH1D("h01b3c","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01b4c = new TH1D("h01b4c","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());
	TH1D* h01sc = new TH1D("h01sc","#Delta#phi",100, -2*TMath::Pi(),2*TMath::Pi());

	TFile* fdist = TFile::Open("../bkg/generator_jets/ptratio_distributions.root");
	TH1F* hpy1 = (TH1F*)(fdist->Get("hpy1"));
	TH1F* hpy2 = (TH1F*)(fdist->Get("hpy2"));
	TH1F* hpy3 = (TH1F*)(fdist->Get("hpy3"));
	TH1F* hpy4 = (TH1F*)(fdist->Get("hpy4"));

	string input = "/data/cpapage/jpsiInclusive-gen100K-pt30-60.root";
	countb1 = bkg_mass(input, 0, h01b1, h01b1c, hpy1, hpy2, hpy3, hpy4);
	input = "/data/cpapage/jpsiInclusive-gen100K-pt60-90.root";
	countb2 = bkg_mass(input, 0, h01b2, h01b2c, hpy1, hpy2, hpy3, hpy4);
	input = "/data/cpapage/jpsiInclusive-gen100K-pt90-120.root";
	countb3 = bkg_mass(input, 0, h01b3, h01b3c, hpy1, hpy2, hpy3, hpy4);
	input = "/data/cpapage/jpsiInclusive-gen100K-pt120.root";
	countb4 = bkg_mass(input, 0, h01b4, h01b4c, hpy1, hpy2, hpy3, hpy4);

	// Insert the .root file and get the tree and its Event branch
        TFile* f2 = TFile::Open("/data/cpapage/higgsTojpsiGamma-gen.root");
        TTree* t2 = dynamic_cast< TTree* >(f2->Get("Pythia"));
        TClonesArray* particles2 = new TClonesArray("Event");
        t2->SetBranchAddress("Event", &particles2);

	// This is like the Event loop
        for (Long64_t i = 0; i < t2->GetEntries(); ++i) {
                t2->GetEntry(i);
		Int_t iHiggs = 0;

		// This is a particle loop
		for (Long64_t l = 0; l < particles2->GetEntries(); ++l) {
			TParticle* part = reinterpret_cast< TParticle* >(particles2->At(l));
			if (part->GetPdgCode() == 25) {
				iHiggs = l;
				break;
			}
		}
		iHiggs = FindLastParticle(iHiggs, 25, particles2);
		Int_t iGamma = iHiggs;
		iGamma = FindLastParticle(iGamma, 22, particles2);
		Int_t iJpsi = iHiggs;
		iJpsi = FindLastParticle(iJpsi, 443, particles2);
		Int_t iMuon = iJpsi;
		iMuon = FindLastParticle(iMuon, 13, particles2);
		Int_t iAntimuon = iJpsi;
		iAntimuon = FindLastParticle(iAntimuon, -13, particles2);

		TParticle* jpsi = reinterpret_cast< TParticle* >(particles2->At(iJpsi));
		TParticle* muon = reinterpret_cast< TParticle* >(particles2->At(iMuon));
		TParticle* antimuon = reinterpret_cast< TParticle* >(particles2->At(iAntimuon));
		TParticle* gamma = reinterpret_cast< TParticle* >(particles2->At(iGamma));
	 	Bool_t cut1 = gamma->Eta() < 2.37 && gamma->Eta() > -2.37;
		Bool_t cut2 = muon->Eta() < 2.5 && muon->Eta() > -2.5;
		Bool_t cut3 = antimuon->Eta() < 2.5 && antimuon->Eta() > -2.5;
		Float_t deltaPhi = jpsi->Phi()-gamma->Phi();
		Bool_t cut7 = muon->Pt() > 10 && antimuon->Pt() > 10;
		Bool_t cut4 = muon->Pt() > 20 || antimuon->Pt() > 20;
		Bool_t cut5 = gamma->Pt() > 36;
		Bool_t cut6 = jpsi->Pt() > 36;
		
		if (cut1 && cut2 && cut3) {
			h01s->Fill(deltaPhi);
			if (cut4 && cut5 && cut6 && cut7) {
				h01sc->Fill(deltaPhi);
			}
		}
	}

	TFile* f = new TFile("dphi.root","recreate");
	h01b1->Write();
	h01b2->Write();
	h01b3->Write();
	h01b4->Write();
	h01b1c->Write();
	h01b2c->Write();
	h01b3c->Write();
	h01b4c->Write();
	h01s->Write();
	h01sc->Write();
	t.Stop();
	t.Print();
}
