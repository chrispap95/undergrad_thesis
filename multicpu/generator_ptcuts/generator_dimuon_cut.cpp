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

using namespace std;

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

void generator_dimuon_cut() {
	TStopwatch t;

	//Define all the histograms
	TCanvas* c2 = new TCanvas("c2", "trial analysis", 1);
	TH1I* hcountb = new TH1I("hcountb","B: counts for bkg;pT cut (GeV/c)",100, 10, 30);
	TH1I* hcounts = new TH1I("hcounts","S: counts for signal;pT cut (GeV/c)",100, 10, 30);
	TH1D* hsb = new TH1D("hsb","S/B vs  p_{T} cut for #dimuon;p_{T} cut (GeV/c); S/B",100, 10, 30);
	TH1D* hsbsq = new TH1D("hsbsq","S/#sqrt{B} vs p_{T} cut for #dimuon;p_{T} cut (GeV/c); S/#sqrt{B}",100, 10, 30);

	//Create array of counts
	Int_t countb[101]; //1-100 bins content
	Int_t counts[101]; //1-100 bins content
	for (Int_t i = 0; i < 101; ++i) {
		countb[i] = 0;
		counts[i] = 0;
	}
	
	///////////////////////////////////////////////////////////////////
	//                                                               //
	//  Insert the .root file and get the tree and its Event branch  //
	//                                                               //
	///////////////////////////////////////////////////////////////////
        TFile* f1 = TFile::Open("/data/cpapage/jpsiInclusive-gen8K.root");
        TTree* t1 = dynamic_cast< TTree* >(f1->Get("Pythia"));
        TClonesArray* particles1 = new TClonesArray("Event");
        t1->SetBranchAddress("Event", &particles1);

	// Set up the jet finder FastJet
	Double_t Rparam = 0.5;
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
	fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rparam, recombScheme, strategy);

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
		vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
		fastjet::ClusterSequence clustSeq(fjInputs, jetDef);

		// Output some info for the first run of the algorithm
		if (!i) {
			cout << "Ran " << jetDef.description() << endl;
			cout << "Strategy adopted by FastJet was "
			     << clustSeq.strategy_string() << endl << endl;
		}

		// Populate the output vector
		inclusiveJets = clustSeq.inclusive_jets(20.0);
		fastjet::Selector select_rapidity = fastjet::SelectorAbsEtaMax(2.4);
		fastjet::Selector select_all = select_rapidity;
		vector <fastjet::PseudoJet> selectedJets = select_all(inclusiveJets);

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
			TParticle* muon = reinterpret_cast< TParticle* >(particles1->At(iJpsi));
			Int_t iAntimuon = FindLastParticle(iJpsi, -13, particles1);
			TParticle* antimuon = reinterpret_cast< TParticle* >(particles1->At(iJpsi));
			Float_t energy = jpsi->Energy()+selectedJets[l].e();
			Float_t px = jpsi->Px()+selectedJets[l].px();
			Float_t py = jpsi->Py()+selectedJets[l].py();
			Float_t pz = jpsi->Pz()+selectedJets[l].pz();
			Float_t mass = sqrt(energy*energy-px*px-py*py-pz*pz);
			Bool_t cut_mass = mass > 122.5 && mass < 127.5;
			if (cut_mass) {
				Float_t pT1 = muon->Pt();
				Float_t pT2 = antimuon->Pt();
				Float_t pTmax;
				if (pT1 > pT2) pTmax = pT1;
				else pTmax = pT2;
				Float_t dbin = (pTmax-10)/0.2;
				for (Int_t ibin = 0; ibin < (int)dbin && ibin < 101; ++ibin) {
					++countb[ibin];
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////
	//                                                               //
	//  Insert the .root file and get the tree and its Event branch  //
	//                                                               //
	///////////////////////////////////////////////////////////////////
        TFile* f2 = TFile::Open("/data/cpapage/higgsTojpsiGamma-gen.root");
        TTree* t2 = dynamic_cast< TTree* >(f2->Get("Pythia"));
        TClonesArray* particles2 = new TClonesArray("Event");
        t2->SetBranchAddress("Event", &particles2);
        for (Long64_t i = 0; i < t2->GetEntries(); ++i) { // This is like the Event loop 
                t2->GetEntry(i);
		
		Int_t iHiggs = 0;
		for (Long64_t l = 0; l < particles2->GetEntries(); ++l) { //This is a particle loop
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

		TParticle* gamma = reinterpret_cast< TParticle* >(particles2->At(iGamma));
		TParticle* muon = reinterpret_cast< TParticle* >(particles2->At(iMuon));
		TParticle* antimuon = reinterpret_cast< TParticle* >(particles2->At(iAntimuon));
		Bool_t cut1 = gamma->Eta() < 2.4 && gamma->Eta() > -2.4 && gamma->Pt() > 20;
		Bool_t cut2 = muon->Eta() < 2.4 && muon->Eta() > -2.4;
		Bool_t cut3 = antimuon->Eta() < 2.4 && antimuon->Eta() > -2.4;
		if (cut1) {
			Float_t pT1 = muon->Pt();
			Float_t pT2 = antimuon->Pt();
			Float_t pTmax;
			if (pT1 > pT2) pTmax = pT1;
			else pTmax = pT2;
			Float_t dbin = (pTmax-10)/0.2;
			for (Int_t ibin = 0; ibin < (int)dbin && ibin < 101; ++ibin) {
				++counts[ibin];
			}	
		}	
	}

	//Normalization for 100K events of bkg and 1K events of signal
	Int_t countb_norm[101];
	Int_t counts_norm[101];
	for (Int_t i = 1; i<101; ++i) {
		countb_norm[i] = countb[i];
		counts_norm[i] = counts[i];
	}
	float normb = 169800*0.0035*100/8;
	float norms = 0.027306;
	for (Int_t i = 1; i<101; ++i) {
		countb_norm[i] *= normb;
		counts_norm[i] *= norms;
	}

	for (Int_t i = 1; i<101; ++i) {
		if (countb_norm[i] >= 1) {
			hsb->SetBinContent(i, (double)counts_norm[i]/(double)countb_norm[i]);
			hsbsq->SetBinContent(i, counts_norm[i]/sqrt(countb_norm[i]));
			hcounts->SetBinContent(i, counts[i]);
			hcountb->SetBinContent(i, countb[i]);
		}
		else {
			hsb->SetBinContent(i, 0);
			hsbsq->SetBinContent(i, 0);
			hcounts->SetBinContent(i, counts[i]);
			hcountb->SetBinContent(i, countb[i]);
		}
	}

	c2->Divide(2,2);
	c2->cd(1); 
	hcounts->Draw();
	c2->cd(2);
	hcountb->Draw();
	c2->cd(3);
	hsb->Draw();
	c2->cd(4);
	hsbsq->Draw();
        c2->Print("generator_dimuon_cut_8K.pdf");

	t.Stop();
	t.Print();
}
