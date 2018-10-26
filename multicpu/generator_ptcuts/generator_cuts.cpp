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

void generator_bkg_cuts(string &input,TH1D* hcountb1,TH1D* hcountb2,TH1D* hcountb3,TH1D* hcountb4,Float_t normb) {
	//Create array of counts
	Double_t countb1[101],countb2[101],countb3[101],countb4[101]; //1-100 bins content
	for (Int_t i = 0; i < 101; ++i) {
		countb1[i] = 0;
		countb2[i] = 0;
		countb3[i] = 0;
		countb4[i] = 0;
	}
	
	///////////////////////////////////////////////////////////////////
	//                                                               //
	//  Insert the .root file and get the tree and its Event branch  //
	//                                                               //
	///////////////////////////////////////////////////////////////////
        TFile* f1 = TFile::Open(input.c_str());
        TTree* t1 = dynamic_cast< TTree* >(f1->Get("Pythia"));
        TClonesArray* particles1 = new TClonesArray("Event");
        t1->SetBranchAddress("Event", &particles1);

	// Set up the jet finder FastJet
	Double_t Rparam = 0.5;
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
	fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, Rparam, recombScheme, strategy);

	TFile* fdist = TFile::Open("../bkg/generator_jets/ptratio_distributions.root");
	TH1F* hpy1 = (TH1F*)(fdist->Get("hpy1"));
	TH1F* hpy2 = (TH1F*)(fdist->Get("hpy2"));
	TH1F* hpy3 = (TH1F*)(fdist->Get("hpy3"));
	TH1F* hpy4 = (TH1F*)(fdist->Get("hpy4"));

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
		fastjet::Selector select_rapidity = fastjet::SelectorAbsEtaMax(2.37);
		fastjet::Selector select_all = select_rapidity;
		vector <fastjet::PseudoJet> selectedJets = select_all(inclusiveJets);

		// This is a particle loop
                for (Long64_t l = 0; l < selectedJets.size(); ++l) {
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
			Int_t iMuon = iJpsi;
			iMuon = FindLastParticle(iMuon, 13, particles1);
			TParticle* muon = reinterpret_cast< TParticle* >(particles1->At(iMuon));
			Int_t iAntimuon = iJpsi;
			iAntimuon = FindLastParticle(iAntimuon, -13, particles1);
			TParticle* antimuon = reinterpret_cast< TParticle* >(particles1->At(iAntimuon));

			Float_t energy = jpsi->Energy()+fraction*selectedJets[l].e();
			Float_t px = jpsi->Px()+fraction*selectedJets[l].px();
			Float_t py = jpsi->Py()+fraction*selectedJets[l].py();
			Float_t pz = jpsi->Pz()+fraction*selectedJets[l].pz();
			Float_t mass = sqrt(energy*energy-px*px-py*py-pz*pz);

			Bool_t cut_mass = mass > 123 && mass < 127;
			Bool_t cut2 = muon->Eta() < 2.5 && muon->Eta() > -2.5;
			Bool_t cut3 = antimuon->Eta() < 2.5 && antimuon->Eta() > -2.5;
			Bool_t cut4 = abs(jpsi->Phi()-selectedJets[l].phi()) > 0.5;

			if (cut_mass && cut2 && cut3 && cut4) {
				Float_t pT1 = fraction*selectedJets[l].pt();
				Float_t pT2,pT3;
				if (muon->Pt() > antimuon->Pt()) {
					pT2 = muon->Pt();
					pT3 = antimuon->Pt();
				}else {
					pT2 = antimuon->Pt();
					pT3 = muon->Pt();
				}
				Float_t pT4 = jpsi->Pt();
				Float_t dbin1 = (pT1-20)/0.4;
				Float_t dbin2 = (pT2)/0.3;
				Float_t dbin3 = (pT3)/0.2;
				Float_t dbin4 = (pT4)/0.5;

				Bool_t cut_photon = pT1 > 58;
				Bool_t cut_muon_low = 1;//pT2 > 8.1 && pT3 > 8.1;
				Bool_t cut_muon_high = pT2 > 15 || pT3 > 15;
				Bool_t cut_jpsi = 1;//pT4 > 30.5;

				if (cut_muon_low && cut_jpsi && cut_muon_high) { 
					for (Int_t ibin = 0; ibin < (int)dbin1 && ibin < 101; ++ibin) {
						countb1[ibin] += weight;
					}
				}
				if (cut_photon && cut_muon_low && cut_jpsi) {
					for (Int_t ibin = 0; ibin < (int)dbin2 && ibin < 101; ++ibin) {
						countb2[ibin] += weight;
					}
				}
				if (cut_photon && cut_jpsi && cut_muon_high) {
					for (Int_t ibin = 0; ibin < (int)dbin3 && ibin < 101; ++ibin) {
						countb3[ibin] += weight;
					}
				}
				if (cut_photon && cut_muon_low && cut_muon_high) { 
					for (Int_t ibin = 0; ibin < (int)dbin4 && ibin < 101; ++ibin) {
						countb4[ibin] += weight;
					}
				}
			}
		}
	}

	Double_t countb1_norm[101],countb2_norm[101],countb3_norm[101],countb4_norm[101];
	for (Int_t i = 1; i<101; ++i) {
		countb1_norm[i] = countb1[i];
		countb2_norm[i] = countb2[i];
		countb3_norm[i] = countb3[i];
		countb4_norm[i] = countb4[i];
	}
	for (Int_t i = 1; i<101; ++i) {
		countb1_norm[i] *= normb;
		countb2_norm[i] *= normb;
		countb3_norm[i] *= normb;
		countb4_norm[i] *= normb;
	}
	for (Int_t i = 1; i<101; ++i) {
		hcountb1->SetBinContent(i, countb1_norm[i]);
		hcountb2->SetBinContent(i, countb2_norm[i]);
		hcountb3->SetBinContent(i, countb3_norm[i]);
		hcountb4->SetBinContent(i, countb4_norm[i]);
	}
}

void generator_cuts() {
	TStopwatch t;

	//Define all the histograms
	TCanvas* c2 = new TCanvas("c2", "trial analysis", 1);
	//Gamma pT cut
	TH1D* hcountb10 = new TH1D("hcountb10","B: counts for bkg;pT cut (GeV/c)",100, 20, 60);
	TH1D* hcountb11 = new TH1D("hcountb11","B: counts for bkg;pT cut (GeV/c)",100, 20, 60);
	TH1D* hcountb12 = new TH1D("hcountb12","B: counts for bkg;pT cut (GeV/c)",100, 20, 60);
	TH1D* hcountb13 = new TH1D("hcountb13","B: counts for bkg;pT cut (GeV/c)",100, 20, 60);
	TH1D* hcounts1 = new TH1D("hcounts1","S: counts for signal;pT cut (GeV/c)",100, 20, 60);

	//High pT Mu pT cut 
	TH1D* hcountb20 = new TH1D("hcountb20","B: counts for bkg;pT cut (GeV/c)",100, 0, 30);
	TH1D* hcountb21 = new TH1D("hcountb21","B: counts for bkg;pT cut (GeV/c)",100, 0, 30);
	TH1D* hcountb22 = new TH1D("hcountb22","B: counts for bkg;pT cut (GeV/c)",100, 0, 30);
	TH1D* hcountb23 = new TH1D("hcountb23","B: counts for bkg;pT cut (GeV/c)",100, 0, 30);
	TH1D* hcounts2 = new TH1D("hcounts2","S: counts for signal;pT cut (GeV/c)",100, 0, 30);

	//Low pT Mu pT cut
	TH1D* hcountb30 = new TH1D("hcountb30","B: counts for bkg;pT cut (GeV/c)",100, 0, 20);
	TH1D* hcountb31 = new TH1D("hcountb31","B: counts for bkg;pT cut (GeV/c)",100, 0, 20);
	TH1D* hcountb32 = new TH1D("hcountb32","B: counts for bkg;pT cut (GeV/c)",100, 0, 20);
	TH1D* hcountb33 = new TH1D("hcountb33","B: counts for bkg;pT cut (GeV/c)",100, 0, 20);
	TH1D* hcounts3 = new TH1D("hcounts3","S: counts for signal;pT cut (GeV/c)",100, 0, 20);

	//J/psi pT cut
	TH1D* hcountb40 = new TH1D("hcountb40","B: counts for bkg;pT cut (GeV/c)",100, 0, 50);
	TH1D* hcountb41 = new TH1D("hcountb41","B: counts for bkg;pT cut (GeV/c)",100, 0, 50);
	TH1D* hcountb42 = new TH1D("hcountb42","B: counts for bkg;pT cut (GeV/c)",100, 0, 50);
	TH1D* hcountb43 = new TH1D("hcountb43","B: counts for bkg;pT cut (GeV/c)",100, 0, 50);
	TH1D* hcounts4 = new TH1D("hcounts4","S: counts for signal;pT cut (GeV/c)",100, 0, 50);

	///////////////////////////////////////////////////////////////////
	//  Analyze Bkg  //////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	string input = "/data/cpapage/jpsiInclusive-gen100K-pt30-60.root";
	generator_bkg_cuts(input,hcountb10,hcountb20,hcountb30,hcountb40,28010.88);
	input = "/data/cpapage/jpsiInclusive-gen100K-pt60-90.root";
	generator_bkg_cuts(input,hcountb11,hcountb21,hcountb31,hcountb41,1371.25);
	input = "/data/cpapage/jpsiInclusive-gen100K-pt90-120.root";
	generator_bkg_cuts(input,hcountb12,hcountb22,hcountb32,hcountb42,205);
	input = "/data/cpapage/jpsiInclusive-gen100K-pt120.root";
	generator_bkg_cuts(input,hcountb13,hcountb23,hcountb33,hcountb43,75.55);

	hcountb10->Add(hcountb11); hcountb10->Add(hcountb12); hcountb10->Add(hcountb13);
	hcountb20->Add(hcountb21); hcountb20->Add(hcountb22); hcountb20->Add(hcountb23);
	hcountb30->Add(hcountb31); hcountb30->Add(hcountb32); hcountb30->Add(hcountb33);
	hcountb40->Add(hcountb41); hcountb40->Add(hcountb42); hcountb40->Add(hcountb43);

	///////////////////////////////////////////////////////////////////
	//  Analyze Signal  ///////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	//                                                               //
	//  Insert the .root file and get the tree and its Event branch  //
	//                                                               //
	///////////////////////////////////////////////////////////////////
	Double_t counts1[101],counts2[101],counts3[101],counts4[101]; //1-100 bins content
	for (Int_t i = 0; i < 101; ++i) {
		counts1[i] = 0;
		counts2[i] = 0;
		counts3[i] = 0;
		counts4[i] = 0;
	}


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

		TParticle* jpsi = reinterpret_cast< TParticle* >(particles2->At(iJpsi));
		TParticle* gamma = reinterpret_cast< TParticle* >(particles2->At(iGamma));
		TParticle* muon = reinterpret_cast< TParticle* >(particles2->At(iMuon));
		TParticle* antimuon = reinterpret_cast< TParticle* >(particles2->At(iAntimuon));

		Bool_t cut1 = gamma->Eta() < 2.37 && gamma->Eta() > -2.37 && gamma->Pt() > 20;
		Bool_t cut2 = muon->Eta() < 2.5 && muon->Eta() > -2.5;
		Bool_t cut3 = antimuon->Eta() < 2.5 && antimuon->Eta() > -2.5;
		Bool_t cut4 = abs(jpsi->Phi()-gamma->Phi()) > 0.5;

		if (cut1 && cut2 && cut3 && cut4) {
			Float_t pT1 = gamma->Pt();
			Float_t pT2,pT3;
			if (muon->Pt() > antimuon->Pt()) {
				pT2 = muon->Pt();
				pT3 = antimuon->Pt();
			}else {
				pT2 = antimuon->Pt();
				pT3 = muon->Pt();
			}
			Float_t pT4 = jpsi->Pt();
			Float_t dbin1 = (pT1-20)/0.4;
			Float_t dbin2 = (pT2)/0.3;
			Float_t dbin3 = (pT3)/0.2;
			Float_t dbin4 = (pT4)/0.5;

			Bool_t cut_photon = pT1 > 58;
			Bool_t cut_muon_low = 1;//pT2 > 8.1 && pT3 > 8.1;
			Bool_t cut_muon_high = pT2 > 15 || pT3 > 15;
			Bool_t cut_jpsi = 1;//pT4 > 30.5;

			if (cut_muon_low && cut_jpsi && cut_muon_high) { 
				for (Int_t ibin = 0; ibin < (int)dbin1 && ibin < 101; ++ibin) {
					++counts1[ibin];
				}
			}
			if (cut_photon && cut_muon_low && cut_jpsi) {
				for (Int_t ibin = 0; ibin < (int)dbin2 && ibin < 101; ++ibin) {
					++counts2[ibin];
				}
			}
			if (cut_photon && cut_jpsi && cut_muon_high) {
				for (Int_t ibin = 0; ibin < (int)dbin3 && ibin < 101; ++ibin) {
					++counts3[ibin];
				}
			}
			if (cut_photon && cut_muon_low && cut_muon_high) { 
				for (Int_t ibin = 0; ibin < (int)dbin4 && ibin < 101; ++ibin) {
					++counts4[ibin];
				}
			}
		}
	}

	//Normalization for 1K events of signal
	Double_t counts1_norm[101],counts2_norm[101],counts3_norm[101],counts4_norm[101];
	for (Int_t i = 1; i<101; ++i) {
		counts1_norm[i] = counts1[i];
		counts2_norm[i] = counts2[i];
		counts3_norm[i] = counts3[i];
		counts4_norm[i] = counts4[i];
	}
	float norms = 0.02913;
	for (Int_t i = 1; i<101; ++i) {
		counts1_norm[i] *= norms;
		counts2_norm[i] *= norms;
		counts3_norm[i] *= norms;
		counts4_norm[i] *= norms;
	}
	for (Int_t i = 1; i<101; ++i) {
		hcounts1->SetBinContent(i, counts1_norm[i]);
		hcounts2->SetBinContent(i, counts2_norm[i]);
		hcounts3->SetBinContent(i, counts3_norm[i]);
		hcounts4->SetBinContent(i, counts4_norm[i]);
	}

	//Save the results
	TFile* fout = new TFile("cuts.root","recreate");
	hcounts1->Write(); //Cut 1
	hcountb10->Write();
	hcounts2->Write(); //Cut 2
	hcountb20->Write();
	hcounts3->Write(); //Cut 3
	hcountb30->Write();
	hcounts4->Write(); //Cut 4
	hcountb40->Write();
	fout->Close();

	t.Stop();
	t.Print();
}
