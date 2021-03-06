#ifndef __CINT__
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParticle.h>
#endif

void generator_routine_1() {
	//Define all the histograms
	TCanvas* c1 = new TCanvas("c1", "trial analysis", 1);
	TH1F* h10 = new TH1F("h10", "p_{T} of H; p_{T} (GeV/c)", 100, 0, 400);
	TH1F* h20 = new TH1F("h20", "y of H; y", 100, -5, 5);
	TH1F* h30 = new TH1F("h30", "#eta of H; #eta", 100, -10, 10);
	TH1F* h40 = new TH1F("h40", "E of H; E(GeV)", 100, 0, 2000);
	TH1F* h11 = new TH1F("h11", "p_{T} of J/#psi; p_{T} (GeV/c)", 100, 0, 200);
	TH1F* h21 = new TH1F("h21", "y of J/#psi; y", 100, -6, 6);
	TH1F* h31 = new TH1F("h31", "#eta of J/#psi; #eta", 100, -8, 8);
	TH1F* h41 = new TH1F("h41", "E of J/#psi; E(GeV)", 100, 0, 1000);
	TH1F* h12 = new TH1F("h12", "p_{T} of #gamma; p_{T} (GeV/c)", 100, 0, 200);
	TH1F* h22 = new TH1F("h22", "y of #gamma; y", 100, -6, 6);
	TH1F* h32 = new TH1F("h32", "#eta of #gamma; #eta", 100, -6, 6);
	TH1F* h42 = new TH1F("h42", "E of #gamma; E(GeV)", 100, 0, 1000);
	TH1F* h13 = new TH1F("h13", "p_{T} of #mu; p_{T} (GeV/c)", 100, 0, 150);
	TH1F* h23 = new TH1F("h23", "y of #mu; y", 100, -8, 8);
	TH1F* h33 = new TH1F("h33", "#eta of #mu; #eta", 100, -8, 8);
	TH1F* h43 = new TH1F("h43", "E of #mu; E (GeV)", 100, 0, 600);
	TH1F* h00 = new TH1F("h00", "P_{J/#psi#gamma}-P_{H}; #Delta P (GeV/c)", 100, -4, 4);
	TH1F* h01 = new TH1F("h01", "m_{J/#psi#gamma}; m_{J/#psi#gamma} (GeV/c^{2})", 100, 0, 250);
	
	//Insert the .root file and get the tree and its Event branch
        TFile* f = TFile::Open("/data/cpapage/higgsTojpsiGamma-gen.root");
        TTree* t = dynamic_cast< TTree* >(f->Get("Pythia"));
        TClonesArray* particles = new TClonesArray("Event");
        t->SetBranchAddress("Event", &particles);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) { // This is like the Event loop 
                t->GetEntry(i);

		Int_t iHiggs = 0;
		Int_t iJpsi = 0;
		Int_t iGamma = 0;
		Int_t iMuon = 0;
		Int_t iAntiMuon = 0;

                for (Long64_t l = 0; l < particles->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles->At(l));
			if (part->GetPdgCode() == 25) {
				iHiggs = l;
				break;
			}
                }
		while (1) {    //Find the last Higgs
			TParticle* part = reinterpret_cast< TParticle* >(particles->At(iHiggs));
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = reinterpret_cast< TParticle* >(particles->At(daughter0id));
				if (daughter0->GetPdgCode() == 25) {
					iHiggs = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = reinterpret_cast< TParticle* >(particles->At(daughter1id));
				if (daughter1->GetPdgCode() == 25) {
					iHiggs = part->GetDaughter(1);
					daughter1validity = false;
					continue;
				}
			}
			if (daughter0validity && daughter1validity) break;
		}

		iJpsi = iHiggs;
		while (1) {    //Find the last J/psi
			TParticle* part = reinterpret_cast< TParticle* >(particles->At(iJpsi));
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = reinterpret_cast< TParticle* >(particles->At(daughter0id));
				if (daughter0->GetPdgCode() == 443) {
					iJpsi = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = reinterpret_cast< TParticle* >(particles->At(daughter1id));
				if (daughter1->GetPdgCode() == 443) {
					iJpsi = part->GetDaughter(1);
					daughter1validity = false;
					continue;
				}
			}
			if (daughter0validity && daughter1validity) break;
		}
		
		iGamma = iHiggs;
		while (1) {    //Find the last gamma
			TParticle* part = reinterpret_cast< TParticle* >(particles->At(iGamma));
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = reinterpret_cast< TParticle* >(particles->At(daughter0id));
				if (daughter0->GetPdgCode() == 22) {
					iGamma = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = reinterpret_cast< TParticle* >(particles->At(daughter1id));
				if (daughter1->GetPdgCode() == 22) {
					iGamma = part->GetDaughter(1);
					daughter1validity = false;
					continue;
				}
			}
			if (daughter0validity && daughter1validity) break;
		}

		iMuon = iJpsi;
		while (1) {    //Find the last Muon
			TParticle* part = reinterpret_cast< TParticle* >(particles->At(iMuon));
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = reinterpret_cast< TParticle* >(particles->At(daughter0id));
				if (daughter0->GetPdgCode() == 13) {
					iMuon = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = reinterpret_cast< TParticle* >(particles->At(daughter1id));
				if (daughter1->GetPdgCode() == 13) {
					iMuon = part->GetDaughter(1);
					daughter1validity = false;
					continue;
				}
			}
			if (daughter0validity && daughter1validity) break;
		}

		iAntiMuon = iJpsi;
		while (1) {    //Find the last Antimuon
			TParticle* part = reinterpret_cast< TParticle* >(particles->At(iAntiMuon));
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = reinterpret_cast< TParticle* >(particles->At(daughter0id));
				if (daughter0->GetPdgCode() == -13) {
					iAntiMuon = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = reinterpret_cast< TParticle* >(particles->At(daughter1id));
				if (daughter1->GetPdgCode() == -13) {
					iAntiMuon = part->GetDaughter(1);
					daughter1validity = false;
					continue;
				}
			}
			if (daughter0validity && daughter1validity) break;
		}

		TParticle* higgs = reinterpret_cast< TParticle* >(particles->At(iHiggs));
		TParticle* jpsi = reinterpret_cast< TParticle* >(particles->At(iJpsi));
		TParticle* gamma = reinterpret_cast< TParticle* >(particles->At(iGamma));
		TParticle* muon = reinterpret_cast< TParticle* >(particles->At(iMuon));
		TParticle* antimuon = reinterpret_cast< TParticle* >(particles->At(iAntiMuon));

		Float_t pT_higgs = higgs->Pt();
		h10->Fill(pT_higgs);
		Float_t y_higgs = higgs->Y();
		h20->Fill(y_higgs);
		Float_t eta_higgs = higgs->Eta();
		h30->Fill(eta_higgs);
		Float_t E_higgs = higgs->Energy();
		h40->Fill(E_higgs);

		Float_t pT_jpsi = jpsi->Pt();
		h11->Fill(pT_jpsi);
		Float_t y_jpsi = jpsi->Y();
		h21->Fill(y_jpsi);
		Float_t eta_jpsi = jpsi->Eta();
		h31->Fill(eta_jpsi);
		Float_t E_jpsi = jpsi->Energy();
		h41->Fill(E_jpsi);

		Float_t pT_gamma = gamma->Pt();
		h12->Fill(pT_gamma);
		Float_t y_gamma = gamma->Y();
		h22->Fill(y_gamma);
		Float_t eta_gamma = gamma->Eta();
		h32->Fill(eta_gamma);
		Float_t E_gamma = gamma->Energy();
		h42->Fill(E_gamma);

		Float_t pT_muon = muon->Pt();
		h13->Fill(pT_muon);
		Float_t y_muon = muon->Y();
		h23->Fill(y_muon);
		Float_t eta_muon = muon->Eta();
		h33->Fill(eta_muon);
		Float_t E_muon = muon->Energy();
		h43->Fill(E_muon);

		Float_t pT_antimuon = antimuon->Pt();
		h13->Fill(pT_antimuon);
		Float_t y_antimuon = antimuon->Y();
		h23->Fill(y_antimuon);
		Float_t eta_antimuon = antimuon->Eta();
		h33->Fill(eta_antimuon);
		Float_t E_antimuon = antimuon->Energy();
		h43->Fill(E_antimuon);
		
		Float_t px_gamma = gamma->Px();
		Float_t py_gamma = gamma->Py();
		Float_t pz_gamma = gamma->Pz();
		Float_t px_jpsi = jpsi->Px();
		Float_t py_jpsi = jpsi->Py();
		Float_t pz_jpsi = jpsi->Pz();
		Float_t p_higgs = higgs->P();
		Float_t px_jpsigamma = px_gamma+px_jpsi;
		Float_t py_jpsigamma = py_gamma+py_jpsi;
		Float_t pz_jpsigamma = pz_gamma+pz_jpsi;
		Float_t p_jpsigamma = sqrt(px_jpsigamma*px_jpsigamma+py_jpsigamma*py_jpsigamma+pz_jpsigamma*pz_jpsigamma);
		h00->Fill(p_jpsigamma-p_higgs);
		h01->Fill(sqrt((E_gamma+E_jpsi)*(E_gamma+E_jpsi)-p_jpsigamma*p_jpsigamma));
        }

	TText* t1 = new TText(120,40,"nEvents = 1000");
	TText* t2 = new TText(120,50,"Pythia generator results");
	TLatex t3;
	
	c1->Divide(2,2);
	c1->cd(1); h10->Draw();
	c1->cd(2); h20->Draw();
	c1->cd(3); h30->Draw();
	c1->cd(4); h40->Draw();
	c1->cd(1); t1->Draw(); t2->Draw();
	c1->Print("check_hists.pdf(");
	c1->cd(1); h11->Draw();
	c1->cd(2); h21->Draw();
	c1->cd(3); h31->Draw();
	c1->cd(4); h41->Draw();
	c1->Print("check_hists.pdf");
	c1->cd(1); h12->Draw();
	c1->cd(2); h22->Draw();
	c1->cd(3); h32->Draw();
	c1->cd(4); h42->Draw();
	c1->Print("check_hists.pdf");
	c1->cd(1); h13->Draw();
	c1->cd(2); h23->Draw();
	c1->cd(3); h33->Draw();
	c1->cd(4); h43->Draw();
	c1->cd(1); t3.DrawLatex(30,120,"Both #mu^{-} and #mu^{+} are counted");
	c1->Print("check_hists.pdf");
	c1->Clear();
	c1->Divide(2,2);
        c1->cd(1); h00->Draw();
	c1->cd(2); h01->Draw();
	c1->Print("check_hists.pdf)");

	TFile* f1 = new TFile("signal_hists.root", "RECREATE");
	h10->Write();
	h20->Write();
	h30->Write();
	h40->Write();
	h11->Write();
	h21->Write();
	h31->Write();
	h41->Write();
	h12->Write();
	h22->Write();
	h32->Write();
	h42->Write();
	h13->Write();
	h23->Write();
	h33->Write();
	h43->Write();
	f1->Close();
}
