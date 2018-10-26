#ifndef __CINT__
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParticle.h>
#endif

void generator_routine_1_bkg() {
	//Define all the histograms
	TCanvas* c1 = new TCanvas("c1", "trial analysis", 1);
	TH1F* h11b = new TH1F("h11b", "p_{T} of J/#psi; p_{T} (GeV/c)", 100, 0, 200);
	TH1F* h21b = new TH1F("h21b", "y of J/#psi; y", 100, -8, 8);
	TH1F* h31b = new TH1F("h31b", "#eta of J/#psi; #eta", 100, -8, 8);
	TH1F* h41b = new TH1F("h41b", "E of J/#psi; E(GeV)", 100, 0, 1000);
	TH1F* h12b = new TH1F("h12b", "p_{T} of #gamma; p_{T} (GeV/c)", 100, 0, 200);
	TH1F* h22b = new TH1F("h22b", "y of #gamma; y", 100, -11, 11);
	TH1F* h32b = new TH1F("h32b", "#eta of #gamma; #eta", 100, -11, 11);
	TH1F* h42b = new TH1F("h42b", "E of #gamma; E(GeV)", 100, 0, 1000);
	TH1F* h13b = new TH1F("h13b", "p_{T} of #mu; p_{T} (GeV/c)", 100, 0, 150);
	TH1F* h23b = new TH1F("h23b", "y of #mu; y", 100, -8, 8);
	TH1F* h33b = new TH1F("h33b", "#eta of #mu; #eta", 100, -8, 8);
	TH1F* h43b = new TH1F("h43b", "E of #mu; E (GeV)", 100, 0, 600);
	
	//Insert the .root file and get the tree and its Event branch
        TFile* f = TFile::Open("/data/cpapage/jpsiInclusive-gen8K.root");
        TTree* t = dynamic_cast< TTree* >(f->Get("Pythia"));
        TClonesArray* particles = new TClonesArray("Event");
        t->SetBranchAddress("Event", &particles);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) { // This is like the Event loop 
                t->GetEntry(i);

		Int_t iJpsi = 0;
		Int_t iMuon = 0;
		Int_t iAntiMuon = 0;

                for (Long64_t l = 0; l < particles->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles->At(l));
			if (part->GetPdgCode() == 443 && iJpsi == 0) iJpsi = l;
			if (part->GetPdgCode() == 22 && part->GetStatusCode() > 0 && part->Pt()>4) {
				Float_t pT_gamma = part->Pt();
				h12b->Fill(pT_gamma);
				Float_t y_gamma = part->Y();
				h22b->Fill(y_gamma);
				Float_t eta_gamma = part->Eta();
				h32b->Fill(eta_gamma);
				Float_t E_gamma = part->Energy();
				h42b->Fill(E_gamma);
			}
                }
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

		TParticle* jpsi = reinterpret_cast< TParticle* >(particles->At(iJpsi));
		TParticle* muon = reinterpret_cast< TParticle* >(particles->At(iMuon));
		TParticle* antimuon = reinterpret_cast< TParticle* >(particles->At(iAntiMuon));

		Float_t pT_jpsi = jpsi->Pt();
		h11b->Fill(pT_jpsi);
		Float_t y_jpsi = jpsi->Y();
		h21b->Fill(y_jpsi);
		Float_t eta_jpsi = jpsi->Eta();
		h31b->Fill(eta_jpsi);
		Float_t E_jpsi = jpsi->Energy();
		h41b->Fill(E_jpsi);

		Float_t pT_muon = muon->Pt();
		h13b->Fill(pT_muon);
		Float_t y_muon = muon->Y();
		h23b->Fill(y_muon);
		Float_t eta_muon = muon->Eta();
		h33b->Fill(eta_muon);
		Float_t E_muon = muon->Energy();
		h43b->Fill(E_muon);

		Float_t pT_antimuon = antimuon->Pt();
		h13b->Fill(pT_antimuon);
		Float_t y_antimuon = antimuon->Y();
		h23b->Fill(y_antimuon);
		Float_t eta_antimuon = antimuon->Eta();
		h33b->Fill(eta_antimuon);
		Float_t E_antimuon = antimuon->Energy();
		h43b->Fill(E_antimuon);
        }

	TText* t1 = new TText(30,500,"nEvents = 8000");
	TText* t4 = new TText(30,450,"J/psi inclusive production");
	TText* t2 = new TText(30,400,"Pythia generator results");
	TLatex t3;
	
	c1->Clear();
	c1->Divide(2,2);
	c1->cd(1); h11b->Draw(); t1->Draw(); t2->Draw(); t4->Draw();
	c1->cd(2); h21b->Draw();
	c1->cd(3); h31b->Draw();
	c1->cd(4); h41b->Draw();
	c1->Print("jpsi_inclusive.pdf(");
	c1->cd(1); h12b->Draw();
	c1->cd(2); h22b->Draw();
	c1->cd(3); h32b->Draw();
	c1->cd(4); h42b->Draw();
	c1->Print("jpsi_inclusive.pdf");
	c1->cd(1); h13b->Draw(); t3.DrawLatex(15,800,"Both #mu^{-} and #mu^{+} are counted");
	c1->cd(2); h23b->Draw();
	c1->cd(3); h33b->Draw();
	c1->cd(4); h43b->Draw();
	c1->Print("jpsi_inclusive.pdf)");

	TFile* f1 = new TFile("bkg_hists.root", "RECREATE");
	h11b->Write();
	h21b->Write();
	h31b->Write();
	h41b->Write();
	h12b->Write();
	h22b->Write();
	h32b->Write();
	h42b->Write();
	h13b->Write();
	h23b->Write();
	h33b->Write();
	h43b->Write();
	f1->Close();
}
