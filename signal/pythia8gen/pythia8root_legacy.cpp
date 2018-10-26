#include "Pythia.h" 
#include "TPythia8.h" 
#include "TParticle.h" 
#include "TClonesArray.h" 
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TError.h"
#include "TLatex.h"
#include "TText.h"

/*This version is up to date as of 
19:00 - 29/4/2018 */

using namespace Pythia8; 	// Let Pythia8:: be implicit.

int main(int argc, char* argv[]) {	// Begin main program.

	if (argc < 2) {
		cerr << "Missing Pythia data card. " << endl;
		return 1;
	}

	if (argc < 3) {
		cerr << "Missing output filename. " << endl;
		return 1;
	}


	const char *p8dataenv = gSystem->Getenv("PYTHIA8DATA");
	if (!p8dataenv) {
		const char *p8env = gSystem->Getenv("PYTHIA8");
		if (!p8env) {
			Error("pythia8.C", "Environment variable PYTHIA8 must contain path to pythia directory!");
			return 1;
		}
		TString p8d = p8env;
		p8d += "/xmldoc";
		gSystem->Setenv("PYTHIA8DATA", p8d);
	}

	const char* path = gSystem->ExpandPathName("$PYTHIA8DATA");
	if (gSystem->AccessPathName(path)) {
		Error("pythia8.C", "Environment variable PYTHIA8DATA must contain path to $PYTHIA8/xmldoc directory !");
		return 1;
	}

	// Set up generation.
	TPythia8* pythia = new TPythia8();
	pythia->ReadConfigFile(argv[1]);	// Read Pythia configuration from data card
	pythia->Initialize(2212, 2212, 13000.); // Initialize Pythia
	Pythia8::Pythia* p8 = pythia->Pythia8();
	int nEvent = p8->mode("Main:numberOfEvents");    // set number of events to generate
	
	// Setup ROOT output
	TFile* f = TFile::Open(argv[2], "RECREATE");
	TClonesArray* particles = new TClonesArray("TParticle", 1000);
	TClonesArray& pref = *particles;
	TTree* t = new TTree("Pythia", "Pythia");
	t->Branch("Event", "TClonesArray", &pref, 256000, 0);
	particles->BypassStreamer();		

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

	// Generate event(s). (aka "Event loop")
	int prog, oldprog = 0;	// progress indicator variables
	for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
		// print out progress 
		prog = (int)(iEvent * 100. / nEvent);
		if (prog > oldprog) {
			oldprog = prog;
			cout << prog << "\% complete..." << endl;
		}
		pythia->GenerateEvent();
		pref.Clear();
		pythia->ImportParticles(particles,"All");
		Int_t np = particles->GetEntriesFast();
		t->Fill();
		
		Int_t iHiggs = 0;
		Int_t iJpsi = 0;
		Int_t iGamma = 0;
		Int_t iMuon = 0;
		Int_t iAntiMuon = 0;
		for (Int_t ip = 0; ip < np; ip++) {    //Find the first Higgs
			TParticle* part = (TParticle*)particles->At(ip);
			if (part->GetPdgCode() == 25) {
				iHiggs = ip;
				break;
			}
		}
		while (1) {    //Find the last Higgs
			TParticle* part = (TParticle*)particles->At(iHiggs);
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = (TParticle*)particles->At(daughter0id);
				if (daughter0->GetPdgCode() == 25) {
					iHiggs = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = (TParticle*)particles->At(daughter1id);
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
			TParticle* part = (TParticle*)particles->At(iJpsi);
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = (TParticle*)particles->At(daughter0id);
				if (daughter0->GetPdgCode() == 443) {
					iJpsi = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = (TParticle*)particles->At(daughter1id);
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
			TParticle* part = (TParticle*)particles->At(iGamma);
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = (TParticle*)particles->At(daughter0id);
				if (daughter0->GetPdgCode() == 22) {
					iGamma = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = (TParticle*)particles->At(daughter1id);
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
			TParticle* part = (TParticle*)particles->At(iMuon);
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = (TParticle*)particles->At(daughter0id);
				if (daughter0->GetPdgCode() == 13) {
					iMuon = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = (TParticle*)particles->At(daughter1id);
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
			TParticle* part = (TParticle*)particles->At(iAntiMuon);
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = (TParticle*)particles->At(daughter0id);
				if (daughter0->GetPdgCode() == -13) {
					iAntiMuon = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = (TParticle*)particles->At(daughter1id);
				if (daughter1->GetPdgCode() == -13) {
					iAntiMuon = part->GetDaughter(1);
					daughter1validity = false;
					continue;
				}
			}
			if (daughter0validity && daughter1validity) break;
		}

		TParticle* higgs = (TParticle*)particles->At(iHiggs);
		TParticle* jpsi = (TParticle*)particles->At(iJpsi);
		TParticle* gamma = (TParticle*)particles->At(iGamma);
		TParticle* muon = (TParticle*)particles->At(iMuon);
		TParticle* antimuon = (TParticle*)particles->At(iAntiMuon);

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


	pythia->PrintStatistics();		 // print out pythia statistics

	// Save ROOT output to file
	t->Print();
	t->Write();
	f->Close();	
	delete f;
	
	return 0;
} // End main program with error-free return.
