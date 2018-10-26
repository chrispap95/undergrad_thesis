#ifndef __CINT__
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParticle.h>
#endif


void generator_routine_2() {
	//Define all the histograms
	TCanvas* c1 = new TCanvas("c1", "trial analysis", 1);
	
	TH1F* h12 = new TH1F("h12", "p_{T} of #gamma; p_{T} (GeV/c)", 100, 0, 200);
	TH1F* h22 = new TH1F("h22", "y of #gamma; y", 100, -11, 11);
	TH1F* h32 = new TH1F("h32", "#eta of #gamma; #eta", 100, -11, 11);
	TH1F* h42 = new TH1F("h42", "E of #gamma; E(GeV)", 100, 0, 1000);
       	TH1F* h13 = new TH1F("h13", "p_{T} of #mu; p_{T} (GeV/c)", 100, 0, 150);
	TH1F* h23 = new TH1F("h23", "y of #mu; y", 100, -8, 8);
	TH1F* h33 = new TH1F("h33", "#eta of #mu; #eta", 100, -8, 8);
	TH1F* h43 = new TH1F("h43", "E of #mu; E (GeV)", 100, 0, 600);
	TH1F* h53 = new TH1F("h53", "p_{T} of #mu max; p_{T} (GeV/c)", 100, 0, 150);
	TH1F* h63 = new TH1F("h63", "y of #mu max; y", 100, -8, 8);
	TH1F* h73 = new TH1F("h73", "#eta of #mu max; #eta", 100, -8, 8);
	TH1F* h83 = new TH1F("h83", "E of #mu max; E (GeV)", 100, 0, 600);
	TH1F* h93 = new TH1F("h93", "p_{T} of #mu min; p_{T} (GeV/c)", 100, 0, 150);
	TH1F* h103 = new TH1F("h103", "y of #mu min; y", 100, -8, 8);
	TH1F* h113 = new TH1F("h113", "#eta of #mu min; #eta", 100, -8, 8);
	TH1F* h123 = new TH1F("h123", "E of #mu min; E (GeV)", 100, 0, 600);
       	TH1F* h01 = new TH1F("h01", "m_{J/#psi#gamma}; m_{J/#psi#gamma} (GeV/c^{2})", 100, 0, 250);
       	TH1F* h02 = new TH1F("h02", "#Delta R_{J/#psi#gamma}; #Delta R", 100, 0, 10);

	TH1F* h12c = new TH1F("h12c", "p_{T} of #gamma accepted; p_{T} (GeV/c)", 100, 0, 200);
	TH1F* h22c = new TH1F("h22c", "y of #gamma accepted; y", 100, -11, 11);
	TH1F* h32c = new TH1F("h32c", "#eta of #gamma accepted; #eta", 100, -11, 11);
	TH1F* h42c = new TH1F("h42c", "E of #gamma accepted; E(GeV)", 100, 0, 1000);
       	TH1F* h13c = new TH1F("h13c", "p_{T} of #mu accepted; p_{T} (GeV/c)", 100, 0, 150);
	TH1F* h23c = new TH1F("h23c", "y of #mu accepted; y", 100, -8, 8);
	TH1F* h33c = new TH1F("h33c", "#eta of #mu accepted; #eta", 100, -8, 8);
	TH1F* h43c = new TH1F("h43c", "E of #mu accepted; E (GeV)", 100, 0, 600);
	TH1F* h53c = new TH1F("h53c", "p_{T} of #mu max accepted; p_{T} (GeV/c)", 100, 0, 150);
	TH1F* h63c = new TH1F("h63c", "y of #mu max accepted; y", 100, -8, 8);
	TH1F* h73c = new TH1F("h73c", "#eta of #mu max accepted; #eta", 100, -8, 8);
	TH1F* h83c = new TH1F("h83c", "E of #mu max accepted; E (GeV)", 100, 0, 600);
	TH1F* h93c = new TH1F("h93c", "p_{T} of #mu min accepted; p_{T} (GeV/c)", 100, 0, 150);
	TH1F* h103c = new TH1F("h103c", "y of #mu min accepted; y", 100, -8, 8);
	TH1F* h113c = new TH1F("h113c", "#eta of #mu min accepted; #eta", 100, -8, 8);
	TH1F* h123c = new TH1F("h123c", "E of #mu min accepted; E (GeV)", 100, 0, 600);
       	TH1F* h01c = new TH1F("h01c", "m_{J/#psi#gamma} accepted; m_{J/#psi#gamma} (GeV/c^{2})", 100, 0, 250);
       	TH1F* h02c = new TH1F("h02c", "#Delta R_{J/#psi#gamma} accepted; #Delta R", 100, 0, 10);
	
	THStack* hs13 = new THStack("hs13", "p_{T} of #mu; p_{T} (GeV/c)");
	THStack* hs23 = new THStack("hs23", "y of #mu; y");
	THStack* hs33 = new THStack("hs33", "#eta of #mu; #eta");
	THStack* hs43 = new THStack("hs43", "E of #mu; E (GeV)");

	THStack* hs13c = new THStack("hs13c", "p_{T} of #mu accepted; p_{T} (GeV/c)");
	THStack* hs23c = new THStack("hs23c", "y of #mu accepted; y");
	THStack* hs33c = new THStack("hs33c", "#eta of #mu accepted; #eta");
	THStack* hs43c = new THStack("hs43c", "E of #mu accepted; E (GeV)");

	//Insert the .root file and get the tree and its Event branch
        TFile* f = TFile::Open("/data/cpapage/higgsTojpsiGamma-gen8K.root");
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
			if (part->GetPdgCode() == 25 && iHiggs == 0) {
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

		TParticle* jpsi = reinterpret_cast< TParticle* >(particles->At(iJpsi));
		Float_t pT_jpsi = jpsi->Pt(); Float_t y_jpsi = jpsi->Y(); Float_t phi_jpsi = jpsi->Phi();
		Float_t eta_jpsi = jpsi->Eta(); Float_t E_jpsi = jpsi->Energy();
		TParticle* gamma = reinterpret_cast< TParticle* >(particles->At(iGamma));
		Float_t pT_gamma = gamma->Pt();	Float_t y_gamma = gamma->Y(); Float_t phi_gamma = gamma->Phi();
		Float_t eta_gamma = gamma->Eta(); Float_t E_gamma = gamma->Energy();
		TParticle* muon = reinterpret_cast< TParticle* >(particles->At(iMuon));
		Float_t pT_muon = muon->Pt(); Float_t y_muon = muon->Y();
		Float_t eta_muon = muon->Eta();	Float_t E_muon = muon->Energy();
		TParticle* antimuon = reinterpret_cast< TParticle* >(particles->At(iAntiMuon));
		Float_t pT_antimuon = antimuon->Pt(); Float_t y_antimuon = antimuon->Y();
		Float_t eta_antimuon = antimuon->Eta(); Float_t E_antimuon = antimuon->Energy();

		//Fill inclusive histograms
		h12->Fill(pT_gamma);
		h22->Fill(y_gamma);
		h32->Fill(eta_gamma);
		h42->Fill(E_gamma);
		h13->Fill(pT_muon);
		h23->Fill(y_muon);
		h33->Fill(eta_muon);
		h43->Fill(E_muon);
		h13->Fill(pT_antimuon);
		h23->Fill(y_antimuon);
		h33->Fill(eta_antimuon);
		h43->Fill(E_antimuon);
		if (pT_muon > pT_antimuon) { //Fill mu max & min histograms
			h53->Fill(pT_muon);
			h63->Fill(y_muon);
			h73->Fill(eta_muon);
			h83->Fill(E_muon);
			h93->Fill(pT_antimuon);
			h103->Fill(y_antimuon);
			h113->Fill(eta_antimuon);
			h123->Fill(E_antimuon);
		} else {
			h53->Fill(pT_antimuon);
			h63->Fill(y_antimuon);
			h73->Fill(eta_antimuon);
			h83->Fill(E_antimuon);
			h93->Fill(pT_muon);
			h103->Fill(y_muon);
			h113->Fill(eta_muon);
			h123->Fill(E_muon);
		}

		//DeltaR
		Float_t DeltaR = sqrt((phi_gamma-phi_jpsi)*(phi_gamma-phi_jpsi)+(eta_gamma-eta_jpsi)*(eta_gamma-eta_jpsi));
		h02->Fill(DeltaR);

		//Apply general cuts for mu and gamma
		Bool_t cut1 = (eta_gamma < 2.4) && (eta_gamma > -2.4);
		Bool_t cut2 = (eta_muon < 2.4) && (eta_muon > -2.4) && (eta_antimuon < 2.4) && (eta_antimuon > -2.4);
		Bool_t cut3 = (pT_gamma > 10) && (pT_muon > 10) && (pT_antimuon > 10);

		//m_jpsigamma
		Float_t px_gamma = gamma->Px();
		Float_t py_gamma = gamma->Py();
		Float_t pz_gamma = gamma->Pz();
		Float_t px_jpsi = jpsi->Px();
		Float_t py_jpsi = jpsi->Py();
		Float_t pz_jpsi = jpsi->Pz();
		Float_t px_jpsigamma = px_gamma+px_jpsi;
		Float_t py_jpsigamma = py_gamma+py_jpsi;
		Float_t pz_jpsigamma = pz_gamma+pz_jpsi;
		Float_t p_jpsigamma = sqrt(px_jpsigamma*px_jpsigamma+py_jpsigamma*py_jpsigamma+pz_jpsigamma*pz_jpsigamma);
		h01->Fill(sqrt((E_gamma+E_jpsi)*(E_gamma+E_jpsi)-p_jpsigamma*p_jpsigamma));

		//Fill the acceptance histograms
		if (cut1 && cut2 && cut3) {
			h12c->Fill(pT_gamma);
			h22c->Fill(y_gamma);
			h32c->Fill(eta_gamma);
			h42c->Fill(E_gamma);
			h13c->Fill(pT_muon);
			h23c->Fill(y_muon);
			h33c->Fill(eta_muon);
			h43c->Fill(E_muon);
			h13c->Fill(pT_antimuon);
			h23c->Fill(y_antimuon);
			h33c->Fill(eta_antimuon);
			h43c->Fill(E_antimuon);
			h01c->Fill(sqrt((E_gamma+E_jpsi)*(E_gamma+E_jpsi)-p_jpsigamma*p_jpsigamma));
			h02c->Fill(DeltaR);
			if (pT_muon > pT_antimuon) { //Fill mu max & min histograms
				h53c->Fill(pT_muon);
				h63c->Fill(y_muon);
				h73c->Fill(eta_muon);
				h83c->Fill(E_muon);
				h93c->Fill(pT_antimuon);
				h103c->Fill(y_antimuon);
				h113c->Fill(eta_antimuon);
				h123c->Fill(E_antimuon);
			} else {
				h53c->Fill(pT_antimuon);
				h63c->Fill(y_antimuon);
				h73c->Fill(eta_antimuon);
				h83c->Fill(E_antimuon);
				h93c->Fill(pT_muon);
				h103c->Fill(y_muon);
				h113c->Fill(eta_muon);
				h123c->Fill(E_muon);
			}
		}
	}

	hs13->Add(h13); hs13->Add(h53); hs13->Add(h93);
	hs23->Add(h23); hs23->Add(h63); hs23->Add(h103);
	hs33->Add(h33); hs33->Add(h73); hs33->Add(h113);
	hs43->Add(h43); hs43->Add(h83); hs43->Add(h123);
	
	auto l13 = new TLegend(0.75,0.7,0.9,0.9); l13->AddEntry(h13,"sum","l"); l13->AddEntry(h53,"max","l"); 
	l13->AddEntry(h93,"min","l"); l13->SetTextSize(0.05);
	auto l23 = new TLegend(0.75,0.7,0.9,0.9); l23->AddEntry(h23,"sum","l"); l23->AddEntry(h63,"max","l");  
	l23->AddEntry(h103,"min","l"); l23->SetTextSize(0.05);
	auto l33 = new TLegend(0.75,0.7,0.9,0.9); l33->AddEntry(h33,"sum","l"); l33->AddEntry(h73,"max","l"); 
	l33->AddEntry(h113,"min","l"); l33->SetTextSize(0.05);
	auto l43 = new TLegend(0.75,0.7,0.9,0.9); l43->AddEntry(h43,"sum","l"); l43->AddEntry(h83,"max","l"); 
	l43->AddEntry(h123,"min","l"); l43->SetTextSize(0.05);

	hs13c->Add(h13c); hs13c->Add(h53c); hs13c->Add(h93c);
	hs23c->Add(h23c); hs23c->Add(h63c); hs23c->Add(h103c);
	hs33c->Add(h33c); hs33c->Add(h73c); hs33c->Add(h113c);
	hs43c->Add(h43c); hs43c->Add(h83c); hs43c->Add(h123c);
	
	auto l13c = new TLegend(0.75,0.7,0.9,0.9); l13c->AddEntry(h13c,"sum","l"); l13c->AddEntry(h53c,"max","l"); 
	l13c->AddEntry(h93c,"min","l"); l13c->SetTextSize(0.05);
	auto l23c = new TLegend(0.75,0.7,0.9,0.9); l23c->AddEntry(h23c,"sum","l"); l23c->AddEntry(h63c,"max","l");  
	l23c->AddEntry(h103c,"min","l"); l23c->SetTextSize(0.05);
	auto l33c = new TLegend(0.75,0.7,0.9,0.9); l33c->AddEntry(h33c,"sum","l"); l33c->AddEntry(h73c,"max","l"); 
	l33c->AddEntry(h113c,"min","l"); l33c->SetTextSize(0.05);
	auto l43c = new TLegend(0.75,0.7,0.9,0.9); l43c->AddEntry(h43c,"sum","l"); l43c->AddEntry(h83c,"max","l"); 
	l43c->AddEntry(h123c,"min","l"); l43c->SetTextSize(0.05);	

	c1->Divide(2,2);
	c1->cd(1); h12->Draw();
	c1->cd(2); h12c->Draw();
	c1->cd(3); h22->Draw();
	c1->cd(4); h22c->Draw();
	c1->Print("generator_signal_2.pdf(");
	c1->cd(1); h32->Draw();
	c1->cd(2); h32c->Draw();
	c1->cd(3); h42->Draw();
	c1->cd(4); h42c->Draw();
        c1->Print("generator_signal_2.pdf");

	h53->SetLineColor(kRed); h93->SetLineColor(kGreen); 
	h63->SetLineColor(kRed); h103->SetLineColor(kGreen);
	h73->SetLineColor(kRed); h113->SetLineColor(kGreen);
	h83->SetLineColor(kRed); h123->SetLineColor(kGreen);
	h13->SetStats(kFALSE); h23->SetStats(kFALSE); h33->SetStats(kFALSE); h43->SetStats(kFALSE);
	h53->SetStats(kFALSE); h63->SetStats(kFALSE); h73->SetStats(kFALSE); h83->SetStats(kFALSE);
	h93->SetStats(kFALSE); h103->SetStats(kFALSE); h113->SetStats(kFALSE); h123->SetStats(kFALSE);

	h53c->SetLineColor(kRed); h93c->SetLineColor(kGreen); 
	h63c->SetLineColor(kRed); h103c->SetLineColor(kGreen);
	h73c->SetLineColor(kRed); h113c->SetLineColor(kGreen);
	h83c->SetLineColor(kRed); h123c->SetLineColor(kGreen);
	h13c->SetStats(kFALSE); h23c->SetStats(kFALSE); h33c->SetStats(kFALSE); h43c->SetStats(kFALSE);
	h53c->SetStats(kFALSE); h63c->SetStats(kFALSE); h73c->SetStats(kFALSE); h83c->SetStats(kFALSE);
	h93c->SetStats(kFALSE); h103c->SetStats(kFALSE); h113c->SetStats(kFALSE); h123c->SetStats(kFALSE);

	c1->cd(1); hs13->Draw("nostack"); l13->Draw();
	c1->cd(2); hs13c->Draw("nostack"); l13c->Draw();
	c1->cd(3); hs23->Draw("nostack"); l23->Draw();
	c1->cd(4); hs23c->Draw("nostack"); l23c->Draw();
        c1->Print("generator_signal_2.pdf");
	c1->cd(1); hs33->Draw("nostack"); l33->Draw();
	c1->cd(2); hs33c->Draw("nostack"); l33c->Draw();
	c1->cd(3); hs43->Draw("nostack"); l43->Draw();
	c1->cd(4); hs43c->Draw("nostack"); l43c->Draw();
        c1->Print("generator_signal_2.pdf");
       	c1->cd(1); h01->Draw();
	c1->cd(2); h01c->Draw();
	c1->cd(3); h02->Draw();
	c1->cd(4); h02c->Draw();
        c1->Print("generator_signal_2.pdf)");


	TFile* f1 = new TFile("generator_signal_2.root", "RECREATE");
	h01->Write();
	h02->Write();
	h12->Write();
	h22->Write();
	h32->Write();
	h42->Write();
	h13->Write();
	h23->Write();
	h33->Write();
	h43->Write();
	h53->Write();
	h63->Write();
	h73->Write();
	h83->Write();
	h93->Write();
	h103->Write();
	h113->Write();
	h123->Write();
	h01c->Write();
	h02c->Write();
	h12c->Write();
	h22c->Write();
	h32c->Write();
	h42c->Write();
	h13c->Write();
	h23c->Write();
	h33c->Write();
	h43c->Write();
	h53c->Write();
	h63c->Write();
	h73c->Write();
	h83c->Write();
	h93c->Write();
	h103c->Write();
	h113c->Write();
	h123c->Write();
	f1->Close();
}
