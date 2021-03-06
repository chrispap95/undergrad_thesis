#ifndef __CINT__
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParticle.h>
#endif

void generator_routine_2_bkg() {
	//Define all the histograms
	TCanvas* c1 = new TCanvas("c1", "trial analysis", 1);
	
	TH1F* h12b = new TH1F("h12b", "p_{T} of #gamma; p_{T} (GeV/c)", 100, 0, 100);
	TH1F* h22b = new TH1F("h22b", "y of #gamma; y", 100, -8, 8);
	TH1F* h32b = new TH1F("h32b", "#eta of #gamma; #eta", 100, -8, 8);
	TH1F* h42b = new TH1F("h42b", "E of #gamma; E(GeV)", 100, 0, 200);	
	TH1F* h13b = new TH1F("h13b", "p_{T} of #mu; p_{T} (GeV/c)", 100, 0, 60);
	TH1F* h23b = new TH1F("h23b", "y of #mu; y", 100, -8, 8);
	TH1F* h33b = new TH1F("h33b", "#eta of #mu; #eta", 100, -8, 8);
	TH1F* h43b = new TH1F("h43b", "E of #mu; E (GeV)", 100, 0, 200);
	TH1F* h53b = new TH1F("h53b", "p_{T} of #mu max; p_{T} (GeV/c)", 100, 0, 60);
	TH1F* h63b = new TH1F("h63b", "y of #mu max; y", 100, -8, 8);
	TH1F* h73b = new TH1F("h73b", "#eta of #mu max; #eta", 100, -8, 8);
	TH1F* h83b = new TH1F("h83b", "E of #mu max; E (GeV)", 100, 0, 200);
	TH1F* h93b = new TH1F("h93b", "p_{T} of #mu min; p_{T} (GeV/c)", 100, 0, 60);
	TH1F* h103b = new TH1F("h103b", "y of #mu min; y", 100, -8, 8);
	TH1F* h113b = new TH1F("h113b", "#eta of #mu min; #eta", 100, -8, 8);
	TH1F* h123b = new TH1F("h123b", "E of #mu min; E (GeV)", 100, 0, 200);
       	TH1F* h01b = new TH1F("h01b", "m_{J/#psi#gamma}; m_{J/#psi#gamma} (GeV/c^{2})", 100, 0, 250);
       	TH1F* h02b = new TH1F("h02b", "#Delta R_{J/#psi#gamma}; #Delta R", 100, 0, 10);
       	TH1F* h03b = new TH1F("h03b", "Isolation of #pi^{0} around #Delta R<0.3; Isolation (Gev/c)", 100, 0, 4000);
       	TH1F* h04b = new TH1F("h04b", "Relative Isolation of #pi^{0} around #Delta R<0.3; Relative Isolation", 100, 0, 300);

	TH1F* h12bc = new TH1F("h12bc", "p_{T} of #gamma; p_{T} (GeV/c)", 100, 0, 100);
	TH1F* h22bc = new TH1F("h22bc", "y of #gamma; y", 100, -8, 8);
	TH1F* h32bc = new TH1F("h32bc", "#eta of #gamma; #eta", 100, -8, 8);
	TH1F* h42bc = new TH1F("h42bc", "E of #gamma; E(GeV)", 100, 0, 200);	
	TH1F* h13bc = new TH1F("h13bc", "p_{T} of #mu; p_{T} (GeV/c)", 100, 0, 60);
	TH1F* h23bc = new TH1F("h23bc", "y of #mu; y", 100, -8, 8);
	TH1F* h33bc = new TH1F("h33bc", "#eta of #mu; #eta", 100, -8, 8);
	TH1F* h43bc = new TH1F("h43bc", "E of #mu; E (GeV)", 100, 0, 200);
	TH1F* h53bc = new TH1F("h53bc", "p_{T} of #mu max; p_{T} (GeV/c)", 100, 0, 60);
	TH1F* h63bc = new TH1F("h63bc", "y of #mu max; y", 100, -8, 8);
	TH1F* h73bc = new TH1F("h73bc", "#eta of #mu max; #eta", 100, -8, 8);
	TH1F* h83bc = new TH1F("h83bc", "E of #mu max; E (GeV)", 100, 0, 200);
	TH1F* h93bc = new TH1F("h93bc", "p_{T} of #mu min; p_{T} (GeV/c)", 100, 0, 60);
	TH1F* h103bc = new TH1F("h103bc", "y of #mu min; y", 100, -8, 8);
	TH1F* h113bc = new TH1F("h113bc", "#eta of #mu min; #eta", 100, -8, 8);
	TH1F* h123bc = new TH1F("h123bc", "E of #mu min; E (GeV)", 100, 0, 200);       
       	TH1F* h01bc = new TH1F("h01bc", "m_{J/#psi#gamma} accepted; m_{J/#psi#gamma} (GeV/c^{2})", 100, 0, 250);
       	TH1F* h02bc = new TH1F("h02bc", "#Delta R_{J/#psi#gamma} accepted; #Delta R", 100, 0, 10);

	TH1F* htemp = new TH1F("htemp","charge",100, -6, 6);
	TH1F* htemp1 = new TH1F("htemp1","pt",100, 0, 20);

	THStack* hs13 = new THStack("hs13", "p_{T} of #mu; p_{T} (GeV/c)");
	THStack* hs23 = new THStack("hs23", "y of #mu; y");
	THStack* hs33 = new THStack("hs33", "#eta of #mu; #eta");
	THStack* hs43 = new THStack("hs43", "E of #mu; E (GeV)");

	THStack* hs13c = new THStack("hs13c", "p_{T} of #mu; p_{T} (GeV/c)");
	THStack* hs23c = new THStack("hs23c", "y of #mu; y");
	THStack* hs33c = new THStack("hs33c", "#eta of #mu; #eta");
	THStack* hs43c = new THStack("hs43c", "E of #mu; E (GeV)");

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
			if (part->GetPdgCode() == 443) iJpsi = l;
			if (part->GetPdgCode() == 111 && part->Pt() >= 10) {
				Float_t isolation = 0;
				for (Long64_t k = 0; k < particles->GetEntries(); ++k) {
					TParticle* candidate = reinterpret_cast< TParticle* >(particles->At(k));
					Float_t DeltaPhi = part->Phi()-candidate->Phi();
					Float_t DeltaEta = part->Eta()-candidate->Eta();
					Float_t DeltaR = sqrt(DeltaPhi*DeltaPhi+DeltaEta*DeltaEta);
					TParticlePDG* pdg = reinterpret_cast< TParticlePDG* >(particles->At(k));
					if (DeltaR < 0.3 && k!=l) {
						isolation += candidate->Pt();
						htemp->Fill(pdg->Charge());
						htemp1->Fill(candidate->Pt());
					}
				}
				h03b->Fill(isolation);
				h04b->Fill(isolation/part->Pt());
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
		Float_t y_jpsi = jpsi->Y();
		Float_t eta_jpsi = jpsi->Eta();
		Float_t E_jpsi = jpsi->Energy();
		Float_t phi_jpsi = jpsi->Phi();
		Float_t pT_muon = muon->Pt();
		Float_t y_muon = muon->Y();
		Float_t eta_muon = muon->Eta();
		Float_t E_muon = muon->Energy();
		Float_t pT_antimuon = antimuon->Pt();
		Float_t y_antimuon = antimuon->Y();
		Float_t eta_antimuon = antimuon->Eta();
		Float_t E_antimuon = antimuon->Energy();

		Bool_t cut1 = (eta_muon < 2.4) && (eta_muon > -2.4) && (eta_antimuon < 2.4) && (eta_antimuon > -2.4);
		Bool_t cut2 = (pT_muon > 10) && (pT_antimuon > 10);

                for (Long64_t l = 0; l < particles->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles->At(l));
			if (part->GetPdgCode() == 22 && part->GetStatusCode() > 0 && part->Pt() > 2) {
				Float_t pT_gamma = part->Pt();
				Float_t y_gamma = part->Y();
				Float_t eta_gamma = part->Eta();
				Float_t E_gamma = part->Energy();
				Float_t phi_gamma = part->Phi();
				h12b->Fill(pT_gamma);
				h22b->Fill(y_gamma);
				h32b->Fill(eta_gamma);
				h42b->Fill(E_gamma);

				//DeltaR
				Float_t DeltaR = sqrt((phi_gamma-phi_jpsi)*(phi_gamma-phi_jpsi)+(eta_gamma-eta_jpsi)*(eta_gamma-eta_jpsi));
				h02b->Fill(DeltaR);
				if (part->Pt()>10 && part->Eta()<2.4 && part->Eta()>-2.4) {
					h12bc->Fill(pT_gamma);
					h22bc->Fill(y_gamma);
					h32bc->Fill(eta_gamma);
					h42bc->Fill(E_gamma);
					if (cut1 && cut2) h02bc->Fill(DeltaR);
				}

				//m_jpsigamma
				Float_t px_gamma = part->Px();
				Float_t py_gamma = part->Py();
				Float_t pz_gamma = part->Pz();
				Float_t px_jpsi = muon->Px()+antimuon->Px();
				Float_t py_jpsi = muon->Py()+antimuon->Py();
				Float_t pz_jpsi = muon->Pz()+antimuon->Pz();
				Float_t px_jpsigamma = px_gamma+px_jpsi;
				Float_t py_jpsigamma = py_gamma+py_jpsi;
				Float_t pz_jpsigamma = pz_gamma+pz_jpsi;
				Float_t p_jpsigamma = sqrt(px_jpsigamma*px_jpsigamma+py_jpsigamma*py_jpsigamma+pz_jpsigamma*pz_jpsigamma);
				h01b->Fill(sqrt((E_gamma+E_muon+E_antimuon)*(E_gamma+E_muon+E_antimuon)-p_jpsigamma*p_jpsigamma));
				if (part->Pt()>10 && part->Eta()<2.4 && part->Eta()>-2.4) {
					if (cut1 && cut2) h01bc->Fill(sqrt((E_gamma+E_muon+E_antimuon)*(E_gamma+E_muon+E_antimuon)-p_jpsigamma*p_jpsigamma));
				}
			}
                }

		h13b->Fill(pT_muon);
		h23b->Fill(y_muon);
		h33b->Fill(eta_muon);
		h43b->Fill(E_muon);
		h13b->Fill(pT_antimuon);
		h23b->Fill(y_antimuon);
		h33b->Fill(eta_antimuon);
		h43b->Fill(E_antimuon);
		if (pT_muon > pT_antimuon) { //Fill mu max & min histograms
			h53b->Fill(pT_muon);
			h63b->Fill(y_muon);
			h73b->Fill(eta_muon);
			h83b->Fill(E_muon);
			h93b->Fill(pT_antimuon);
			h103b->Fill(y_antimuon);
			h113b->Fill(eta_antimuon);
			h123b->Fill(E_antimuon);
		} else {
			h53b->Fill(pT_antimuon);
			h63b->Fill(y_antimuon);
			h73b->Fill(eta_antimuon);
			h83b->Fill(E_antimuon);
			h93b->Fill(pT_muon);
			h103b->Fill(y_muon);
			h113b->Fill(eta_muon);
			h123b->Fill(E_muon);
		}

		//Apply general cuts for mu
		if (cut1 && cut2) {
			h13bc->Fill(pT_muon);
			h23bc->Fill(y_muon);
			h33bc->Fill(eta_muon);
			h43bc->Fill(E_muon);
			h13bc->Fill(pT_antimuon);
			h23bc->Fill(y_antimuon);
			h33bc->Fill(eta_antimuon);
			h43bc->Fill(E_antimuon);
			if (pT_muon > pT_antimuon) { //Fill mu max & min histograms
				h53bc->Fill(pT_muon);
				h63bc->Fill(y_muon);
				h73bc->Fill(eta_muon);
				h83bc->Fill(E_muon);
				h93bc->Fill(pT_antimuon);
				h103bc->Fill(y_antimuon);
				h113bc->Fill(eta_antimuon);
				h123bc->Fill(E_antimuon);
			} else {
				h53bc->Fill(pT_antimuon);
				h63bc->Fill(y_antimuon);
				h73bc->Fill(eta_antimuon);
				h83bc->Fill(E_antimuon);
				h93bc->Fill(pT_muon);
				h103bc->Fill(y_muon);
				h113bc->Fill(eta_muon);
				h123bc->Fill(E_muon);
			}
		}
	}

	hs13->Add(h13b); hs13->Add(h53b); hs13->Add(h93b);
	hs23->Add(h23b); hs23->Add(h63b); hs23->Add(h103b);
	hs33->Add(h33b); hs33->Add(h73b); hs33->Add(h113b);
	hs43->Add(h43b); hs43->Add(h83b); hs43->Add(h123b);
	hs13c->Add(h13bc); hs13c->Add(h53bc); hs13c->Add(h93bc);
	hs23c->Add(h23bc); hs23c->Add(h63bc); hs23c->Add(h103bc);
	hs33c->Add(h33bc); hs33c->Add(h73bc); hs33c->Add(h113bc);
	hs43c->Add(h43bc); hs43c->Add(h83bc); hs43c->Add(h123bc);

	
	auto l13 = new TLegend(0.75,0.7,0.9,0.9); l13->AddEntry(h13b,"sum","l"); l13->AddEntry(h53b,"max","l"); 
	l13->AddEntry(h93b,"min","l"); l13->SetTextSize(0.05);
	auto l23 = new TLegend(0.75,0.7,0.9,0.9); l23->AddEntry(h23b,"sum","l"); l23->AddEntry(h63b,"max","l");  
	l23->AddEntry(h103b,"min","l"); l23->SetTextSize(0.05);
	auto l33 = new TLegend(0.75,0.7,0.9,0.9); l33->AddEntry(h33b,"sum","l"); l33->AddEntry(h73b,"max","l"); 
	l33->AddEntry(h113b,"min","l"); l33->SetTextSize(0.05);
	auto l43 = new TLegend(0.75,0.7,0.9,0.9); l43->AddEntry(h43b,"sum","l"); l43->AddEntry(h83b,"max","l"); 
	l43->AddEntry(h123b,"min","l"); l43->SetTextSize(0.05);
	auto l13c = new TLegend(0.75,0.7,0.9,0.9); l13c->AddEntry(h13bc,"sum","l"); l13c->AddEntry(h53bc,"max","l"); 
	l13c->AddEntry(h93bc,"min","l"); l13c->SetTextSize(0.05);
	auto l23c = new TLegend(0.75,0.7,0.9,0.9); l23c->AddEntry(h23bc,"sum","l"); l23c->AddEntry(h63bc,"max","l");  
	l23c->AddEntry(h103bc,"min","l"); l23c->SetTextSize(0.05);
	auto l33c = new TLegend(0.75,0.7,0.9,0.9); l33c->AddEntry(h33bc,"sum","l"); l33c->AddEntry(h73bc,"max","l"); 
	l33c->AddEntry(h113bc,"min","l"); l33c->SetTextSize(0.05);
	auto l43c = new TLegend(0.75,0.7,0.9,0.9); l43c->AddEntry(h43bc,"sum","l"); l43c->AddEntry(h83bc,"max","l"); 
	l43c->AddEntry(h123bc,"min","l"); l43c->SetTextSize(0.05);
	

	c1->Divide(2,2);
	c1->cd(1); h12b->Draw();
	c1->cd(2); h12bc->Draw();
	c1->cd(3); h22b->Draw();
	c1->cd(4); h22bc->Draw();
	c1->Print("generator_bkg_2.pdf(");
	c1->cd(1); h32b->Draw();
	c1->cd(2); h32bc->Draw();
	c1->cd(3); h42b->Draw();
	c1->cd(4); h42bc->Draw();
	c1->Print("generator_bkg_2.pdf");
	// c1->cd(1); h13b->Draw();
	// c1->cd(2); h13bc->Draw();
	// c1->cd(3); h23b->Draw();
	// c1->cd(4); h23bc->Draw();
        // c1->Print("generator_bkg_2.pdf");
	// c1->cd(1); h33b->Draw();
	// c1->cd(2); h33bc->Draw();
	// c1->cd(3); h43b->Draw();
	// c1->cd(4); h43bc->Draw();
        // c1->Print("generator_bkg_2.pdf");
	// c1->cd(1); h53b->Draw();
	// c1->cd(2); h53bc->Draw();
	// c1->cd(3); h63b->Draw();
	// c1->cd(4); h63bc->Draw();
        // c1->Print("generator_bkg_2.pdf");
	// c1->cd(1); h73b->Draw();
	// c1->cd(2); h73bc->Draw();
	// c1->cd(3); h83b->Draw();
	// c1->cd(4); h83bc->Draw();
        // c1->Print("generator_bkg_2.pdf");
	// c1->cd(1); h93b->Draw();
	// c1->cd(2); h93bc->Draw();
	// c1->cd(3); h103b->Draw();
	// c1->cd(4); h103bc->Draw();
        // c1->Print("generator_bkg_2.pdf");
	// c1->cd(1); h113b->Draw();
	// c1->cd(2); h113bc->Draw();
	// c1->cd(3); h123b->Draw();
	// c1->cd(4); h123bc->Draw();
        // c1->Print("generator_bkg_2.pdf");

	h53b->SetLineColor(kRed); h93b->SetLineColor(kGreen); 
	h63b->SetLineColor(kRed); h103b->SetLineColor(kGreen);
	h73b->SetLineColor(kRed); h113b->SetLineColor(kGreen);
	h83b->SetLineColor(kRed); h123b->SetLineColor(kGreen);
	h13b->SetStats(kFALSE); h23bc->SetStats(kFALSE); h33b->SetStats(kFALSE); h43b->SetStats(kFALSE);
	h53b->SetStats(kFALSE); h63bc->SetStats(kFALSE); h73b->SetStats(kFALSE); h83b->SetStats(kFALSE);
	h93b->SetStats(kFALSE); h103bc->SetStats(kFALSE); h113b->SetStats(kFALSE); h123b->SetStats(kFALSE);
	h53bc->SetLineColor(kRed); h93bc->SetLineColor(kGreen); 
	h63bc->SetLineColor(kRed); h103bc->SetLineColor(kGreen);
	h73bc->SetLineColor(kRed); h113bc->SetLineColor(kGreen);
	h83bc->SetLineColor(kRed); h123bc->SetLineColor(kGreen);
	h13bc->SetStats(kFALSE); h23bc->SetStats(kFALSE); h33bc->SetStats(kFALSE); h43bc->SetStats(kFALSE);
	h53bc->SetStats(kFALSE); h63bc->SetStats(kFALSE); h73bc->SetStats(kFALSE); h83bc->SetStats(kFALSE);
	h93bc->SetStats(kFALSE); h103bc->SetStats(kFALSE); h113bc->SetStats(kFALSE); h123bc->SetStats(kFALSE);

	c1->cd(1); hs13->Draw("nostack"); l13->Draw();
	c1->cd(2); hs13c->Draw("nostack"); l13c->Draw();
	c1->cd(3); hs23->Draw("nostack"); l23->Draw();
	c1->cd(4); hs23c->Draw("nostack"); l23c->Draw();
        c1->Print("generator_bkg_2.pdf");
	c1->cd(1); hs33->Draw("nostack"); l33->Draw();
	c1->cd(2); hs33c->Draw("nostack"); l33c->Draw();
	c1->cd(3); hs43->Draw("nostack"); l43->Draw();
	c1->cd(4); hs43c->Draw("nostack"); l43c->Draw();
        c1->Print("generator_bkg_2.pdf");
	c1->cd(1); h01b->Draw();
	c1->cd(2); h01bc->Draw();
	c1->cd(3); h02b->Draw();
	c1->cd(4); h02bc->Draw();
        c1->Print("generator_bkg_2.pdf");
	c1->Clear(); c1->Divide(2,2);
	c1->cd(1); h03b->Draw();
	c1->cd(2); h04b->Draw();
	c1->cd(3); htemp->Draw();
	c1->cd(4); htemp1->Draw();
        c1->Print("generator_bkg_2.pdf)");

	TFile* f1 = new TFile("generator_bkg_2.root", "RECREATE");
	h01b->Write();
	h02b->Write();
	h03b->Write();
	h04b->Write();
	h12b->Write();
	h22b->Write();
	h32b->Write();
	h42b->Write();
	h13b->Write();
	h23b->Write();
	h33b->Write();
	h43b->Write();
	h53b->Write();
	h63b->Write();
	h73b->Write();
	h83b->Write();
	h93b->Write();
	h103b->Write();
	h113b->Write();
	h123b->Write();
	h01bc->Write();
	h02bc->Write();
	h12bc->Write();
	h22bc->Write();
	h32bc->Write();
	h42bc->Write();
	h13bc->Write();
	h23bc->Write();
	h33bc->Write();
	h43bc->Write();
	h53bc->Write();
	h63bc->Write();
	h73bc->Write();
	h83bc->Write();
	h93bc->Write();
	h103bc->Write();
	h113bc->Write();
	h123bc->Write();
	f1->Close();
}
