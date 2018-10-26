#ifndef __CINT__
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParticle.h>
#endif

Float_t DeltaR(TParticle* part, TParticle* candidate) {
	Float_t DeltaPhi;
	Float_t DeltaPhip = part->Phi()-candidate->Phi();
	Float_t DeltaPhim = abs(part->Phi()-candidate->Phi()-TMath::TwoPi());
	Float_t DeltaEta = part->Eta()-candidate->Eta();
	if (DeltaPhip > DeltaPhim) DeltaPhi = DeltaPhim;
	else DeltaPhi = DeltaPhip;
	return sqrt(DeltaPhi*DeltaPhi+DeltaEta*DeltaEta);
}

Float_t Isolation(Int_t l, TParticle* part, TClonesArray* particles) {
	Float_t isol = 0;
	for (Long64_t k = 0; k < particles->GetEntries(); ++k) {
		TParticle* candidate = reinterpret_cast< TParticle* >(particles->At(k));
		Bool_t cut2 = candidate->Pt() > 0.1 && candidate->Eta() < 2.4 && candidate->Eta() > -2.4;
		if (cut2 && candidate->GetStatusCode() > 0 && k!=l && k!=part->GetDaughter(0) && k!=part->GetDaughter(1)) {
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

void generator_mass() {
	TStopwatch t;

	//Define all the histograms
	TCanvas* c2 = new TCanvas("c2", "trial analysis", 1);
	TH1D* h01b = new TH1D("h01b","m_{J/#psi#gamma} background; m_{J/#psi#gamma} (GeV/c^{2})",100, 0, 300);
	TH1D* h01s = new TH1D("h01s","m_{J/#psi#gamma} signal; m_{J/#psi#gamma} (GeV/c^{2})",100, 0, 300);
	THStack* hmass = new THStack("hmass","m_{J/#psi#gamma} signal and background; m_{J/#psi#gamma} (GeV/c^{2})");

	//Insert the .root file and get the tree and its Event branch
        TFile* f1 = TFile::Open("/data/cpapage/jpsiInclusive-gen100K.root");
        TTree* t1 = dynamic_cast< TTree* >(f1->Get("Pythia"));
        TClonesArray* particles1 = new TClonesArray("Event");
        t1->SetBranchAddress("Event", &particles1);
        for (Long64_t i = 0; i < t1->GetEntries(); ++i) { // This is the Event loop 
                t1->GetEntry(i);
                for (Long64_t l = 0; l < particles1->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles1->At(l));
			Bool_t cut1 = part->Pt() > 10 && part->Eta() < 2.4 && part->Eta() > -2.4;
			if (part->GetPdgCode() == 111 && cut1) {
				Float_t isolation = Isolation(l, part, particles1);
				if (isolation/part->Pt() < 0.5) {
					Int_t iJpsi = 0;
					for (Int_t k = 0; k < particles1->GetEntries(); ++k) {
						TParticle* temp = reinterpret_cast< TParticle* >(particles1->At(k));
						if (temp->GetPdgCode() == 443) {
							iJpsi = k;
							break;
						}
					}
					iJpsi = FindLastParticle(iJpsi, 443, particles1);
					Int_t iMuon = iJpsi;
					iMuon = FindLastParticle(iMuon, 13, particles1);
					Int_t iAntimuon = iJpsi;
					iAntimuon = FindLastParticle(iAntimuon, -13, particles1);

					TParticle* muon = reinterpret_cast< TParticle* >(particles1->At(iMuon));
					TParticle* antimuon = reinterpret_cast< TParticle* >(particles1->At(iAntimuon));
					if (1) {//cut_muon && cut_antimuon) {
						Float_t px = part->Px()+muon->Px()+antimuon->Px();
						Float_t py = part->Py()+muon->Py()+antimuon->Py();
						Float_t pz = part->Pz()+muon->Pz()+antimuon->Pz();
						Float_t E = part->Energy()+muon->Energy()+antimuon->Energy();
						h01b->Fill(sqrt(E*E-px*px-py*py-pz*pz));
					}
				}
			}
		}
	}

	//Insert the .root file and get the tree and its Event branch
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
		TParticle* muon = reinterpret_cast< TParticle* >(particles2->At(iMuon));
		TParticle* antimuon = reinterpret_cast< TParticle* >(particles2->At(iAntimuon));
		TParticle* gamma = reinterpret_cast< TParticle* >(particles2->At(iGamma));
		Bool_t cut1 = gamma->Eta() < 2.4 && gamma->Eta() > -2.4 && gamma->Pt() > 10;
		Bool_t cut2 = muon->Eta() < 2.4 && muon->Eta() > -2.4 && muon->Pt() > 10;
		Bool_t cut3 = antimuon->Eta() < 2.4 && antimuon->Eta() > -2.4 && antimuon->Pt() > 10;
		if (cut1 && cut2 && cut3) {
			Float_t px = gamma->Px()+muon->Px()+antimuon->Px();
			Float_t py = gamma->Py()+muon->Py()+antimuon->Py();
			Float_t pz = gamma->Pz()+muon->Pz()+antimuon->Pz();
			Float_t E = gamma->Energy()+muon->Energy()+antimuon->Energy();
			h01s->Fill(sqrt(E*E-px*px-py*py-pz*pz));
		}	
	}

	TLatex l1;
	TLatex l2; 

	h01s->SetLineColor(kRed);

	Int_t nEntries_sgn = h01s->GetEntries();
	Int_t nEntries_bkg = h01b->GetEntries();

	TH1D h01b_norm = 56600*(*h01b);
	TH1D h01s_norm = 3.7*0.00246*(*h01s);

	hmass->Add(&h01b_norm);
	hmass->Add(&h01s_norm);

	auto legend1 = new TLegend(0.75,0.7,0.9,0.9); legend1->AddEntry(&h01b_norm,"bkg","l"); legend1->AddEntry(&h01s_norm,"signal","l"); 
	legend1->SetTextSize(0.05);

	c2->Divide(2,2);
	c2->cd(1); h01b->Draw(); l1.DrawLatex(80,10, "bkg using nEvents = 100k"); 
	l1.DrawLatex(80,8, "Cuts applied only to fake #gamma ");
	l1.DrawLatex(80,7, "p_{T} > 10 GeV/c"); l1.DrawLatex(80,6, "and relative isolation < 0.5");
	c2->cd(2); h01s->Draw(); l2.DrawLatex(140,300, "signal using nEvents = 1k");
	l2.DrawLatex(140,250, "Cuts applied only to #gamma ");
	l2.DrawLatex(140,220, "p_{T} > 10 GeV/c");
        c2->Print("generator_mass_all_100K.pdf(");
	c2->Clear();
	hmass->Draw("hist"); legend1->Draw();
        c2->Print("generator_mass_all_100K.pdf)");

	TFile* f3 = new TFile("hmass.root", "RECREATE");
	h01b->Write();
	h01s->Write();
	h01b_norm.Write();
	h01s_norm.Write();
	hmass->Write();
	f3->Close();
	t.Stop();
	t.Print();
}
