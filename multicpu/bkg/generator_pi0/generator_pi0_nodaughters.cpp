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

Float_t Isolation(Float_t R, Int_t l, TParticle* part, TClonesArray* particles) {
	Float_t isol = 0;
	for (Long64_t k = 0; k < particles->GetEntries(); ++k) {
		TParticle* candidate = reinterpret_cast< TParticle* >(particles->At(k));
		Bool_t cut2 = 1;//candidate->Eta() < 2.37 && candidate->Eta() > -2.37;
		if (cut2 && candidate->GetStatusCode() > 0 && k != l && k != part->GetDaughter(0) && k != part->GetDaughter(1)) {
			Float_t deltaR = DeltaR(part, candidate);
			TParticlePDG* candpdg = candidate->GetPDG(0);
			if (deltaR < R && candpdg) {
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

void generator_pi0_nodaughters() {
	TStopwatch t;

	//Define all the histograms
	TCanvas* c2 = new TCanvas("c2", "trial analysis", 1);
	TH1F* h03b = new TH1F("h03b","isolation p^{cone0.4}_{T}+E^{cone0.4}_{T} for bkg;GeV",100, 0, 10);
	TH1F* h04b = new TH1F("h04b","relative isolation (p^{cone0.4}_{T}+E^{cone0.4}_{T})/p_{T} for bkg",100, 0, 2);
	TH1F* h03s = new TH1F("h03s","isolation p^{cone0.4}_{T}+E^{cone0.4}_{T} for signal;GeV",100, 0, 10);
	TH1F* h04s = new TH1F("h04s","relative isolation (p^{cone0.4}_{T}+E^{cone0.4}_{T})/p_{T} for signal",100, 0, 2);

	//Insert the .root file and get the tree and its Event branch
        TFile* f1 = TFile::Open("/data/cpapage/jpsiInclusive-gen8K.root");
        TTree* t1 = dynamic_cast< TTree* >(f1->Get("Pythia"));
        TClonesArray* particles1 = new TClonesArray("Event");
        t1->SetBranchAddress("Event", &particles1);
        for (Long64_t i = 0; i < t1->GetEntries(); ++i) { // This is like the Event loop 
                t1->GetEntry(i);
                for (Long64_t l = 0; l < particles1->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles1->At(l));
			Bool_t cut1 = part->Eta() < 2.37 && part->Eta() > -2.37 && part->Pt() > 10;
			if (part->GetPdgCode() == 111) {
			        h03b->Fill(Isolation(0.2,l,part,particles1));
				h04b->Fill(Isolation(0.2,l,part,particles1)/part->Pt());
			}
		}
	}

	//Insert the .root file and get the tree and its Event branch
        TFile* f2 = TFile::Open("/data/cpapage/higgsTojpsiGamma-gen8K.root");
        TTree* t2 = dynamic_cast< TTree* >(f2->Get("Pythia"));
        TClonesArray* particles2 = new TClonesArray("Event");
        t2->SetBranchAddress("Event", &particles2);
        for (Long64_t i = 0; i < t2->GetEntries(); ++i) { // This is like the Event loop 
                t2->GetEntry(i);

		Int_t iHiggs = 0;
		Int_t iGamma = 0;
                for (Long64_t l = 0; l < particles2->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles2->At(l));
			if (part->GetPdgCode() == 25) {
				iHiggs = l;
				break;
			}
		}
		iHiggs = FindLastParticle(iHiggs,25,particles2);
		iGamma = FindLastParticle(iHiggs,22,particles2);

		TParticle* part = reinterpret_cast< TParticle* >(particles2->At(iGamma));
		if (part->Pt()>10 && part->Eta() > -2.37 && part->Eta() < 2.37) {
			h03s->Fill(Isolation(0.2,iGamma,part,particles2));
			h04s->Fill(Isolation(0.2,iGamma,part,particles2)/part->Pt());
		}
	}

	c2->Divide(2,2);
	c2->cd(1); h03b->Draw();
	c2->cd(2); h04b->Draw();
	c2->cd(3); h03s->Draw();
	c2->cd(4); h04s->Draw();

	TFile* fnew = new TFile("hisol.root","RECREATE");
	h03b->Write();
	h04b->Write();
	h03s->Write();
	h04s->Write();
	fnew->Close();

	t.Stop();
	t.Print();
}
