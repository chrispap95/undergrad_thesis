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

void generator_pi0_jets() {
	TStopwatch t;

	//Define all the histograms
	TCanvas* c2 = new TCanvas("c2", "trial analysis", 1);
	TH1D* h1 = new TH1D("h1","#eta of isolated #pi^{0}, #eta",100, -6, 6);

	//Insert the .root file and get the tree and its Event branch
        TFile* f1 = TFile::Open("/data/cpapage/jpsiInclusive-gen100K.root");
        TTree* t1 = dynamic_cast< TTree* >(f1->Get("Pythia"));
        TClonesArray* particles1 = new TClonesArray("Event");
        t1->SetBranchAddress("Event", &particles1);
        for (Long64_t i = 0; i < t1->GetEntries(); ++i) { // This is the Event loop 
                t1->GetEntry(i);
                for (Long64_t l = 0; l < particles1->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles1->At(l));
			Bool_t cut1 = part->Pt() > 25;
			if (part->GetPdgCode() == 111 && cut1) {
				// Float_t isolation = Isolation(l, part, particles1);
				if (1) {//isolation/part->Pt() < 0.5) {
					h1->Fill(part->Eta());
				}
			}
		}
	}

	h1->Draw();

	t.Stop();
	t.Print();
}
