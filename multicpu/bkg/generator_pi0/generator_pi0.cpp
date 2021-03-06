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

void generator_pi0() {
	TStopwatch t;

	//Define all the histograms
	TCanvas* c2 = new TCanvas("c2", "trial analysis", 1);
	TH1F* h03b = new TH1F("h03b","isolation p^{cone0.4}_{T}+E^{cone0.4}_{T} for bkg;GeV",50, 0, 50);
	TH1F* h04b = new TH1F("h04b","relative isolation (p^{cone0.4}_{T}+E^{cone0.4}_{T})/p_{T} for bkg",50, 0, 5);
	TH1F* h03s = new TH1F("h03s","isolation p^{cone0.4}_{T}+E^{cone0.4}_{T} for signal;GeV",100, 0, 5);
	TH1F* h04s = new TH1F("h04s","relative isolation (p^{cone0.4}_{T}+E^{cone0.4}_{T})/p_{T} for signal",100, 0, 1);

	//Insert the .root file and get the tree and its Event branch
        TFile* f1 = TFile::Open("/data/cpapage/jpsiInclusive-gen8K.root");
        TTree* t1 = dynamic_cast< TTree* >(f1->Get("Pythia"));
        TClonesArray* particles1 = new TClonesArray("Event");
        t1->SetBranchAddress("Event", &particles1);
        for (Long64_t i = 0; i < t1->GetEntries(); ++i) { // This is like the Event loop 
                t1->GetEntry(i);
                for (Long64_t l = 0; l < particles1->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles1->At(l));
			Bool_t cut1 = part->Eta() < 2.4 && part->Eta() > -2.4;
			if (part->GetPdgCode() == 111 && part->Pt()>10 && cut1) {
				Float_t isolation = 0;
				for (Long64_t k = 0; k < particles1->GetEntries(); ++k) {
					TParticle* candidate = reinterpret_cast< TParticle* >(particles1->At(k));
					// Bool_t cut2 = candidate->Eta() < 2.4 && candidate->Eta() > -2.4;
					if (candidate->GetStatusCode() > 0) {
						Float_t DR = DeltaR(part,candidate); 
						TParticlePDG* candpdg = candidate->GetPDG(0);
						if (DR < 0.2 && k!=l && candpdg && k!=part->GetDaughter(0) && k!=part->GetDaughter(1)) {
							Float_t pt = candidate->Pt();
							Float_t mass = candidate->GetMass();
							Float_t Et = sqrt(pt*pt+mass*mass);
							isolation += Et;
							if (candpdg->Charge()!=0) isolation += pt;
						}
					}
				}
			        h03b->Fill(isolation);
				h04b->Fill(isolation/part->Pt());
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
		while (1) {    //Find the last Higgs
			TParticle* part = reinterpret_cast< TParticle* >(particles2->At(iHiggs));
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = reinterpret_cast< TParticle* >(particles2->At(daughter0id));
				if (daughter0->GetPdgCode() == 25) {
					iHiggs = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = reinterpret_cast< TParticle* >(particles2->At(daughter1id));
				if (daughter1->GetPdgCode() == 25) {
					iHiggs = part->GetDaughter(1);
					daughter1validity = false;
					continue;
				}
			}
			if (daughter0validity && daughter1validity) break;
		}
	
		iGamma = iHiggs;
		while (1) {    //Find the last gamma
			TParticle* part = reinterpret_cast< TParticle* >(particles2->At(iGamma));
		        Int_t daughter0id = part->GetDaughter(0);
		        Int_t daughter1id = part->GetDaughter(1);
			Bool_t daughter0validity = true;
			Bool_t daughter1validity = true;
			if (daughter0id >= 0) {
				TParticle* daughter0 = reinterpret_cast< TParticle* >(particles2->At(daughter0id));
				if (daughter0->GetPdgCode() == 22) {
					iGamma = part->GetDaughter(0);
					daughter0validity = false;
					continue;
				}
			}
			if (daughter1id >= 0) {
				TParticle* daughter1 = reinterpret_cast< TParticle* >(particles2->At(daughter1id));
				if (daughter1->GetPdgCode() == 22) {
					iGamma = part->GetDaughter(1);
					daughter1validity = false;
					continue;
				}
			}
			if (daughter0validity && daughter1validity) break;
		}

		TParticle* part = reinterpret_cast< TParticle* >(particles2->At(iGamma));
		Bool_t cut1 = part->Eta() < 2.4 && part->Eta() > -2.4;
		if (part->Pt()>10 && cut1) {
			Float_t isolation = 0;
			for (Long64_t k = 0; k < particles2->GetEntries(); ++k) {
				TParticle* candidate = reinterpret_cast< TParticle* >(particles2->At(k));
				if (candidate->GetStatusCode() > 0) {
					Float_t DeltaPhi;
					Float_t DeltaPhip = part->Phi()-candidate->Phi();
					Float_t DeltaPhim = abs(part->Phi()-candidate->Phi()-TMath::TwoPi());
					Float_t DeltaEta = part->Eta()-candidate->Eta();
					if (DeltaPhip > DeltaPhim) DeltaPhi = DeltaPhim;
					else DeltaPhi = DeltaPhip;
					Float_t DeltaR = sqrt(DeltaPhi*DeltaPhi+DeltaEta*DeltaEta);
					// Bool_t cut2 = candidate->Eta() < 2.4 && candidate->Eta() > -2.4;
					TParticlePDG* candpdg = candidate->GetPDG(0);
					if (DeltaR < 0.2 && k!=iGamma && candpdg) {
						Float_t pt = candidate->Pt();
						Float_t mass = candidate->GetMass();
						Float_t Et = sqrt(pt*pt+mass*mass);
						isolation += Et;
						if (candpdg->Charge()!=0) isolation += pt;
					}
				}
			}
			h03s->Fill(isolation);
			h04s->Fill(isolation/part->Pt());
		}	
	}

	c2->Divide(2,2);
	c2->cd(1); h03b->Draw();
	c2->cd(2); h04b->Draw();
	c2->cd(3); h03s->Draw();
	c2->cd(4); h04s->Draw();

	// TFile* fnew = new TFile("hisol02.root","RECREATE");
	// h03b->Write();
	// h03s->Write();
	// h04b->Write();
	// h04s->Write();
	// fnew->Close();

	t.Stop();
	t.Print();
}
