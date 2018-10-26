#ifndef __CINT__
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParticle.h>
#endif

void generator_check_gamma() {
	//Define all the histograms
	TCanvas* c1 = new TCanvas("c1", "trial analysis", 1);
	TH1F* h011b = new TH1F("h011b","#pi^{0} daughters' PDG codes",1000, 0, 1000);

	//Insert the .root file and get the tree and its Event branch
        TFile* f = TFile::Open("/data/cpapage/jpsiInclusive-gen8K.root");
        TTree* t = dynamic_cast< TTree* >(f->Get("Pythia"));
        TClonesArray* particles = new TClonesArray("Event");
        t->SetBranchAddress("Event", &particles);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) { // This is like the Event loop 
                t->GetEntry(i);
                for (Long64_t l = 0; l < particles->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles->At(l));
			if (part->GetPdgCode() == 111) {
				Float_t id1 = part->GetDaughter(0);
				Float_t id2 = part->GetDaughter(1);
				if (id1 >= 0) {
					TParticle* daughter1 = reinterpret_cast< TParticle* >(particles->At(id1));
					if(daughter1->GetStatusCode()<0) h011b->Fill(abs(daughter1->GetPdgCode()));
				}
				if (id2 >= 0) {
					TParticle* daughter2 = reinterpret_cast< TParticle* >(particles->At(id2));
					if(daughter2->GetStatusCode()<0) h011b->Fill(abs(daughter2->GetPdgCode()));
				}
			}
       		}
	}

	h011b->Draw();
	c1->Print("generator_check_pi0.pdf");
}
