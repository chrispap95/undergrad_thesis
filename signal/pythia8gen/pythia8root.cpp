#include "Pythia.h" 
#include "TPythia8.h" 
#include "TParticle.h" 
#include "TClonesArray.h" 
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TError.h"

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
	}

	pythia->PrintStatistics();		 // print out pythia statistics

	// Save ROOT output to file
	t->Print();
	t->Write();
	f->Close();	
	delete f;
	
	return 0;
} // End main program with error-free return.
