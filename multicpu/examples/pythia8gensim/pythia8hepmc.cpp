#include "Pythia.h" 
#include "Pythia8Plugins/HepMC2.h"
#include "Event.h" 
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
	Pythia pythia;	// Declare Pythia object
	pythia.readFile(argv[1]);	// Read Pythia configuration from data card
	int nEvent = pythia.mode("Main:numberOfEvents");	// set number of events to generate
	pythia.init(); // Initialize Pythia
	
	// Book Pythia histograms
	Hist pT("top transverse momentum", 100, 0., 200.);
	Hist eta("top pseudorapidity", 100, -5., 5.);

	// Setup HepMC output	
	HepMC::GenEvent* hepmcevt;
	HepMC::Pythia8ToHepMC ToHepMC; // Interface for conversion from Pythia8::Event to HepMC event.
	HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out); // Specify file where HepMC events will be stored.
	
	// Generate event(s). (aka "Event loop")
	int prog, oldprog = 0;	// progress indicator variables
	for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
		// print out progress 
		prog = (int)(iEvent * 100. / nEvent);
		if (prog > oldprog) {
			oldprog = prog;
			cout << prog << "\% complete..." << endl;
		}
		
		
		if (! pythia.next()) {
			cout << "skipping event due to error" << endl;
			continue; // Generate an(other) event. Fill event record.
		}
		int iTop = 0;
		for (int i = 0; i < pythia.event.size(); ++i) {	// Particle loop
//			cout << "i = " << i << ", id = " << pythia.event[i].id() << endl;
			if (pythia.event[i].id() == 6) {
				iTop = i;
			}
		}
		if (iTop != 0) {			
//			cout << "top  pT: " << pythia.event[iTop].pT() << endl;
//			cout << "top eta: " << pythia.event[iTop].eta() << endl;
			// Fill histograms
			pT.fill(pythia.event[iTop].pT());
			eta.fill(pythia.event[iTop].eta());			
		}
		
		// Save HepMC event to file
		// Construct new empty HepMC event and fill it.
		// Units will be as chosen for HepMC build, but can be changed
		// by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
		hepmcevt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(pythia, hepmcevt);		

		ascii_io << hepmcevt;	// Write the HepMC event to file. Done with it.
		delete hepmcevt;
	}	
	cout << pT << eta;	 // print out pythia histograms
	pythia.stat();		 // print out pythia statistics
	
	return 0;
} // End main program with error-free return.
