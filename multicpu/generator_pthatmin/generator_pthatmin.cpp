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

void generator_pthatmin() {
	TStopwatch t;

	//Define all the histograms
	TCanvas* c2 = new TCanvas("c2", "trial analysis", 1);
	// TH1I* h20 = new TH1I("h20","p_{T} of J/#psi for pTHatMin = 20;pT (GeV/c)",100, 0, 140);
	// TH1I* h60 = new TH1I("h60","p_{T} of J/#psi for pTHatMin = 60;pT (GeV/c)",100, 0, 140);
	// TH1I* h80 = new TH1I("h80","p_{T} of J/#psi for pTHatMin = 80;pT (GeV/c)",100, 0, 140);
	// THStack* hpt = new THStack("hpt","p_{T} of J/#psi");
	// TH1I* h21 = new TH1I("h21","p_{T} of #gamma for pTHatMin = 20;pT (GeV/c)",100, 5, 35);
	// TH1I* h61 = new TH1I("h61","p_{T} of #gamma for pTHatMin = 60;pT (GeV/c)",100, 5, 35);
	// TH1I* h81 = new TH1I("h81","p_{T} of #gamma for pTHatMin = 80;pT (GeV/c)",100, 5, 35);
	// THStack* hptph = new THStack("hptph","p_{T} of #gamma");
	TH1I* h02 = new TH1I("h02","#Delta#phi_{J/#psi#gamma} for signal;#Delta#phi",100, -10, 10);
	TH1I* h22 = new TH1I("h22","#Delta#phi_{J/#psi#gamma} for pTHatMin = 20;#Delta#phi",100, -10, 10);
	TH1I* h52 = new TH1I("h52","#Delta#phi_{J/#psi#gamma} for pTHatMin = 50;#Delta#phi",100, -10, 10);
	TH1I* h82 = new TH1I("h82","#Delta#phi_{J/#psi#gamma} for pTHatMin = 80;#Delta#phi",100, -10, 10);
	THStack* hdphi = new THStack("hdphi","#Delta#phi_{J/#psi#gamma};#Delta#phi");
	TH1I* h03 = new TH1I("h03","#Delta R_{J/#psi#gamma} for signal;#Delta R",100, 0, 10);
	TH1I* h23 = new TH1I("h23","#Delta R_{J/#psi#gamma} for pTHatMin = 20;#Delta R",100, 0, 10);
	TH1I* h53 = new TH1I("h53","#Delta R_{J/#psi#gamma} for pTHatMin = 50;#Delta R",100, 0, 10);
	TH1I* h83 = new TH1I("h83","#Delta R_{J/#psi#gamma} for pTHatMin = 80;#Delta R",100, 0, 10);
	THStack* hdr = new THStack("hdr","#Delta R_{J/#psi#gamma};#Delta R");
	
	//pTHatMin = 20
	//Insert the .root file and get the tree and its Event branch
        TFile* f1 = TFile::Open("/data/cpapage/jpsiInclusive-gen8K.root");
        TTree* t1 = dynamic_cast< TTree* >(f1->Get("Pythia"));
        TClonesArray* particles1 = new TClonesArray("Event");
        t1->SetBranchAddress("Event", &particles1);
        for (Long64_t i = 0; i < t1->GetEntries(); ++i) { // This is the Event loop 
		t1->GetEntry(i);
		Int_t iJpsi = 0;
                for (Long64_t l = 0; l < particles1->GetEntries(); ++l) { //This is a particle loop
			TParticle* part = reinterpret_cast< TParticle* >(particles1->At(l));
			if (part->GetPdgCode() == 443 && iJpsi == 0) {
				iJpsi = l;
				break;
			}
		}
		iJpsi = FindLastParticle(iJpsi, 443, particles1);
		TParticle* jpsi = reinterpret_cast< TParticle* >(particles1->At(iJpsi));
		Bool_t cut1 = jpsi->Eta() < 2.4 && jpsi->Eta() > -2.4;
                for (Long64_t l = 0; l < particles1->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles1->At(l));
			Bool_t cut2 = part->Eta() < 2.4 && part->Eta() > -2.4;
			Float_t px = part->Px() + jpsi->Px();
			Float_t py = part->Py() + jpsi->Py();
			Float_t pz = part->Pz() + jpsi->Pz();
			Float_t energy = part->Energy() + jpsi->Energy();
			Float_t mass = sqrt(energy*energy-px*px-py*py-pz*pz);
			if(part->GetPdgCode() == 111 && cut1 && cut2 && mass > 90) {
				h22->Fill(jpsi->Phi()-part->Phi());
				h23->Fill(DeltaR(jpsi, part));
			}
		}
	}

	//pTHatMin = 50
	//Insert the .root file and get the tree and its Event branch
        TFile* f3 = TFile::Open("/data/cpapage/jpsiInclusive-gen-pt50.root");
        TTree* t3 = dynamic_cast< TTree* >(f3->Get("Pythia"));
        TClonesArray* particles3 = new TClonesArray("Event");
        t3->SetBranchAddress("Event", &particles3);
        for (Long64_t i = 0; i < t3->GetEntries(); ++i) { // This is the Event loop 
                t3->GetEntry(i);
		Int_t iJpsi = 0;
                for (Long64_t l = 0; l < particles3->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles3->At(l));
			if (part->GetPdgCode() == 443 && iJpsi == 0) {
				iJpsi = l;
				break;
			}
		}
		iJpsi = FindLastParticle(iJpsi, 443, particles3);
		TParticle* jpsi = reinterpret_cast< TParticle* >(particles3->At(iJpsi));
		Bool_t cut1 = jpsi->Eta() < 2.4 && jpsi->Eta() > -2.4;
                for (Long64_t l = 0; l < particles3->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles3->At(l));
			Bool_t cut2 = part->Eta() < 2.4 && part->Eta() > -2.4;
			Float_t px = part->Px() + jpsi->Px();
			Float_t py = part->Py() + jpsi->Py();
			Float_t pz = part->Pz() + jpsi->Pz();
			Float_t energy = part->Energy() + jpsi->Energy();
			Float_t mass = sqrt(energy*energy-px*px-py*py-pz*pz);
			if(part->GetPdgCode() == 111 && cut1 && cut2 && mass > 90) {
				h52->Fill(jpsi->Phi()-part->Phi());
				h53->Fill(DeltaR(jpsi, part));
			}
		} 
	}

	//pTHatMin = 80
	//Insert the .root file and get the tree and its Event branch
        TFile* f4 = TFile::Open("/data/cpapage/jpsiInclusive-gen-pt80.root");
        TTree* t4 = dynamic_cast< TTree* >(f4->Get("Pythia"));
        TClonesArray* particles4 = new TClonesArray("Event");
        t4->SetBranchAddress("Event", &particles4);
        for (Long64_t i = 0; i < t4->GetEntries(); ++i) { // This is the Event loop 
                t4->GetEntry(i);
		Int_t iJpsi = 0;
                for (Long64_t l = 0; l < particles4->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles4->At(l));
			if (part->GetPdgCode() == 443 && iJpsi == 0) {
				iJpsi = l;
				break;
			}
		}
		iJpsi = FindLastParticle(iJpsi, 443, particles4);
		TParticle* jpsi = reinterpret_cast< TParticle* >(particles4->At(iJpsi));
		Bool_t cut1 = jpsi->Eta() < 2.4 && jpsi->Eta() > -2.4;
                for (Long64_t l = 0; l < particles4->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles4->At(l));
			Bool_t cut2 = part->Eta() < 2.4 && part->Eta() > -2.4;
			Float_t px = part->Px() + jpsi->Px();
			Float_t py = part->Py() + jpsi->Py();
			Float_t pz = part->Pz() + jpsi->Pz();
			Float_t energy = part->Energy() + jpsi->Energy();
			Float_t mass = sqrt(energy*energy-px*px-py*py-pz*pz);
			if(part->GetPdgCode() == 111 && cut1 && cut2 && mass > 90) {
				h82->Fill(jpsi->Phi()-part->Phi());
				h83->Fill(DeltaR(jpsi, part));
			}
		}
	}

	//Signal
	//Insert the .root file and get the tree and its Event branch
        TFile* f5 = TFile::Open("/data/cpapage/higgsTojpsiGamma-gen.root");
        TTree* t5 = dynamic_cast< TTree* >(f5->Get("Pythia"));
        TClonesArray* particles5 = new TClonesArray("Event");
        t5->SetBranchAddress("Event", &particles5);
        for (Long64_t i = 0; i < t5->GetEntries(); ++i) { // This is the Event loop 
                t5->GetEntry(i);
		Int_t iHiggs = 0;
                for (Long64_t l = 0; l < particles5->GetEntries(); ++l) { //This is a particle loop
                        TParticle* part = reinterpret_cast< TParticle* >(particles5->At(l));
			if (part->GetPdgCode() == 25 && iHiggs == 0) {
				iHiggs = l;
				break;
			}
		}
		iHiggs = FindLastParticle(iHiggs, 25, particles5);
		Int_t iJpsi = FindLastParticle(iHiggs, 443, particles5);
		Int_t iGamma = FindLastParticle(iHiggs, 22, particles5);
		TParticle* jpsi = reinterpret_cast< TParticle* >(particles5->At(iJpsi));
		TParticle* gamma = reinterpret_cast< TParticle* >(particles5->At(iGamma));
		Bool_t cut1 = jpsi->Eta() < 2.4 && jpsi->Eta() > -2.4;
		Bool_t cut2 = gamma->Eta() < 2.4 && gamma->Eta() > -2.4;
		if(cut1 && cut2) {
			h02->Fill(jpsi->Phi()-gamma->Phi());
			h03->Fill(DeltaR(jpsi, gamma));
		}
	}

	// h22->Scale(3000/h22->GetEntries());
	// h52->Scale(3000/h52->GetEntries());
	// h82->Scale(3000/h82->GetEntries());
	// h23->Scale(3000/h23->GetEntries());
	// h53->Scale(3000/h53->GetEntries());
	// h83->Scale(3000/h83->GetEntries());
	h02->SetLineColor(kOrange); //h02->SetFillColor(kMagenta);
	h52->SetLineColor(kGreen); 
	h82->SetLineColor(kViolet);
	h03->SetLineColor(kOrange); //h03->SetFillColor(kMagenta);
	h53->SetLineColor(kGreen);
	h83->SetLineColor(kViolet);
	hdphi->Add(h02);
	hdphi->Add(h22);
	hdphi->Add(h52);
	hdphi->Add(h82);
	hdr->Add(h03);
	hdr->Add(h23);
	hdr->Add(h53);
	hdr->Add(h83);

	auto legend = new TLegend(0.7,0.65,0.9,0.9); legend->SetHeader("pTHatmin", "C");
	legend->AddEntry(h22, "20", "l"); legend->AddEntry(h52, "50", "l"); legend->AddEntry(h82, "80", "l");
	legend->AddEntry(h02, "signal", "l");
	legend->SetTextSize(0.04); legend->SetTextAlign(22); legend->SetMargin(0.3);

	TLatex text1; text1.SetTextSize(0.03);

	hdphi->Draw("nostack hist"); legend->Draw(); 
	// text1.DrawLatex(-9,140,"all 3 bkg hists are scaled to 3000 entries");
	// text1.DrawLatex(-9,130,"signal has 1000 entries");
	c2->Print("generator_pthatmin_angular_check_n.pdf(");
	hdr->Draw("nostack hist"); legend->Draw();
	c2->Print("generator_pthatmin_angular_check_n.pdf)");

	t.Stop();
	t.Print();
}
