{
	TCanvas* c2 = new TCanvas("c2","c2",1);
	TFile* f = TFile::Open("hmass_2.root");
	TH1D* h01b1 = (TH1D*)(f->Get("h01b1"));
	TH1D* h01b2 = (TH1D*)(f->Get("h01b2"));
	TH1D* h01b3 = (TH1D*)(f->Get("h01b3"));
	TH1D* h01b4 = (TH1D*)(f->Get("h01b4"));
	TH1D* h01s = (TH1D*)(f->Get("h01s"));

	float normb1 = 27865.21;
	float normb2 = 1371.25;
	float normb3 = 205;
	float normb4 = 75.55;
	float norms = 0.029131;

	TH1D h01b1_norm = normb1*(*h01b1);
	TH1D h01b2_norm = normb2*(*h01b2);
	TH1D h01b3_norm = normb3*(*h01b3);
	TH1D h01b4_norm = normb4*(*h01b4);
	TH1D h01s_norm = norms*(*h01s);

	TH1D* h01b_norm = (TH1D*)(h01b1_norm.Clone());
	h01b_norm->Add(&h01b2_norm);
	h01b_norm->Add(&h01b3_norm);
	h01b_norm->Add(&h01b4_norm);

	//Find The S/B and S/sqrt(B) ratios
	Double_t S = 0;
	Double_t B = 0;
	for (int i = 275-10; i <= 275+10; ++i) {
		S += h01s_norm.GetBinContent(i);
		B += h01b_norm->GetBinContent(i);
	}

	cout << "B = " << B << endl;
	cout << "S = " << S << endl;
	cout << "S/B = " << S/B << endl;
	cout << "S/sqrt(B) = " << S/sqrt(B) << endl;
	cout << "S_cL = " << sqrt(2*(S+B)*log(1+S/B)-2*S) << endl;

	//Plot the invariant mass of mumugamma
	h01b_norm->Rebin(10);
	h01s_norm.Rebin(10);
	h01s_norm.SetFillColor(kRed);

	THStack* hmass = new THStack("hmass","m_{J/#psi#gamma};m_{J/#psi#gamma} (GeV/c^{2})");
	hmass->Add(h01b_norm);
	hmass->Add(&h01s_norm);

	auto legend1 = new TLegend(0.75,0.7,0.9,0.9); legend1->AddEntry(h01b_norm,"bkg","l"); legend1->AddEntry(&h01s_norm,"signal","f"); 
	legend1->SetTextSize(0.05);

	hmass->Draw("hist"); legend1->Draw();
        c2->Print("generator_mass_2_100K.eps");
}
