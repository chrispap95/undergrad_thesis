void sim_hists_overlay() {
	TFile* f1 = TFile::Open("signal/pythia8gensim/signal_sim_hists.root");
	TFile* f2 = TFile::Open("bkg/pythia8gensim/bkg_sim_hists.root");
	TCanvas* c = new TCanvas("c", "hists_overlay",1);

	TH1F* h12s = (TH1F*)(f1->Get("h12"));
	TH1F* h22s = (TH1F*)(f1->Get("h22"));
	TH1F* h32s = (TH1F*)(f1->Get("h32"));
	TH1F* h13s = (TH1F*)(f1->Get("h13"));
	TH1F* h23s = (TH1F*)(f1->Get("h23"));
	TH1F* h33s = (TH1F*)(f1->Get("h33"));
	TH1F* h43s = (TH1F*)(f1->Get("h43"));

	TH1F* h12b = (TH1F*)(f2->Get("h12b"));
	TH1F* h22b = (TH1F*)(f2->Get("h22b"));
	TH1F* h32b = (TH1F*)(f2->Get("h32b"));
	TH1F* h13b = (TH1F*)(f2->Get("h13b"));
	TH1F* h23b = (TH1F*)(f2->Get("h23b"));
	TH1F* h33b = (TH1F*)(f2->Get("h33b"));
	TH1F* h43b = (TH1F*)(f2->Get("h43b"));

	// TH1F* h12b = 4.89*(*h12bo);
	// TH1F* h22b = 4.89*(*h22bo);
	// TH1F* h32b = 4.89*(*h32bo);
	// TH1F* h13b = 4.89*(*h13bo);
	// TH1F* h23b = 8.4*(*h23bo);
	// TH1F* h33b = 8.4*(*h33bo);
	// TH1F* h43b = 8.4*(*h43bo);
        
	THStack* hs12 = new THStack("hs12", "p_{T} of #mu");
	THStack* hs22 = new THStack("hs22", "y of #mu");
	THStack* hs32 = new THStack("hs32", "#eta of #mu");
       	THStack* hs13 = new THStack("hs13", "p_{T} of #gamma");
	THStack* hs23 = new THStack("hs23", "y of #gamma");
	THStack* hs33 = new THStack("hs33", "#eta of #gamma");
	THStack* hs43 = new THStack("hs43", "E of #gamma");

	hs12->Add(h12s); hs12->Add(h12b);
	hs22->Add(h22s); hs22->Add(h22b);
	hs32->Add(h32s); hs32->Add(h32b);
	hs13->Add(h13s); hs13->Add(h13b);
	hs23->Add(h23s); hs23->Add(h23b);
	hs33->Add(h33s); hs33->Add(h33b);
	hs43->Add(h43s); hs43->Add(h43b);
	
	
	// h11b->SetBins(100, 0, 200);
	// h21s->SetBins(100, -8, 8);
	// h22s->SetBins(100, -11, 11);
	// h12b->SetBins(100, 0, 200);
	// h32s->SetBins(100, -11, 11);
	// h42b->SetBins(100, 0, 1000);
	// h13b->SetBins(100, 0, 150);
	// h43b->SetBins(100, 0, 600);

	// TH1F* h11 = (TH1F*)h11s->Clone();
	// TH1F* h21 = (TH1F*)h21s->Clone();
	// TH1F* h31 = (TH1F*)h31s->Clone();
	// TH1F* h41 = (TH1F*)h41s->Clone();
	// TH1F* h12 = (TH1F*)h12s->Clone();
	// TH1F* h22 = (TH1F*)h22s->Clone();
	// TH1F* h32 = (TH1F*)h32s->Clone();
	// TH1F* h42 = (TH1F*)h42s->Clone();
	// TH1F* h13 = (TH1F*)h13s->Clone();
	// TH1F* h23 = (TH1F*)h23s->Clone();
	// TH1F* h33 = (TH1F*)h33s->Clone();
	// TH1F* h43 = (TH1F*)h43s->Clone();

	// h11->Add(h11b);
	// h21->Add(h21b);
	// h31->Add(h31b);
	// h41->Add(h41b);
	// h12->Add(h12b);
	// h22->Add(h22b);
	// h32->Add(h32b);
	// h42->Add(h42b);
	// h13->Add(h13b);
	// h23->Add(h23b);
	// h33->Add(h33b);
	// h43->Add(h43b);

	h12s->SetLineColor(kRed); h12s->SetStats(kFALSE); h12b->SetStats(kFALSE);
	h22s->SetLineColor(kRed); h22s->SetStats(kFALSE); h22b->SetStats(kFALSE);
	h32s->SetLineColor(kRed); h32s->SetStats(kFALSE); h32b->SetStats(kFALSE);
	h13s->SetLineColor(kRed); h13s->SetStats(kFALSE); h13b->SetStats(kFALSE);
	h23s->SetLineColor(kRed); h23s->SetStats(kFALSE); h23b->SetStats(kFALSE);
	h33s->SetLineColor(kRed); h33s->SetStats(kFALSE); h33b->SetStats(kFALSE);
	h43s->SetLineColor(kRed); h43s->SetStats(kFALSE); h43b->SetStats(kFALSE);

	auto l12 = new TLegend(0.65,0.7,0.9,0.9); l12->AddEntry(h12s,"Signal","l"); l12->AddEntry(h12b,"Background","l"); l12->SetTextSize(0.05);
	auto l22 = new TLegend(0.65,0.7,0.9,0.9); l22->AddEntry(h22s,"Signal","l"); l22->AddEntry(h22b,"Background","l"); l22->SetTextSize(0.05);
	auto l32 = new TLegend(0.65,0.7,0.9,0.9); l32->AddEntry(h32s,"Signal","l"); l32->AddEntry(h32b,"Background","l"); l32->SetTextSize(0.05);
	auto l13 = new TLegend(0.65,0.7,0.9,0.9); l13->AddEntry(h13s,"Signal","l"); l13->AddEntry(h13b,"Background","l"); l13->SetTextSize(0.05);
	auto l23 = new TLegend(0.65,0.7,0.9,0.9); l23->AddEntry(h23s,"Signal","l"); l23->AddEntry(h23b,"Background","l"); l23->SetTextSize(0.05);
	auto l33 = new TLegend(0.65,0.7,0.9,0.9); l33->AddEntry(h33s,"Signal","l"); l33->AddEntry(h33b,"Background","l"); l33->SetTextSize(0.05);
	auto l43 = new TLegend(0.65,0.7,0.9,0.9); l43->AddEntry(h43s,"Signal","l"); l43->AddEntry(h43b,"Background","l"); l43->SetTextSize(0.05);
	

	c->Divide(2,2);
	c->cd(1); hs12->Draw("nostack"); l12->Draw();
	c->cd(2); hs22->Draw("nostack"); l22->Draw();
	c->cd(3); hs32->Draw("nostack"); l32->Draw();
	c->Print("sim_hists_overlay.pdf(");
	c->cd(1); hs13->Draw("nostack"); l13->Draw();
	c->cd(2); hs23->Draw("nostack"); l23->Draw();
	c->cd(3); hs33->Draw("nostack"); l33->Draw();
	c->cd(4); hs43->Draw("nostack"); l43->Draw();
	c->Print("sim_hists_overlay.pdf)");
}
