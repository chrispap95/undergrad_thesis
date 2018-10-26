void hists_overlay() {
	TFile* f1 = TFile::Open("signal/signal_hists.root");
	TFile* f2 = TFile::Open("bkg/bkg_hists.root");

	TCanvas* c = new TCanvas("c", "hists_overlay",1);

	TH1F* h11s = (TH1F*)(f1->Get("h11"));
	TH1F* h21s = (TH1F*)(f1->Get("h21"));
	TH1F* h31s = (TH1F*)(f1->Get("h31"));
	TH1F* h41s = (TH1F*)(f1->Get("h41"));
	TH1F* h12s = (TH1F*)(f1->Get("h12"));
	TH1F* h22s = (TH1F*)(f1->Get("h22"));
	TH1F* h32s = (TH1F*)(f1->Get("h32"));
	TH1F* h42s = (TH1F*)(f1->Get("h42"));
	TH1F* h13s = (TH1F*)(f1->Get("h13"));
	TH1F* h23s = (TH1F*)(f1->Get("h23"));
	TH1F* h33s = (TH1F*)(f1->Get("h33"));
	TH1F* h43s = (TH1F*)(f1->Get("h43"));

	TH1F* h11b = (TH1F*)(f2->Get("h11b"));
	TH1F* h21b = (TH1F*)(f2->Get("h21b"));
	TH1F* h31b = (TH1F*)(f2->Get("h31b"));
	TH1F* h41b = (TH1F*)(f2->Get("h41b"));
	TH1F* h12b = (TH1F*)(f2->Get("h12b"));
	TH1F* h22b = (TH1F*)(f2->Get("h22b"));
	TH1F* h32b = (TH1F*)(f2->Get("h32b"));
	TH1F* h42b = (TH1F*)(f2->Get("h42b"));
	TH1F* h13b = (TH1F*)(f2->Get("h13b"));
	TH1F* h23b = (TH1F*)(f2->Get("h23b"));
	TH1F* h33b = (TH1F*)(f2->Get("h33b"));
	TH1F* h43b = (TH1F*)(f2->Get("h43b"));

	h11b->SetBins(100, 0, 200);
	h21s->SetBins(100, -8, 8);
	h22s->SetBins(100, -11, 11);
	h12b->SetBins(100, 0, 200);
	h32s->SetBins(100, -11, 11);
	h42b->SetBins(100, 0, 1000);
	h13b->SetBins(100, 0, 150);
	h43b->SetBins(100, 0, 600);

	TH1F* h11 = (TH1F*)h11s->Clone();
	TH1F* h21 = (TH1F*)h21s->Clone();
	TH1F* h31 = (TH1F*)h31s->Clone();
	TH1F* h41 = (TH1F*)h41s->Clone();
	TH1F* h12 = (TH1F*)h12s->Clone();
	TH1F* h22 = (TH1F*)h22s->Clone();
	TH1F* h32 = (TH1F*)h32s->Clone();
	TH1F* h42 = (TH1F*)h42s->Clone();
	TH1F* h13 = (TH1F*)h13s->Clone();
	TH1F* h23 = (TH1F*)h23s->Clone();
	TH1F* h33 = (TH1F*)h33s->Clone();
	TH1F* h43 = (TH1F*)h43s->Clone();

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

	c->Divide(2,2);
	c->cd(1); h11->Draw(); h11b->Draw("Same"); h11b->SetLineColor(kRed); //h11s->Draw("Same"); h11s->SetLineColor(kBlue);
	c->cd(2); h21->Draw(); h21b->Draw("SAME"); h21b->SetLineColor(kRed); //h21s->Draw("SAME"); h21s->SetLineColor(kBlue);
	c->cd(3); h31->Draw(); h31b->Draw("SAME"); h31b->SetLineColor(kRed); //h31s->Draw("SAME"); h31s->SetLineColor(kBlue);
	c->cd(4); h41->Draw(); h41b->Draw("SAME"); h41b->SetLineColor(kRed); //h41s->Draw("SAME"); h41s->SetLineColor(kBlue);
	c->Print("hists_overlay.pdf(");
	c->cd(1); h12->Draw(); h12b->Draw("SAME"); h12b->SetLineColor(kRed); //h12s->Draw("SAME"); h12s->SetLineColor(kBlue);
	c->cd(2); h22->Draw(); h22b->Draw("SAME"); h22b->SetLineColor(kRed); //h22s->Draw("SAME"); h22s->SetLineColor(kBlue);
	c->cd(3); h32->Draw(); h32b->Draw("SAME"); h32b->SetLineColor(kRed); //h32s->Draw("SAME"); h32s->SetLineColor(kBlue);
	c->cd(4); h42->Draw(); h42b->Draw("SAME"); h42b->SetLineColor(kRed); //h42s->Draw("SAME"); h42s->SetLineColor(kBlue);
	c->Print("hists_overlay.pdf");
	c->cd(1); h13->Draw(); h13b->Draw("SAME"); h13b->SetLineColor(kRed); //h13s->Draw("SAME"); h13s->SetLineColor(kBlue);
	c->cd(2); h23->Draw(); h23b->Draw("SAME"); h23b->SetLineColor(kRed); //h23s->Draw("SAME"); h23s->SetLineColor(kBlue);
	c->cd(3); h33->Draw(); h33b->Draw("SAME"); h33b->SetLineColor(kRed); //h33s->Draw("SAME"); h33s->SetLineColor(kBlue);
	c->cd(4); h43->Draw(); h43b->Draw("Same"); h43b->SetLineColor(kRed); //h43s->Draw("Same"); h43s->SetLineColor(kBlue);
	c->Print("hists_overlay.pdf)");
}
