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

	THStack* hs11 = new THStack("hs11", "p_{T} of J/#psi");
	THStack* hs21 = new THStack("hs21", "y of J/#psi");
	THStack* hs31 = new THStack("hs31", "#eta of J/#psi");
	THStack* hs41 = new THStack("hs41", "E of J/#psi");
	THStack* hs12 = new THStack("hs12", "p_{T} of #gamma");
	THStack* hs22 = new THStack("hs22", "y of #gamma");
	THStack* hs32 = new THStack("hs32", "#eta of #gamma");
	THStack* hs42 = new THStack("hs42", "E of #gamma");
	THStack* hs13 = new THStack("hs13", "p_{T} of #mu");
	THStack* hs23 = new THStack("hs23", "y of #mu");
	THStack* hs33 = new THStack("hs33", "#eta of #mu");
	THStack* hs43 = new THStack("hs43", "E of #mu");

	hs11->Add(h11s); hs11->Add(h11b);
	hs21->Add(h21s); hs21->Add(h21b);
	hs31->Add(h31s); hs31->Add(h31b);
	hs41->Add(h41s); hs41->Add(h41b);
	hs12->Add(h12s); hs12->Add(h12b);
	hs22->Add(h22s); hs22->Add(h22b);
	hs32->Add(h32s); hs32->Add(h32b);
	hs42->Add(h42s); hs42->Add(h42b);
	hs13->Add(h13s); hs13->Add(h13b);
	hs23->Add(h23s); hs23->Add(h23b);
	hs33->Add(h33s); hs33->Add(h33b);
	hs43->Add(h43s); hs43->Add(h43b);
	
	h11s->SetLineColor(kRed); h11s->SetStats(kFALSE); h11b->SetStats(kFALSE);
	h21s->SetLineColor(kRed); h21s->SetStats(kFALSE); h21b->SetStats(kFALSE);
	h31s->SetLineColor(kRed); h31s->SetStats(kFALSE); h31b->SetStats(kFALSE);
	h41s->SetLineColor(kRed); h41s->SetStats(kFALSE); h41b->SetStats(kFALSE);
	h12s->SetLineColor(kRed); h12s->SetStats(kFALSE); h12b->SetStats(kFALSE);
	h22s->SetLineColor(kRed); h22s->SetStats(kFALSE); h22b->SetStats(kFALSE);
	h32s->SetLineColor(kRed); h32s->SetStats(kFALSE); h32b->SetStats(kFALSE);
	h42s->SetLineColor(kRed); h42s->SetStats(kFALSE); h42b->SetStats(kFALSE);
	h13s->SetLineColor(kRed); h13s->SetStats(kFALSE); h13b->SetStats(kFALSE);
	h23s->SetLineColor(kRed); h23s->SetStats(kFALSE); h23b->SetStats(kFALSE);
	h33s->SetLineColor(kRed); h33s->SetStats(kFALSE); h33b->SetStats(kFALSE);
	h43s->SetLineColor(kRed); h43s->SetStats(kFALSE); h43b->SetStats(kFALSE);

	auto l11 = new TLegend(0.65,0.7,0.9,0.9); l11->AddEntry(h11s,"Signal","l"); l11->AddEntry(h11b,"Background","l"); l11->SetTextSize(0.05);
	auto l21 = new TLegend(0.65,0.7,0.9,0.9); l21->AddEntry(h21s,"Signal","l"); l21->AddEntry(h21b,"Background","l"); l21->SetTextSize(0.05);
	auto l31 = new TLegend(0.65,0.7,0.9,0.9); l31->AddEntry(h31s,"Signal","l"); l31->AddEntry(h31b,"Background","l"); l31->SetTextSize(0.05);
	auto l41 = new TLegend(0.65,0.7,0.9,0.9); l41->AddEntry(h41s,"Signal","l"); l41->AddEntry(h41b,"Background","l"); l41->SetTextSize(0.05);
	auto l12 = new TLegend(0.65,0.7,0.9,0.9); l12->AddEntry(h12s,"Signal","l"); l12->AddEntry(h12b,"Background","l"); l12->SetTextSize(0.05);
	auto l22 = new TLegend(0.65,0.7,0.9,0.9); l22->AddEntry(h22s,"Signal","l"); l22->AddEntry(h22b,"Background","l"); l22->SetTextSize(0.05);
	auto l32 = new TLegend(0.65,0.7,0.9,0.9); l32->AddEntry(h32s,"Signal","l"); l32->AddEntry(h32b,"Background","l"); l32->SetTextSize(0.05);
	auto l42 = new TLegend(0.65,0.7,0.9,0.9); l42->AddEntry(h42s,"Signal","l"); l42->AddEntry(h42b,"Background","l"); l42->SetTextSize(0.05);
	auto l13 = new TLegend(0.65,0.7,0.9,0.9); l13->AddEntry(h13s,"Signal","l"); l13->AddEntry(h13b,"Background","l"); l13->SetTextSize(0.05);
	auto l23 = new TLegend(0.65,0.7,0.9,0.9); l23->AddEntry(h23s,"Signal","l"); l23->AddEntry(h23b,"Background","l"); l23->SetTextSize(0.05);
	auto l33 = new TLegend(0.65,0.7,0.9,0.9); l33->AddEntry(h33s,"Signal","l"); l33->AddEntry(h33b,"Background","l"); l33->SetTextSize(0.05);
	auto l43 = new TLegend(0.65,0.7,0.9,0.9); l43->AddEntry(h43s,"Signal","l"); l43->AddEntry(h43b,"Background","l"); l43->SetTextSize(0.05);
	

	c->Divide(2,2);
	c->cd(1); hs11->Draw("nostack"); l11->Draw(); 
	c->cd(2); hs21->Draw("nostack"); l21->Draw();
	c->cd(3); hs31->Draw("nostack"); l31->Draw();
	c->cd(4); hs41->Draw("nostack"); l41->Draw();
	c->Print("hists_overlay.pdf(");
	c->cd(1); hs12->Draw("nostack"); l12->Draw();
	c->cd(2); hs22->Draw("nostack"); l22->Draw();
	c->cd(3); hs32->Draw("nostack"); l32->Draw();
	c->cd(4); hs42->Draw("nostack"); l42->Draw();
	c->Print("hists_overlay.pdf");
	c->cd(1); hs13->Draw("nostack"); l13->Draw();
	c->cd(2); hs23->Draw("nostack"); l23->Draw();
	c->cd(3); hs33->Draw("nostack"); l33->Draw();
	c->cd(4); hs43->Draw("nostack"); l43->Draw();
	c->Print("hists_overlay.pdf)");
}
