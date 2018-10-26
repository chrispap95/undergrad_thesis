{
	TFile* f = TFile::Open("dphi.root");
	TH1D* hb1 = (TH1D*)(f->Get("h01b1"));
	TH1D* hb2 = (TH1D*)(f->Get("h01b2"));
	TH1D* hb3 = (TH1D*)(f->Get("h01b3"));
	TH1D* hb4 = (TH1D*)(f->Get("h01b4"));
	TH1D* hs = (TH1D*)(f->Get("h01s"));

	TH1D hb1_norm = 27865.21*(*hb1);
	TH1D hb2_norm = 1371*(*hb2);
	TH1D hb3_norm = 205*(*hb3);
	TH1D hb4_norm = 75.55*(*hb4);
	TH1D hs_norm = 0.02913*(*hs);

	TH1D* hb_norm = (TH1D*)(hb1_norm.Clone());
	hb_norm->Add(&hb2_norm);
	hb_norm->Add(&hb3_norm);
	hb_norm->Add(&hb4_norm);

	TLegend* legend = new TLegend(0.5,0.6,0.7,0.8); 
	legend->AddEntry(hb_norm,"bkg","l");
	legend->AddEntry(&hs_norm,"signal","l");
	legend->SetTextSize(0.05);

	hs_norm.SetTitle("#Delta#phi(#mu#mu,#gamma); #Delta#phi");
        hb_norm->SetTitle("#Delta#phi(#mu#mu,#gamma); #Delta#phi");

        gStyle->SetOptStat(0000);
        gStyle->SetPalette(kBird);
        TCanvas* c = new TCanvas("c","c",1);
        hs_norm.Draw("hist plc");
        hs_norm.GetYaxis()->SetTitle("Signal");
        c->Update();
        Float_t rightmax = 1.1*hb_norm->GetMaximum();
        Float_t scale = gPad->GetUymax()/rightmax;
        hb_norm->Scale(scale);
        hb_norm->Draw("hist same plc");
        TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
        axis->SetTextColor(kYellow+2);
        axis->SetLineColor(kYellow);
        axis->SetLabelColor(kYellow+1);
        axis->SetLabelSize(0.033);
        axis->SetTextSize(0.033);
        axis->SetTitle("Background");
        axis->Draw();
	TLine line;
	line.SetLineColor(kRed);
	line.SetLineStyle(2);
	line.DrawLine(0.5,0,0.5,1.16);
	line.DrawLine(-0.5,0,-0.5,1.16);
	TBox* box = new TBox(-0.5,0,0.5,1.16);
	box->SetFillColorAlpha(kRed-7,0);
	box->SetFillStyle(3751);
	box->Draw();
	TLatex text;
	text.SetTextColor(kRed);
	text.SetTextSize(0.035);
	text.DrawLatex(-0.85,-0.05,"-0.5");
	text.DrawLatex(0.25,-0.05,"0.5");
	c->Print("dphi.eps");
}
