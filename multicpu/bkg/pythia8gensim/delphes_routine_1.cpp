void delphes_routine_1() {
  TFile* f1 = TFile::Open("/data/cpapage/jpsiInclusive-sim8x.root");
  TTree* t1 = (TTree*)(f1->Get("Delphes"));

  TCanvas* c1 = new TCanvas("c1","Delphes Results",1);
  TH1F* h12b = new TH1F("h12b", "p_{T} of #mu; p_{T} (GeV/c)", 100, 0, 100);
  TH1F* h22b = new TH1F("h22b", "#phi of #mu; #phi", 100, -4, 4);
  TH1F* h32b = new TH1F("h32b", "#eta of #mu; #eta", 100, -3, 3);
  TH1F* h13b = new TH1F("h13b", "p_{T} of #gamma; p_{T} (GeV/c)", 100, 0, 150);
  TH1F* h23b = new TH1F("h23b", "#phi of #gamma; #phi", 100, -4, 4);
  TH1F* h33b = new TH1F("h33b", "#eta of #gamma; #eta", 100, -3, 3);
  TH1F* h43b = new TH1F("h43b", "E of #gamma; E (GeV)", 100, 0, 300);

  TText* text1 = new TText(40,100,"nEvents = 8000");
  TText* text2 = new TText(40,120,"Delphes Simulation");
  TLatex text3;

  c1->Divide(2,2);
  c1->cd(1); t1->Draw("Muon.PT>>h12b"); text1->Draw(); text2->Draw(); text3.DrawLatex(40,80,"pp#rightarrow J/#psi#rightarrow#mu^{+}#mu^{-}");
  c1->cd(2); t1->Draw("Muon.Phi>>h22b"); 
  c1->cd(3); t1->Draw("Muon.Eta>>h32b");
  c1->cd(4); 
  c1->Print("sim_jpsi_inclusive.pdf(");
  c1->cd(1); t1->Draw("Photon.PT>>h13b");
  c1->cd(2); t1->Draw("Photon.Phi>>h23b");
  c1->cd(3); t1->Draw("Photon.Eta>>h33b");
  c1->cd(4); t1->Draw("Photon.E>>h43b");
  c1->Print("sim_jpsi_inclusive.pdf)");
  TFile* f2 = new TFile("bkg_sim_hists.root", "RECREATE");
	h12b->Write();
	h22b->Write();
	h32b->Write();
	h13b->Write();
	h23b->Write();
	h33b->Write();
	h43b->Write();
	f2->Close();
}
