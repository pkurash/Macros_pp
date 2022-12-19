void draw_acceptances()
{

  std::vector<TString> mc2016 = {"LHC16g", "LHC16h", "LHC16i", "LHC16j", "LHC16k", "LHC16l", "LHC16o", "LHC16p"};
  std::vector<TString> mc2017 = {"LHC17c", "LHC17e", "LHC17f", "LHC17h", "LHC17i", "LHC17k",  "LHC17l",  
                                 "LHC17m", "LHC17o", "LHC17r"};
  std::vector<TString> mc2018 = {"LHC18b", "LHC18d", "LHC18e", "LHC18f", "LHC18g", "LHC18h", "LHC18i", "LHC18j", 
                                 "LHC18k", "LHC18l", "LHC18m", "LHC18n", "LHC18o", "LHC18p"};

  const Int_t nPeriods = mc2016.size() +  mc2017.size() +  mc2018.size();

  auto *hAccPeriods = new TH1D("hAccPeriods", "PHOS #gamma acceptance, E_{#gamma} > 0.3 GeV",
                                  nPeriods, 0, 1.0*nPeriods);

  hAccPeriods->GetYaxis()->SetRangeUser(0, 0.032);

  Int_t iPer = 0;

  for (auto mc : mc2016) {
    auto f = TFile::Open(Form("ROOT/Efficiencies/%s/Efficiencies.root", mc.Data()));
    auto hacc = (TH1D*)f->Get("hacc");
    Double_t acc = hacc->GetBinContent(1);
    hAccPeriods->SetBinContent(iPer + 1, acc);
    hAccPeriods->GetXaxis()->SetBinLabel(iPer + 1, Form("%s", mc.Data()));
    iPer++;
    delete f;
  }
  for (auto mc : mc2017) {
    auto f = TFile::Open(Form("ROOT/Efficiencies/%s/Efficiencies.root", mc.Data()));
    auto hacc = (TH1D*)f->Get("hacc");
    Double_t acc = hacc->GetBinContent(1);
    hAccPeriods->SetBinContent(iPer + 1, acc);
    hAccPeriods->GetXaxis()->SetBinLabel(iPer + 1, Form("%s", mc.Data()));
    iPer++;
    delete f;
  }
  for (auto mc : mc2018) {
    auto f = TFile::Open(Form("ROOT/Efficiencies/%s/Efficiencies.root", mc.Data()));
    auto hacc = (TH1D*)f->Get("hacc");
    Double_t acc = hacc->GetBinContent(1);
    hAccPeriods->SetBinContent(iPer + 1, acc);
    hAccPeriods->GetXaxis()->SetBinLabel(iPer + 1, Form("%s", mc.Data()));
    iPer++;
    delete f;
  }

  auto cc = new TCanvas("cc", "cc", 1200, 400);
  cc->cd();
  cc->SetGridx();
  cc->SetGridy();
  hAccPeriods->Draw();
  
  auto fout = TFile::Open("ROOT/AccPeriods.root", "recreate");
  
  hAccPeriods->Write();
  cc->Write();
  
  printf("\nsaving output to gile %s\n", fout->GetName());

  fout->Close();

}

