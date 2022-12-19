std::vector<double> xbins;
std::vector<std::tuple<TString, TString, EColor, EMarkerStyle>> cuts;
std::vector<std::tuple<TString, TString, EMarkerStyle>> modules;

TList *hlist, *hlist_cc;

void TransformTodNdPt(TH1* h); 
void setXbins(TString);
void setStyle();
void DrawPurities(TString mcPeriod);

void mc_purities(TString mcPeriod = "LHC16g", Bool_t isHardBins = kFALSE, TString sqrts = "13TeV")
{
  setXbins(sqrts);
  setStyle();

  auto f = TFile::Open(Form("../MergedResults_MC/%s/Merged.root", mcPeriod.Data()));
  auto data  = (THashList*)f->Get("Data");
  auto data2 = (THashList*)f->Get("Data2");
  hlist = new TList();
  const Int_t nBins = xbins.size() - 1;

  for (auto cut : cuts ) {
     TString cutName  = get<0>(cut);
     TString cutTitle = get<1>(cut);
     EColor  cutColor = get<2>(cut);
     EMarkerStyle cutMarkerStyle = get<3>(cut);

     auto hGamma2D = (TH2F*)data2->FindObject(Form("hCaloPhotonPdgvsPt_%s", cutName.Data()));
     auto hGamma1D = (TH1D*)hGamma2D->ProjectionX(Form("hGamma1D_%s", cutName.Data()), 4000 + 22 + 1, 4000 + 22 + 1);
     auto hAll1D = (TH1D*)hGamma2D->ProjectionX(Form("hAll1D_%s", cutName.Data()), 1, hGamma2D->GetNbinsY());

     auto hGamma = (TH1D*)hGamma1D->Rebin(nBins, hGamma1D->GetName(), &xbins.front());
     auto hAll    = (TH1D*)hAll1D ->Rebin(nBins, hAll1D->GetName(), &xbins.front());

     TEfficiency *pEff;

     if (TEfficiency::CheckConsistency(*hGamma, *hAll)) {
       pEff = new TEfficiency(*hGamma, *hAll);
     } else {
       printf("histograms %s and %s are not consistent!", hGamma->GetName(), hAll->GetName()); 
       continue;
     }

     auto hPur = (TH1F*)hGamma->Clone(Form("hpurity_%s", cutName.Data()));
     hPur->Divide(hPur, hAll, 1, 1, "b");
     //hEff->Scale(1./acc);
     hlist->Add(hPur);

     pEff->SetName(Form("pEff_%s", cutName.Data()));
     
     pEff->SetTitle(Form("%s;p_{T}, GeV/c;#gamma purity, P_{#gamma}=N_{#gamma}/N_{all}", cutTitle.Data()));
     pEff->SetLineColor(cutColor);
     pEff->SetMarkerColor(pEff->GetLineColor());
     pEff->SetMarkerStyle(cutMarkerStyle);
     pEff->SetLineWidth(1);
     pEff->SetMarkerSize(1.);
     pEff->SetFillColor(0);
     pEff->SetFillStyle(0);
     hlist->Add(pEff);
     
  }

  DrawPurities(Form("%s", mcPeriod.Data()));
  
  gSystem->Exec(Form("mkdir -p figures/Purities/%s/", mcPeriod.Data()));

  gSystem->Exec(Form("mkdir -p ROOT/Purities/%s/", mcPeriod.Data()));

  auto fout = TFile::Open(Form("ROOT/Purities/%s/Purities.root", mcPeriod.Data()), "recreate");

  hlist   ->Write();
  hlist_cc->Write();

  printf("\nSaving output to file  %s\n", fout->GetName());
  fout->Close();

}

void TransformTodNdPt(TH1* h) {
   for (Int_t i = 0; i < h->GetNbinsX(); i++) {
      h->SetBinContent(i + 1, h->GetBinContent(i+1)/h->GetBinWidth(i+1));
      h->SetBinError(i + 1,   h->GetBinError(i+1)/h->GetBinWidth(i+1));
   }
}

void setXbins(TString sqrts)
{
 if (sqrts.CompareTo("13TeV") == 0) {
   xbins = {0.8, 1, 1.2, 1.4, 1.6, 1.8,
            2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 
            4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 
            10, 11, 12, 13, 15, 20, 25};
  } else if (sqrts.CompareTo("7TeV") == 0) {
    xbins = {0.8, 1, 1.2, 1.4, 1.6, 1.8, 
             2, 2.2, 2.4, 2.6, 2.8, 
             3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 
             8, 8.5, 9, 9.5, 10, 11, 12, 13, 15, 20};
  }
}

void setStyle()
{
    cuts.emplace_back("all", "no cuts", kRed, kFullCircle);
    cuts.emplace_back("cpv", "cpv", kBlue, kFullSquare);
    cuts.emplace_back("disp", "disp", kGreen, kFullTriangleUp);
    cuts.emplace_back("both", "cpv+disp", kMagenta, kFullTriangleDown);
       
    modules.emplace_back("", "", kOpenCircle);
    modules.emplace_back("_mod1", "M1", kOpenCircle);
    modules.emplace_back("_mod2", "M2", kOpenSquare);
    modules.emplace_back("_mod3", "M3", kOpenTriangleUp);
    modules.emplace_back("_mod4", "M4", kOpenTriangleDown);
}

void DrawPurities(TString mcPeriod)
{
  gStyle->SetOptTitle(0);  
  gStyle->SetOptStat(0);  

  hlist_cc = new TList();

  auto cc_pur = new TCanvas("cc_pur", "cc_pur", 1200, 600);
  
  cc_pur->SetFixedAspectRatio();
  //cc_pur->SetLogy();
  cc_pur->SetLogx();
  cc_pur->SetGridy();
  cc_pur->SetGridx();  
  hlist->Add(cc_pur);
  
  for (auto cut : cuts ) {
     cc_pur->cd();
     auto pEff2 = (TEfficiency*)hlist->FindObject(Form("pEff_%s", get<0>(cut).Data()));
     pEff2->Draw(Form("%s", cc_pur->GetListOfPrimitives()->GetSize()==0 ? "" : "same"));
     if (cc_pur->GetListOfPrimitives()->GetSize() == 1) {
        gPad->Update();
	pEff2->GetPaintedGraph()->GetHistogram()->SetAxisRange(0.3, 1.3, "Y");
	pEff2->GetPaintedGraph()->GetHistogram()->SetAxisRange(0.95, 25, "X");
	pEff2->Draw();
     }
  }  

  cc_pur->BuildLegend(0.72, 0.62, 0.9, 0.85); 

  auto t2 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
  t2->SetBorderSize(0);
  t2->SetFillColor(kWhite);
  t2->AddText("Purity of the reconstructed photon sample");
  cc_pur->cd();   
  t2->Draw(); 

  hlist_cc->Add(cc_pur);

  gSystem->Exec(Form("mkdir -p figures/Purities/%s/", mcPeriod.Data()));
  cc_pur->Print(Form("figures/Purities/%s/Purity.pdf",  mcPeriod.Data()));
  cc_pur->Print(Form("figures/Purities/%s/Purity.root", mcPeriod.Data())); 
}
