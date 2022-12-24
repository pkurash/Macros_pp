#include "cuts.h"

TList *hlist, *hlist_cc;
Double_t getAcceptance(TString mcPeriod);
void DrawEffs(TString mcPeriod);

void detection_efficiency(TString mcPeriod = "LHC16g", Bool_t isHardBins = kFALSE, TString sqrts = "13TeV")
{ 

  setXbins(sqrts);
  setStyle();
  
  auto f = TFile::Open(Form("../MergedResults_MC/%s/Merged.root", mcPeriod.Data()));
  auto data  = (THashList*)f->Get("Data");
  auto data2 = (THashList*)f->Get("Data2");
 
  hlist = new TList();
  const Int_t nBins = xbins.size() - 1;

  Double_t acc = getAcceptance(Form("%s", mcPeriod.Data()));
  auto hacc =  new TH1F("hacc", "acceptance", 1, 0, 1);
  hlist->Add(hacc);
  hacc->SetBinContent(1, acc);

  auto hGammaMC_2D = (TH2F*)data2->FindObject("hGammaMC_true");
  //hGammaMC_2D->Scale(acc);

  auto hGammaMC_1D = (TH1D*)hGammaMC_2D->ProjectionX("hGammaMC_1D", 71, 170);
  auto hGammaMC    = (TH1D*)hGammaMC_1D->Rebin(nBins, hGammaMC_2D->GetName(), &xbins.front());
  TransformTodNdPt(hGammaMC);

  for (auto cut : cuts) {
     TString cutName  = get<0>(cut);
     TString cutTitle = get<1>(cut);
     EColor  cutColor = get<2>(cut);
     EMarkerStyle cutMarkerStyle = get<3>(cut);

     auto hGamma2D = (TH2F*)data2->FindObject(Form("hCaloPhotonPdgvsPt_%s", cutName.Data()));
     auto hGamma1D = (TH1D*)hGamma2D->ProjectionX(Form("hGamma1D_%s", cutName.Data()), 4000 + 22 + 1, 4000 + 22 + 1);
     auto hGamma   = (TH1D*) hGamma1D->Rebin(nBins, hGamma1D->GetName(), &xbins.front());
     TransformTodNdPt(hGamma);

     TEfficiency *pEff;

     if (TEfficiency::CheckConsistency(*hGamma, *hGammaMC)) {
       pEff = new TEfficiency(*hGamma, *hGammaMC);
     } else {
       printf("histograms %s and %s are not consistent!", hGamma->GetName(), hGammaMC->GetName()); 
       continue;
     }

     auto hEff = (TH1F*)hGamma->Clone(Form("heff_%s", cutName.Data()));
     hEff->Divide(hEff, hGammaMC, 1, 1, "b");
     //hEff->Scale(1./acc);
     hlist->Add(hEff);

     pEff->SetName(Form("pEff_%s", cutName.Data()));
     
     pEff->SetTitle(Form("%s;p_{T}, GeV/c;#gamma detection efficiency #varepsilon_{#gamma}", cutTitle.Data()));
     pEff->SetLineColor(cutColor);
     pEff->SetMarkerColor(pEff->GetLineColor());
     pEff->SetMarkerStyle(cutMarkerStyle);
     pEff->SetLineWidth(1);
     pEff->SetMarkerSize(1.);
     pEff->SetFillColor(0);
     pEff->SetFillStyle(0);
     hlist->Add(pEff);
  }

  DrawEffs(Form("%s", mcPeriod.Data()));
  
  gSystem->Exec(Form("mkdir -p figures/Efficiencies/%s/", mcPeriod.Data()));

  gSystem->Exec(Form("mkdir -p ROOT/Efficiencies/%s/", mcPeriod.Data()));

  auto fout = TFile::Open(Form("ROOT/Efficiencies/%s/Efficiencies.root", mcPeriod.Data()), "recreate");

  hlist   ->Write();
  hlist_cc->Write();

  printf("\nSaving output to file  %s\n", fout->GetName());
  fout->Close();
}

Double_t getAcceptance(TString mcPeriod)
{

 auto mcFile = TFile::Open(Form("ROOT/Efficiency/Acceptance_%s.root", mcPeriod.Data()));
 auto accHist = (TH1F*)mcFile->Get("histAcc");
 auto fitFun = new TF1("fitFun", "[0]", 0, 1.e3);

 fitFun->SetParameter(0, 0.025);

 accHist->Fit(fitFun);

 return fitFun->GetParameter(0);

}


void DrawEffs(TString mcPeriod)
{
  gStyle->SetOptTitle(0);  
  gStyle->SetOptStat(0);  

  hlist_cc = new TList();

  auto cc_eff = new TCanvas("cc_eff", "cc_eff", 1200, 600);
  
  cc_eff->SetFixedAspectRatio();
  //cc_eff->SetLogy();
  cc_eff->SetLogx();
  cc_eff->SetGridy();
  cc_eff->SetGridx();  
  hlist->Add(cc_eff);
  
  for (auto cut : cuts ) {
     cc_eff->cd();
     auto pEff2 = (TEfficiency*)hlist->FindObject(Form("pEff_%s", get<0>(cut).Data()));
     pEff2->Draw(Form("%s", cc_eff->GetListOfPrimitives()->GetSize()==0 ? "" : "same"));
     if (cc_eff->GetListOfPrimitives()->GetSize() == 1) {
        gPad->Update();
	pEff2->GetPaintedGraph()->GetHistogram()->SetAxisRange(0.01, 0.043, "Y");
	pEff2->GetPaintedGraph()->GetHistogram()->SetAxisRange(0.95, 25, "X");
	pEff2->Draw();
     }
  }  

  cc_eff->BuildLegend(0.72, 0.62, 0.9, 0.85); 

  auto t2 = new TPaveText(0.2, 0.92, 0.8, 1.0, "brNDC"); // left-up   
  t2->SetBorderSize(0);
  t2->SetFillColor(kWhite);
  t2->AddText("PHOS #gamma efficiency, N^{#gamma}_{PHOS}/(N^{#gamma}_{gen}(|y| < 0.5, 0 < #phi < 2#pi))");
  cc_eff->cd();   
  t2->Draw(); 

  hlist_cc->Add(cc_eff);

  gSystem->Exec(Form("mkdir -p figures/Efficiencies/%s/", mcPeriod.Data()));
  cc_eff->Print(Form("figures/Efficiencies/%s/Efficiency.pdf",  mcPeriod.Data()));
  cc_eff->Print(Form("figures/Efficiencies/%s/Efficiency.root", mcPeriod.Data())); 
}
