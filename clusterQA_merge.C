void CopyHistBinValues(TH1 *h1, TH1 *h2, Int_t iBB);
void clusterQA_merge()
{
  std::vector<TString> periods = {"LHC16g", "LHC16h", "LHC16i", "LHC16j", "LHC16k", "LHC16l", "LHC16o", "LHC16p"};
  const TString pass = "pass2";

  TList *outList = new TList();
  const Int_t nPer = periods.size();
  const Int_t nMod = 4;

  TFile *f;
  Int_t na = 0;
  Int_t nl;

  ifstream ff;
 
  for (auto period : periods) {
    ff.open(Form("../../PHOSGamma_PP/datasets/%s-%s.txt", period.Data(), pass.Data()));
    nl = 0;
    std::string line;
    while ( std::getline(ff, line) ) {
       nl ++;
    }   
    ff.close();
    na = na + nl;
  }

  const Int_t nAll = na;
  cout << "nAll = " << nAll << endl;

  TH1F *hNEvents = new TH1F("hNEvents", "number of events in run;run index;Number of events", nAll, 0, 1.*nAll);
  outList->Add(hNEvents);

  TH1F *hNClusters[nMod], *hEClusters[nMod];
  for (Int_t iMod = 0; iMod < nMod + 1; iMod ++) {
    hNClusters[iMod] = new TH1F(Form("hNClusters%d", iMod), Form("%s;run index;<N_{cl}>", iMod == 0 ? "all modules" : Form("module %d", iMod + 1)), nAll, 0, 1.*nAll);
    hEClusters[iMod] = new TH1F(Form("hEClusters%d", iMod), Form("%s;run index;<E_{cl}>", iMod == 0 ? "all modules" : Form("module %d", iMod + 1)), nAll, 0, 1.*nAll);
    hNClusters[iMod]->SetAxisRange(0, 0.07, "Y");
    hNClusters[iMod]->SetLineColor(iMod + 1);
    hNClusters[iMod]->SetMarkerColor(iMod + 1);
    hNClusters[iMod]->SetMarkerStyle(kFullCircle);
    hNClusters[iMod]->SetMarkerSize(0.5);
  
    hEClusters[iMod]->SetAxisRange(0, 0.7, "Y");
    hEClusters[iMod]->SetLineColor(iMod + 1);
    hEClusters[iMod]->SetMarkerColor(iMod + 1);
    hEClusters[iMod]->SetMarkerStyle(kFullCircle);
    hEClusters[iMod]->SetMarkerSize(0.5);

    outList->Add(hNClusters[iMod]);
    outList->Add(hEClusters[iMod]);
  }

  Int_t iBin = 0;
  TH1F *histN, *histE, *histNev;
  for (auto per : periods) {
    f = TFile::Open(Form("ROOT/QA/Cluster_QA_%s_all.root", per.Data()));
    cout << Form("ROOT/QA/Cluster_QA_%s_all.root", per.Data()) << endl;
    histNev= (TH1F*)f->Get("hEvents");
    Int_t nBinOffset = histNev->GetNbinsX();
    printf("111!\n");
    CopyHistBinValues(hNEvents, histNev, iBin);
    printf("222!\n");
    for (Int_t iMod = 0; iMod < nMod + 1; iMod ++){
       histN = (TH1F*)f->Get(Form("hNClust%d", iMod));
       histE = (TH1F*)f->Get(Form("hEn%d", iMod));
       CopyHistBinValues(hNClusters[iMod], histN, iBin);
       CopyHistBinValues(hEClusters[iMod], histE, iBin);
    }
    f->Close();
    iBin = iBin + nBinOffset;    
  }

  
  for (Int_t iMod = 0; iMod < nMod + 1; iMod ++) {
  }
  hNEvents->Sumw2(0);

  TCanvas *cc_clust = new TCanvas("cc_clust", "cc_clust", 700, 600);
  cc_clust->Divide(1, 2);

  cc_clust->cd(1);
  cc_clust->cd(1)->SetGridx();
  cc_clust->cd(1)->SetGridy();
  hNEvents->SetAxisRange(0, 1e7, "Y");
  hNEvents->Draw();


  cc_clust->cd(1);
  cc_clust->cd(1)->SetGridx();
  cc_clust->cd(1)->SetGridy();
  for (Int_t iMod = 0; iMod < nMod + 1; iMod ++ ) {
     hNClusters[iMod]->Draw(Form("%s", iMod == 0 ? "" : "same"));
  }

  cc_clust->cd(2);
  cc_clust->cd(2)->SetGridx();
  cc_clust->cd(2)->SetGridy();
  hEClusters[0]->Draw();

 // hNEvents->Draw();

  TFile *fout = TFile::Open(Form("ROOT/QA/%s_and_others.root", periods.at(0).Data()), "recreate");

  cc_clust->Write();
  outList->Write();

  printf("\n Output saved to the file %s \n", fout->GetName());
  fout->Close();
}

void CopyHistBinValues(TH1 *h1, TH1 *h2, Int_t iBB) 
{
   for (Int_t iBin = 0; iBin < h2->GetNbinsX(); iBin ++) {
     h1->SetBinContent(iBin + 1 + iBB, h2->GetBinContent(iBin + 1)); 
     h1->SetBinError(iBin + 1 + iBB, h2->GetBinError(iBin + 1)); 
   }
}
