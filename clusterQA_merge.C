void CopyHistBinValues(TH1 *h1, TH1 *h2, Int_t iBB);
void clusterQA_merge(Int_t year = 2016)
{

 
 std::vector<TString> periods;

 if (year == 2016)
   periods = {"LHC16g", "LHC16h", "LHC16i", "LHC16j", "LHC16k", 
              "LHC16l", "LHC16o", "LHC16p"};
 else if (year == 2017)
   periods = {"LHC17c", "LHC17e", "LHC17f", "LHC17h", "LHC17i", 
              "LHC17k", "LHC17l", "LHC17m", "LHC17o", "LHC17r"};
 else if (year == 2018)
   periods = {"LHC18b", "LHC18d", "LHC18e", "LHC18f", "LHC18g", 
              "LHC18h", "LHC18i", "LHC18j", "LHC18k", "LHC18l", 
	      "LHC18m", "LHC18n", "LHC18o", "LHC18p"};
 else 
  {printf("please provide either 2016 2017 or 2018 as an argument\n"); 
   return;}

  TList *outList = new TList();
  const Int_t nPer = periods.size();
  const Int_t nMod = 4;

  TFile *f;
  Int_t na = 0;
  Int_t nl;

  ifstream ff;
 
  for (auto period : periods) {
    ff.open(Form("./datasets/%s_good.txt", period.Data()));
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
  
    hEClusters[iMod]->SetAxisRange(0.1, 0.9, "Y");
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
    f = TFile::Open(Form("ROOT/QA/Cluster_QA_%s_all_good.root", per.Data()));
    histNev= (TH1F*)f->Get("hEvents");
    Int_t nBinOffset = histNev->GetNbinsX();

    printf ("%s& %d & %e & %d-%d \\\\ \n", per.Data(),  histNev->GetNbinsX(),histNev->Integral(), iBin + 1, iBin + nBinOffset);

    CopyHistBinValues(hNEvents, histNev, iBin);
    for (Int_t iMod = 0; iMod < nMod + 1; iMod ++){
       histN = (TH1F*)f->Get(Form("hNClust%d", iMod));
       histE = (TH1F*)f->Get(Form("hEn%d", iMod));
       CopyHistBinValues(hNClusters[iMod], histN, iBin);
       CopyHistBinValues(hEClusters[iMod], histE, iBin);
       hEClusters[iMod]->Sumw2(0);
    }
    f->Close();
    iBin = iBin + nBinOffset;    
  }

  
  printf ("%s & %d & %e & 1-%d \\\\ \n", "Total", hNEvents->GetNbinsX(), hNEvents->Integral(), iBin);
  hNEvents->Sumw2(0);

  TCanvas *cc_clust = new TCanvas("cc_clust", "cc_clust", 700, 600);
  cc_clust->Divide(1, 2);

  cc_clust->cd(1);
//  cc_clust->cd(1)->SetGridx();
//  cc_clust->cd(1)->SetGridy();
  hNEvents->SetAxisRange(0, 1e7, "Y");
  hNEvents->Draw();


  cc_clust->cd(1);
 // cc_clust->cd(1)->SetGridx();
 // cc_clust->cd(1)->SetGridy();
  for (Int_t iMod = 0; iMod < nMod + 1; iMod ++ ) {
     hNClusters[iMod]->Draw(Form("%s", iMod == 0 ? "" : "same"));
  }

  cc_clust->cd(2);
 // cc_clust->cd(2)->SetGridx();
 // cc_clust->cd(2)->SetGridy();
  hEClusters[0]->Draw();

 // hNEvents->Draw();

  auto fout = TFile::Open(Form("ROOT/QA/%s_and_others.root", periods.at(0).Data()), "recreate");

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
