void selectGoodRuns(TString period) 
{
  ifstream allff, badff;
  ofstream goodff;
  
  Int_t nn;

  std::vector<Int_t> allNumbers, badNumbers, goodNumbers;

  badff.open(Form("./datasets/%s_bad.txt", period.Data()));
  while (badff >> nn) {
    badNumbers.emplace_back(nn);
  }
  badff.close();

  allff.open(Form("./datasets/%s.txt", period.Data()));
  while (allff >> nn) {
    Bool_t isGood = kTRUE;
    for (auto nb : badNumbers) {
      if (nb - nn == 0) {
        isGood = kFALSE;
      }
    }
    if (!isGood) {
      printf("%d: BAD number!\n", nn);
      continue;
    } else {
      printf("%d: GOOD number!\n", nn);
      goodNumbers.emplace_back(nn);
    } 
  }
  allff.close();

  goodff.open(Form("./datasets/%s_good.txt", period.Data()));
  for (auto ng : goodNumbers) {
    goodff << ng << endl;
  }
  goodff.close();
}
