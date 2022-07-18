
struct PInfo
{
  int FileNum;
  int pdg;
  int pId;
  int evId;
  float rapid;
  int motherId;
  int fDaughterId;
  int lDaughterId;
  int ITSTrackInd = -1;
  int TPCTrackInd = -1;
  int ITSTPCTrackInd = -1;
  int ITSTPCABTrackInd = -1;
  unsigned char isFakeITS = 0u;
  unsigned char isFakeTPC = 0u;
  unsigned char isFakeITSTPC = 0u;
  unsigned char isFakeITSTPCAB = 0u;
  unsigned short clusters = 0u;
  unsigned short clustersITS = 0u;
  std::vector<double> Vcoor;
  float E;
  float M;
  float P;
  std::vector<double> Vmom;
  float T;
  float pt;
};

void pion_detect(std::string fileName_mc = "tf1/sgn_1_Kine.root", std::string fileName_its = "tf1/o2trac_its.root", std::string fileName_tpc = "tf1/tpctracks.root", std::string fileName_itstpc = "tf1/o2match_itstpc.root")
{
  using namespace o2::its;
  using namespace o2::tpc;
  using namespace o2::dataformats;
  using namespace o2::itsmft;

  /////
  o2::base::GeometryManager::loadGeometry("tf1/");
  auto gman = o2::its::GeometryTGeo::Instance();
  TFile::Open("tf1/o2clus_its.root");
  TTree *clusTree = (TTree *)gFile->Get("o2sim");
  std::vector<CompClusterExt> *clusArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterComp", &clusArr);

  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);
  /////LOADING TREES///////////////////////////
  TChain TPCTree("tpcrec");
  TChain ITSTPCTree("matchTPCITS");
  TChain ITSTree("o2sim");
  TChain mcTree("o2sim");

  // auto fileName_mc = "tf1/sgn_1_Kine.root";
  // auto fileName_its = "tf1/o2trac_its.root";
  // auto fileName_tpc = "tf1/tpctracks.root";
  // auto fileName_itstpc = "tf1/o2match_itstpc.root";
  mcTree.AddFile(fileName_mc.data());
  ITSTree.AddFile(fileName_its.data());
  TPCTree.AddFile(fileName_tpc.data());
  ITSTPCTree.AddFile(fileName_itstpc.data());
  std::cout << fileName_mc << std::endl;
  std::cout << fileName_its << std::endl;
  std::cout << fileName_tpc << std::endl;
  std::cout << fileName_itstpc << std::endl;

  std::vector<o2::MCTrack> *mcArr = nullptr;
  std::vector<TrackITS> *itsArr = nullptr;
  std::vector<o2::MCCompLabel> *itsLabArr = nullptr;
  std::vector<TrackTPC> *tpcArr = nullptr;
  std::vector<o2::MCCompLabel> *tpcLabArr = nullptr;
  // std::vector<o2::TrackTPCITS> *itstpcabArr = nullptr;
  std::vector<o2::MCCompLabel> *itstpcabLabArr = nullptr;
  std::vector<TrackTPCITS> *itstpcArr = nullptr;
  std::vector<TrackTPCITS> *itstpcabArr = nullptr;
  std::vector<o2::MCCompLabel> *itstpcLabArr = nullptr;
  std::vector<o2::itsmft::TrkClusRef> *clusRefArr = nullptr;
  std::vector<int> *clusIdArr = nullptr;

  ITSTPCTree.SetBranchAddress("TPCITSABRefs", &clusRefArr);
  ITSTPCTree.SetBranchAddress("TPCITSABCLID", &clusIdArr);
  ITSTree.SetBranchAddress("ITSTrack", &itsArr);
  ITSTree.SetBranchAddress("ITSTrackMCTruth", &itsLabArr);
  TPCTree.SetBranchAddress("TPCTracks", &tpcArr);
  TPCTree.SetBranchAddress("TPCTracksMCTruth", &tpcLabArr);
  ITSTPCTree.SetBranchAddress("TPCITS", &itstpcArr);
  // ITSTPCTree.SetBranchAddress("TPCITSABRefs", &itstpcabArr);
  ITSTPCTree.SetBranchAddress("MatchMCTruth", &itstpcLabArr);
  // ITSTPCTree.SetBranchAddress("TPCITS", &itstpcabArr);
  ITSTPCTree.SetBranchAddress("MatchABMCTruth", &itstpcabLabArr);
  mcTree.SetBranchAddress("MCTrack", &mcArr);
  //////////////////////////////////////////////////

  ////////////Printing Initial Info/////////////////
  int nev = mcTree.GetEntries();
  ITSTPCTree.GetEntry(0);
  TPCTree.GetEntry(0);
  ITSTree.GetEntry(0);
  clusTree->GetEntry(0);

  std::cout << "Number of Events:" << nev;
  std::cout << '\n';
  std::cout << "Number of ITS tracks:" << itsArr->size();
  std::cout << '\n';
  std::cout << "Number of TPC tracks:" << tpcArr->size();
  std::cout << '\n';
  std::cout << "Number of ITSTPC tracks:" << itstpcArr->size() - itstpcabLabArr->size();
  std::cout << '\n';
  std::cout << "Number of ITSTPCAB tracks:" << itstpcabLabArr->size();
  std::cout << '\n';
  ////////////////////////////////////////////////////////////////
  TH1F h1("Check1", "Sigma radius (Track Total); cm; Number", 100, 0, 100);
  // TH1F h1("Check2", "Sigma momentum (Track Pt)", 60, 0, 15);
  // TH1F h2("Check3", "Sigma Eta (Track)", 60, -1.1, 1.1);
  // TH1F h3("Check4", "Sigma Resolution mom. (Track)", 120, -15, 15);
  // TH1F h4("Check5", "Sigma momentum (MC)", 60, 0, 15);
  // TH1F h5("Check6", "Sigma momentum (MC)", 60, 0, 5);
  //////////////Defining variables///////////////////////////////
  std::vector<std::vector<PInfo>> info(nev); //[Number of Files, number of events in that file, PInfo struct]

  //////////////////////////////////////////////////////////////

  ///////////////////Filling Pinfo structure///////////////////
  int counter = 0;
  for (int n = 0; n < 1000; ++n)
  {
    mcTree.GetEntry(n);
    info[n].resize(mcArr->size());
    for (int m = 0; m < mcArr->size(); ++m)
    {
      auto part = mcArr->at(m);
      info[n][m].T = part.T();
      info[n][m].motherId = part.getMotherTrackId();
      info[n][m].pId = m;
      info[n][m].rapid = part.GetRapidity();
      info[n][m].evId = n;
      info[n][m].P = part.GetP();
      info[n][m].pdg = part.GetPdgCode();
      info[n][m].Vmom = {part.Px(), part.Py(), part.Pz()};
      info[n][m].Vcoor = {part.Vx(), part.Vy(), part.Vz()};
      info[n][m].M = part.GetMass();
      info[n][m].E = part.GetEnergy();
      info[n][m].pt = part.GetPt();
      info[n][m].fDaughterId = part.getFirstDaughterTrackId();
      info[n][m].lDaughterId = part.getLastDaughterTrackId();

      if (std::abs(part.GetPdgCode()) == 3112)
        counter++;
    }
  }
  std::cout << " Total Sigma:"
            << " " << counter << std::endl;

  for (int n = 0; n < itsLabArr->size(); ++n) // saves relevant tracks
  {
    auto lab = itsLabArr->at(n);
    if (!lab.isSet())
      continue;

    int trackID, evID, srcID;
    bool fake;
    lab.get(trackID, evID, srcID, fake);

    if (evID < 0 || evID > nev)
      continue;

    if (trackID < 0 || trackID >= (int)info[evID].size())
      continue;

    info[evID][trackID].isFakeITS += fake;
    info[evID][trackID].ITSTrackInd = n;

    // if (std::abs(info[evID][trackID].pdg) == 211 && fake == false){

    // }
    // if (std::abs(info[evID][trackID].pdg) == 3112 && fake == false){

    //   h.Fill(itsArr->at(n).getP());
    //   h1.Fill(itsArr->at(n).getPt());
    //   h2.Fill(itsArr->at(n).getEta());
    //   h3.Fill(info[evID][trackID].P - itsArr->at(n).getP());
    //   h4.Fill(info[evID][trackID].P);

    //   h5.Fill(sqrt(info[evID][trackID].E*info[evID][trackID].E-itsArr->at(n).getP()*itsArr->at(n).getP()));
    // }
  }
  std::cout << "ITS tracks done." << endl;
  for (int n = 0; n < tpcLabArr->size(); ++n) // saves relevant tracks
  {
    auto lab = tpcLabArr->at(n);
    if (!lab.isSet())
      continue;

    int trackID, evID, srcID;
    bool fake;
    lab.get(trackID, evID, srcID, fake);

    if (evID < 0 || evID > nev)
      continue;
    if (trackID < 0 || trackID >= (int)info[evID].size())
      continue;
    info[evID][trackID].isFakeTPC += fake;
    info[evID][trackID].TPCTrackInd = n;
  }
  std::cout << "TPC tracks done." << endl;
  int temo = 0;
  int offset = itstpcLabArr->size() - itstpcabLabArr->size();
  for (int n = 0; n < itstpcArr->size(); ++n) // saves relevant tracks
  {
    auto track = itstpcArr->at(n);
    auto lab = itstpcLabArr->at(n);

    if (!lab.isSet())
      continue;
    int trackID, evID, srcID;
    bool fake;
    lab.get(trackID, evID, srcID, fake);
    if (evID < 0 || evID > nev)
      continue;
    if (trackID < 0 || trackID >= (int)info[evID].size())
      continue;

    auto bad = false;
    if (track.getRefITS().getSource() == 24)
    {
      auto lab0 = itstpcabLabArr->at(n - offset);
      if (!lab0.isSet())
        continue;
      int trackID0, evID0, srcID0;
      bool fake0;
      lab0.get(trackID0, evID0, srcID0, fake0);

      auto clusRef = clusRefArr->at(n - offset);
      // int it = clusRef.getFirstEntry(), itLim = it + clusRef.getEntries();
      // auto bad = false;
      // for (; it < itLim; it++)
      // {
      //   auto tvid = clusIdArr->at(it);
      //   auto labCl = (clusLabArr->getLabels(tvid))[0];
      //   if (!labCl.isValid() || labCl.getSourceID() != 0 || !labCl.isCorrect())
      //     continue;

      //   int trackIDCl, evIDCl, srcIDCl;
      //   bool fakeCl;
      //   labCl.get(trackIDCl, evIDCl, srcIDCl, fakeCl);

      // if (evID0!=evID || trackID != trackID0)//trackID0 != trackID || trackID != trackIDCl evIDCl!=evID || trackID != trackIDCl
      // {
      //   std::cout << "Bad: " <<it<< std::endl;
      //   temo++;
      //   bad = true;
      //   break;
      // }
      //}
      // if (evID0!=evID || trackID != trackID0)//trackID0 != trackID || trackID != trackIDCl evIDCl!=evID || trackID != trackIDCl
      //   {
      //     bad = true;
      //     std::cout << "Bad: " <<n<< std::endl;
      //     break;

      //   }
      if (bad == false)
      {
        info[evID][trackID].isFakeITSTPCAB += fake;
        info[evID][trackID].ITSTPCABTrackInd = n;
      }
    }
    else if (track.getRefITS().getSource() == 0)
    {

      info[evID][trackID].isFakeITSTPC += fake;
      info[evID][trackID].ITSTPCTrackInd = n;
    }
  }
  std::cout << "Bad: " << temo << std::endl;
  std::cout << "ITSTPC tracks done." << endl;

  for (int frame = 0; frame < clusTree->GetEntriesFast(); frame++)
  { // Cluster frames
    if (!clusTree->GetEvent(frame))
      continue;

    for (unsigned int iClus{0}; iClus < clusArr->size(); ++iClus)
    {
      auto lab = (clusLabArr->getLabels(iClus))[0];
      if (!lab.isValid() || lab.getSourceID() != 0 || !lab.isCorrect())
        continue;

      int trackID, evID, srcID;
      bool fake;
      lab.get(trackID, evID, srcID, fake);
      if (evID < 0 || evID >= (int)info.size())
      {
        std::cout << "Cluster MC label eventID out of range" << std::endl;
        continue;
      }
      if (trackID < 0 || trackID >= (int)info[evID].size())
      {
        std::cout << "Cluster MC label trackID out of range" << std::endl;
        continue;
      }

      const CompClusterExt &c = (*clusArr)[iClus];
      auto layer = gman->getLayer(c.getSensorID());
      // if (evID == n && trackID == k)
      //{
      info[evID][trackID].clustersITS |= 1 << layer;
      //}
    }
  }


  TFile file("/home/justas_t/SigmaP/pion_detect_007.root", "recreate");
  TTree tree("PInfo", "PInfo");
  PInfo pInfo;
  tree.Branch("particle", &pInfo);
  for (auto &event : info)
  {
    for (auto &part : event)
    {

      pInfo = part;
      if (std::abs(pInfo.pdg) != 211 && std::abs(pInfo.pdg) != 3222 &&std::abs(pInfo.pdg) != 3112) //&& std::abs(pInfo.pdg) != 2112) //save only pion and sigma
        continue;
      tree.Fill();
    }
  }
  tree.Write();
  // h.Write();
  // h1.Write();
  //  h2.Write();
  //  h3.Write();
  //  h4.Write();
  //  h5.Write();
  file.Close();
  int tempcounter123 = 0;
  int tempcounter123fake = 0;
  // int tempcounter123 = 0;
  for (int n = 0; n < 1000; n++)
  {
    for (int m = 0; m < info[n].size(); m++)
    {

      if (info[n][m].ITSTPCTrackInd >= 0 && info[n][m].ITSTPCABTrackInd >= 0)
      {
        tempcounter123++;
        // std::cout << "Particle: " << info[n][m].pdg << std::endl;

        auto clusRef = clusRefArr->at(info[n][m].ITSTPCABTrackInd - offset);
        int it = clusRef.getFirstEntry(), itLim = it + clusRef.getEntries();

        for (; it < itLim; it++)
        {
          auto tvid = clusIdArr->at(it);
          auto labCl = (clusLabArr->getLabels(tvid))[0];
          if (!labCl.isValid() || labCl.getSourceID() != 0 || !labCl.isCorrect())
            continue;

          int trackID, evID, srcID;
          bool fake;
          labCl.get(trackID, evID, srcID, fake);

          if (evID < 0 || evID >= (int)info.size() || n != evID)
          {
            std::cout << "Cluster MC label eventID out of range" << std::endl;
            continue;
          }
          if (trackID < 0 || trackID >= (int)info[evID].size() || m != trackID)
          {
            std::cout << "Cluster MC label trackID out of range" << std::endl;
            continue;
          }
          const CompClusterExt &c = (*clusArr)[tvid];
          auto layer = gman->getLayer(c.getSensorID());
          // if (evID == n && trackID == k)
          //{
          info[n][m].clusters |= 1 << layer;
        }
        int nCl{0};
        for (unsigned int bit{0}; bit < sizeof(info[n][m].clusters) * 8; ++bit)
        {
          nCl += bool(info[n][m].clusters & (1 << bit));
        }
        std::cout << "Cluster hits " << nCl << std::endl;
        if (info[n][m].isFakeITSTPC > 0 || info[n][m].isFakeITSTPCAB > 0)
        {
          tempcounter123fake++;
        }
      }
    }
  }
  std::cout << "Overlap " << tempcounter123 << std::endl;
  std::cout << "Overlap Fake: " << tempcounter123fake << std::endl;
  // std::cout<<"Particele Coor: "<< info[0][15].Vcoor[0]<<" "<< info[0][15].Vcoor[1]<<" "<< info[0][15].Vcoor[2]<<endl;
}
