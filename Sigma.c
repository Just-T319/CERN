#include <cmath>
using namespace o2::its;
using namespace std;

double exp(double *t, double *par)
{
  double arg = 0;
  arg = t[0] / par[1];
  double fitval = par[0] * TMath::Exp(-1.0 * arg);
  return fitval;
}

void Sigma()
{

  TChain trackTree("o2sim");
  TChain mcTree("o2sim");

  for (int l{1}; l <= 29; ++l)
  {
    auto fileName_mc = Form("tf%i/sgn_%i_Kine.root", l, l);
    auto fileName = Form("tf%i/o2trac_its.root", l);
    mcTree.AddFile(fileName_mc);
    trackTree.AddFile(fileName);
    std::cout << fileName_mc << std::endl;
    std::cout << fileName << std::endl;
  }

  std::vector<TrackITS> *trackArr = nullptr;
  std::vector<o2::MCCompLabel> *trackLabArr = nullptr;
  std::vector<o2::MCTrack> *mcArr = nullptr;

  trackTree.SetBranchAddress("ITSTrack", &trackArr);
  trackTree.SetBranchAddress("ITSTrackMCTruth", &trackLabArr);
  mcTree.SetBranchStatus("*", 0);
  mcTree.SetBranchStatus("MCTrack*", 1);
  mcTree.SetBranchAddress("MCTrack", &mcArr);

  int nev = mcTree.GetEntries();
  int ntrack = trackTree.GetEntries();

  std::cout << "Number of events:" << nev;
  std::cout << '\n';
  std::cout << "Number of tracks sets:" << ntrack;
  std::cout << '\n';

  std::vector<double> ct;
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> R;
  std::vector<double> ctt;
  std::vector<double> rapid;
  std::vector<double> PT;

  int goodTr = 0;
  int nsigmam = 0;
  int nsigmamanti = 0;
  int fileNum = 0;
  int fake_tr = 0;
  int fake_sigma = 0;
  int FileNumPrev = -1;
  int trace = 0;
  int tempcounter = 0;

  TH1F h1("CT", "ct spectrum; cm", 180, 0, 90);
  TH1F h2("Y", "Rapidity; y", 250, -2, 2);
  TH1F h3("Pt", "Transverse P; GeV/c", 20, 0, 5);
  TH1F h5("CTT", "ct spectrum (time dil); cm", 180, 0, 90);
  TH1F h7("R", "Radius distribution (xy plane); cm", 250, 0, 250);
  TH1F h8("Hits", "Hits on Layers; Layers; Number of hits", 8, 0, 8);
  TH2F h6("Secondary V", "Vertex distribution; cm ; cm", 60, -40, 40, 60, -40, 40);

  for (int n = 0; n < nev; ++n) // over events
  {
    fileNum = n / 1000; // get file number for tracks
    std::vector<double> trackInfo;
    std::vector<double> trackIndex;

    mcTree.GetEntry(n);
    trackTree.GetEntry(fileNum);

    for (int iTrack = 0; iTrack < trackLabArr->size(); ++iTrack) // saves relevant tracks
    {
      auto lab = trackLabArr->at(iTrack);
      if (!lab.isSet())
        continue;

      int trackID, evID, srcID;
      bool fake;
      lab.get(trackID, evID, srcID, fake);

      if (evID < 0 || (evID + fileNum * 1000) != n) // +fileNum*1000 to account for different files
        continue;
      if (trackID < 0 || trackID >= mcArr->size())
        continue;
      if (fake)
      {
        fake_tr++;
        // continue;
      }
      trackInfo.push_back(trackID);
      trackIndex.push_back(iTrack);
    }

    for (int m = 0; m < mcArr->size(); m++) // over particles
    {
      float sigmaP = 0., sigmaM = 0., sigmaE = 0., sigmaT = 0.;
      int neutron = 0, pion = 0;
      auto track = mcArr->at(m);

      //if (std::abs(track.GetRapidity()) > 1) continue; // reject not in central rapidity
      if (std::abs(track.GetPdgCode()) == 3112) nsigmamanti++;
      if (std::abs(track.GetPdgCode()) == 3222) // check for sigma plus
      {
        float dX = 0., dY = 0., dZ = 0., L = 0., daughterT = 0.;
        int counter = 0;

        counter = track.getLastDaughterTrackId() - track.getFirstDaughterTrackId() + 1;
        auto p = std::find(std::begin(trackInfo), std::end(trackInfo), m);
        int index = std::distance(trackInfo.begin(), p);
        nsigmam++;

        if (track.leftTrace(0))
        {
          trace++;
        }

        sigmaP = track.GetP();
        sigmaT = track.T();
        sigmaM = track.GetMass();
        sigmaE = track.GetEnergy();

        h2.Fill(track.GetRapidity());
        h3.Fill(track.GetPt());

        if (p != std::end(trackInfo))
        {
          // std::cout << "Hits reconstructed Sigma, index:  " << trackInfo[index] << endl;

          int trackID, evID, srcID;
          bool fake;
          trackLabArr->at(trackIndex[index]).get(trackID, evID, srcID, fake);

          int nhits = 0;
          for (int q = 0; q < 7; ++q)
          {
            if (trackArr->at(trackIndex[index]).hasHitOnLayer(q))
            {
              h8.Fill(q + 1);
              nhits++;
              // std::cout << "Layer" << q + 1 << ": true" << endl;
            }
          }
          if (fake)
          {
            //std::cout << "Particle Id: " << mcArr->at(trackID).GetPdgCode() << endl;
            fake_sigma++;
          }
          if (nhits > 3)
          {
            goodTr++;
          }
        }

        for (int i = track.getFirstDaughterTrackId(); i <= track.getLastDaughterTrackId(); ++i)
        {

          if (i < 0)
            break;
          if (mcArr->at(i).GetPdgCode() == 11)
            counter--;

          // std::cout << "n: "<< n << " Particle Id: " << mcArr->at(i).GetPdgCode() << " Entry Sigma: " << m << " Entry Daughter: " << i << endl;
          if (std::abs(mcArr->at(i).GetPdgCode()) == 211)
          {
            dX = mcArr->at(i).Vx() - track.Vx();
            dY = mcArr->at(i).Vy() - track.Vy();
            dZ = mcArr->at(i).Vz() - track.Vz();
            L = sqrt(dX * dX + dY * dY + dZ * dZ);
            daughterT = mcArr->at(i).T();
            pion++;
            // std::cout<< ((daughterT-sigmaT)*sigmaM)/sigmaE;
            // std::cout << '\n';
          }
          else if (std::abs(mcArr->at(i).GetPdgCode()) == 2112)
          {
            neutron++;
          }
        }
        if (neutron == pion && neutron == 1 && counter == 2)
        {
          h1.Fill((L * sigmaM) / sigmaP);
          auto vx = dX + track.Vx();
          auto vy = dY + track.Vy();
          h6.Fill(vx, vy);
          h5.Fill(((daughterT - sigmaT) * sigmaM * 100) / sigmaE);
          h7.Fill(sqrt(vx * vx + vy * vy));
          tempcounter++;
        }
      }
    }
  }

  std::cout << "Total number of tracks: " << trackLabArr->size() << endl;
  std::cout << "Total number of trace by sigma: " << trace << endl;
  std::cout << "Total Sigma minus tracks:" << goodTr;
  std::cout << '\n';
  std::cout << "Sigma fake tracks:" << fake_sigma;
  std::cout << '\n';
  std::cout << "Total fake tracks:" << fake_tr;
  std::cout << '\n';
  std::cout << "Total Sigma plus Baryons:" << nsigmam;
  std::cout << '\n';
  std::cout << "Total Sigma minus anti Baryons:" << nsigmamanti;
  std::cout << '\n';
  std::cout << "Total Sigma resulting in corect decay:" << tempcounter;
  std::cout << '\n';

  TF1 *func = new TF1("fit", exp, 0, 90, 2);
  func->SetParameters(20, h1.GetMean());
  func->SetParNames("Pre_Factor", "Mean_value");
  h1.Fit("fit");
  h1.Chisquare(func);

  TFile outputFile("/home/justas_t/SigmaP/sigmaplus_anlz_007.root", "recreate");

  h1.Write();
  h2.Write();
  h3.Write();
  h5.Write();
  h6.Write();
  h7.Write();
  h8.Write();

  outputFile.Close();

  /*

        h2.Fit("expo");
        h2.DrawCopy();
  */
}