#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "TGeoGlobalMagField.h"
#define HomogeneousField
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"

#include "Framework/ConfigParamRegistry.h"

#include "Field/MagneticField.h"
#include "DataFormatsParameters/GRPObject.h"
#include <TClonesArray.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TPad.h>
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "ITSBase/GeometryTGeo.h"
#include "SimulationDataFormat/TrackReference.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include <string>

using namespace std;

struct PInfo
{
  int FileNum;
  int pdg;
  int pId;
  int evId;
  float rapid;
  float eta;
  float phi;
  int motherId;
  int fDaughterId;
  int lDaughterId;
  int ITSTrackInd = -1;
  int TPCTrackInd = -1;
  int ITSTPCTrackInd = -1;
  int ITSTPCABTrackInd = -1;
  int ITSPVInd = -1;
  int TPCPVInd = -1;
  int ITSTPCPVInd = -1;
  int ITSTPCABPVInd = -1;
  unsigned char isFakeITS = 0u;
  unsigned char isFakeTPC = 0u;
  unsigned char isFakeITSTPC = 0u;
  unsigned char isFakeITSTPCAB = 0u;
  unsigned short clusters = 0u;
  std::vector<double> Vcoor;
  float E;
  float R;
  float M;
  float P;
  std::vector<double> Vmom;
  float T;
  float pt;
};

void DCAFit_Random()
{
  // namespace o2
  // {
  //   namespace vertexing
  //   {
  using namespace o2;
  using namespace vertexing;
  using namespace o2::its;
  using namespace o2::tpc;
  using namespace o2::dataformats;
  using GIndex = o2::dataformats::VtxTrackIndex;
  using namespace o2::itsmft;
  using namespace o2::gpu::gpustd;

  string FITTEROPTION = "DCA"; // "DCA_false" or "KFParticle"

  const auto grp = o2::parameters::GRPObject::loadFrom("tf1/o2sim_grp.root");
  o2::base::GeometryManager::loadGeometry("tf1/");
  auto gman = o2::its::GeometryTGeo::Instance();
  TFile::Open("tf1/o2clus_its.root");
  TTree *clusTree = (TTree *)gFile->Get("o2sim");
  std::vector<CompClusterExt> *clusArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterComp", &clusArr);

  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);
  // // LOG(info) << "Magnetic field loaded from GRP " <<

  TChain TPCTree("tpcrec");
  TChain ITSTPCTree("matchTPCITS");
  TChain ITSTree("o2sim");
  TChain pTree("PInfo");
  TChain mcTree("o2sim");
  TChain pvTree("o2sim");

  auto fileName = "~/pion_detect_007.root";
  auto fileName_pv = "tf1/o2_primary_vertex.root";
  auto fileName_mc = "tf1/sgn_1_Kine.root";
  auto fileName_its = "tf1/o2trac_its.root";
  auto fileName_tpc = "tf1/tpctracks.root";
  auto fileName_itstpc = "tf1/o2match_itstpc.root";
  ITSTree.AddFile(fileName_its);
  TPCTree.AddFile(fileName_tpc);
  ITSTPCTree.AddFile(fileName_itstpc);
  pTree.AddFile(fileName);
  mcTree.AddFile(fileName_mc);
  pvTree.AddFile(fileName_pv);

  std::vector<TrackITS> *itsArr = nullptr;
  std::vector<o2::MCCompLabel> *itsLabArr = nullptr;
  std::vector<TrackTPC> *tpcArr = nullptr;
  std::vector<o2::MCCompLabel> *tpcLabArr = nullptr;
  std::vector<o2::MCCompLabel> *itstpcabLabArr = nullptr;
  std::vector<TrackTPCITS> *itstpcArr = nullptr;
  std::vector<o2::MCCompLabel> *itstpcLabArr = nullptr;
  std::vector<o2::MCTrack> *mcArr = nullptr;
  std::vector<o2::dataformats::VtxTrackIndex> *pvIndexArr = nullptr;
  std::vector<o2::dataformats::VtxTrackRef> *pvRefArr = nullptr;
  std::vector<o2::MCEventLabel> *pvMcArr = nullptr;
  PInfo *pArr = nullptr;

  pTree.SetBranchAddress("particle", &pArr);
  mcTree.SetBranchAddress("MCTrack", &mcArr);
  ITSTree.SetBranchAddress("ITSTrack", &itsArr);
  ITSTree.SetBranchAddress("ITSTrackMCTruth", &itsLabArr);
  TPCTree.SetBranchAddress("TPCTracks", &tpcArr);
  TPCTree.SetBranchAddress("TPCTracksMCTruth", &tpcLabArr);
  ITSTPCTree.SetBranchAddress("TPCITS", &itstpcArr);
  ITSTPCTree.SetBranchAddress("MatchMCTruth", &itstpcLabArr);
  ITSTPCTree.SetBranchAddress("MatchABMCTruth", &itstpcabLabArr);
  pvTree.SetBranchAddress("PVTrackIndices", &pvIndexArr);
  pvTree.SetBranchAddress("PV2TrackRefs", &pvRefArr);
  pvTree.SetBranchAddress("PVMCTruth", &pvMcArr);

  ITSTPCTree.GetEntry(0);
  TPCTree.GetEntry(0);
  ITSTree.GetEntry(0);
  pvTree.GetEntry(0);

  int nev = 1000;
  int nparticles = pTree.GetEntries();
  int nIndex = pvIndexArr->size();
  int nRef = pvRefArr->size();
  int nMc = pvMcArr->size();
  std::cout << "PV Index: " << nIndex << std::endl;
  std::cout << "PV Ref: " << nRef << std::endl;
  std::cout << "PV Mc: " << nMc << std::endl;
  // std::vector<std::vector<int>> motherList(nev);
  int motherList = -1;
  int motherTempID = -1;
  std::vector<std::vector<float>> motherfake(nev);

  // std::vector<std::vector<float>> motherTempID(nev);
  std::vector<std::vector<int>> motherListTrackType(nev);
  // std::vector<PInfo> info(nparticles);
  std::vector<std::vector<PInfo>> info(1000);

  int tempcount = 0;

  TH1F h("DCA", "; #Chi^{2} of s+b associations ; Counts", 65, -15, 50);
  TH1F h1("R", "; |Gen Radius vector- Recon Radius vector|/Gen Radius of s+b; Counts", 40, 0, 0.2); // rnd associations
  TH2F h2("Resol-Chi", "; #Chi^{2} of rnd associations; (Gen Radius - Recon Radius)/gen Radius of s+b associations", 16, 0, 8, 30, -1.1, 1.5);
  TH2F h2R("Resol-R", "; Generated Radius of random; (Gen Radius - Recon Radius)/gen Radius of s+b", 100, 0, 100, 30, -1.1, 1.5);
  TH1F hmass("mass_sigma", ";Invariant #Sigma^{-} M (Gev/#it{c^{2}}) of s+b; Counts per 50 Mev/#it{c^{2}}", 60, 1, 1.3);
  // TH1F h3("Hits_gen", "; Pion ITS hits generated; Count", 8, -0.5, 7.5);
  TH1F h3("Pion R", "; Reconstructed pion radius; Count", 100, 0, 1000);

  TH1F h4("Hits_rec", "; Pion ITS hits reconstructed; Count", 8, -0.5, 7.5);
  TH1F hmassResol("Mass Resol", ";Reconstructed - Generated #Sigma^{-} M (Gev/#it{c^{2}}) of s+b; Counts", 100, -0.6, 0.6);
  TH1F hnResol("Neutron Resol", ";Generated - Reconstructed neutron momentum Gev/c^2; Counts", 40, -2, 2);

  // TH1F hmassResol("mass_sigma_resol", ";Generated - Reconstructed #Sigma^{-} M (Gev/#it{c^{2}}) of rnd associations; Counts per 100 Mev/#it{c^{2}}", 100, 0, 10);
  TH2F h7("Eta-Phi", "; Eta Sigma - Pion; Phi Sigma - Pion", 50, -5, 5, 30, -3, 3);
  TH2F h6("Resol-Gen", "; Resolution; Generated R", 20, -10, 1, 50, 0, 50);
  TH1F h10("CheckMom", ";Track rec - KF rec momentum GeV/c of s+b; Counts", 80, -10, 10);
  TH2F h11("Chi-momentum", "; Chi2; Track rec - KF rec momentum GeV/c of s+b", 55, -5, 50, 40, -2, 2);
  TH2F h12("Mass Resol - chi", ";Chi2; Reconstructed - Generated #Sigma^{-} M (Gev/#it{c^{2}}) of s+b", 55, -5, 50, 40, -1, 5);

  for (int n = 0; n < 1000; ++n)
  {
    mcTree.GetEntry(n);
    info[n].resize(mcArr->size());
    for (unsigned int m = 0; m < mcArr->size(); ++m)
    {
      auto part = mcArr->at(m);
      info[n][m].motherId = part.getMotherTrackId();
      info[n][m].pId = m;
      info[n][m].evId = n;
      info[n][m].phi = part.GetPhi();
      info[n][m].eta = part.GetEta();
      info[n][m].P = part.GetP();
      info[n][m].M = part.GetMass();
      info[n][m].Vcoor = {part.Vx(), part.Vy(), part.Vz()};
      info[n][m].R = sqrt(part.Vx() * part.Vx() + part.Vy() * part.Vy());
      info[n][m].pdg = part.GetPdgCode();
      info[n][m].fDaughterId = part.getFirstDaughterTrackId();
      info[n][m].lDaughterId = part.getLastDaughterTrackId();
    }
  }

  ///////////////////Cluster/////////////////////////
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
      info[evID][trackID].clusters |= 1 << layer;
      //}
    }
  }
  ///////////////////////////////////////////////////

  for (unsigned int n = 0; n < itsLabArr->size(); ++n) // saves relevant tracks
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
    if (std::abs(info[evID][trackID].pdg) == 3112 || std::abs(info[evID][trackID].pdg) == 211)
    {
      info[evID][trackID].isFakeITS += fake;
      info[evID][trackID].ITSTrackInd = n;
    }
  }
  // std::cout << "ITS tracks done." << endl;
  // for (unsigned int n = 0; n < tpcLabArr->size(); ++n) // saves relevant tracks
  // {
  //   auto lab = tpcLabArr->at(n);
  //   if (!lab.isSet())
  //     continue;

  //   int trackID, evID, srcID;
  //   bool fake;
  //   lab.get(trackID, evID, srcID, fake);

  //   if (evID < 0 || evID > nev)
  //     continue;
  //   if (trackID < 0 || trackID >= (int)info[evID].size())
  //     continue;

  //   if (std::abs(info[evID][trackID].pdg) == 211)
  //   {
  //     info[evID][trackID].isFakeTPC += fake;
  //     info[evID][trackID].TPCTrackInd = n;
  //   }
  // }
  // std::cout << "TPC tracks done." << endl;
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
    if (std::abs(info[evID][trackID].pdg) != 211)
      continue;
    if (track.getRefITS().getSource() == 24)
    {
      info[evID][trackID].isFakeITSTPCAB += fake;
      info[evID][trackID].ITSTPCABTrackInd = n;
    }
    else if (track.getRefITS().getSource() == 0)
    {

      info[evID][trackID].isFakeITSTPC += fake;
      info[evID][trackID].ITSTPCTrackInd = n;
    }
  }
  int counterPion = 0;
  int diffID = 0;
  int sameID = 0;
  for (int iv = 0; iv < nRef - 1; iv++) // Ignore the last entry as it is for unassigned tracks
  {
    const auto &vtref = pvRefArr->at(iv);
    int evId = pvMcArr->at(iv).getEventID();
    int it = vtref.getFirstEntry(), itLim = it + vtref.getEntries();
    for (; it < itLim; it++)
    {
      auto tvid = pvIndexArr->at(it);
      int trackID, evID, srcID;
      bool fake;
      if (!tvid.isPVContributor())
        continue;

      if (tvid.getSource() == GIndex::TPC)
      {
        auto lab = tpcLabArr->at(tvid.getIndex());
        lab.get(trackID, evID, srcID, fake);
        // if (evID < 0 || evID > nev) continue;
        // if (evId != evID) std::cout<<"Event Id: "<< evID<<std::endl;
        if (evId != evID)
        {
          diffID++;
          continue;
        }
        sameID++;
        info[evID][trackID].TPCPVInd = it;
        if (std::abs(info[evID][trackID].pdg) == 211)
          counterPion++;

        // continue;
      }
      else if (tvid.getSource() == GIndex::ITS)
      {
        auto lab = itsLabArr->at(tvid.getIndex());
        lab.get(trackID, evID, srcID, fake);
        // if (evID < 0 || evID > nev) continue;
        // if (evId != evID) std::cout<<"Event Id: "<< evID<<std::endl;
        if (evId != evID)
        {
          diffID++;
          continue;
        }
        sameID++;
        info[evID][trackID].ITSPVInd = it;
        if (std::abs(info[evID][trackID].pdg) == 211)
          counterPion++;

        // continue;
      }
      else if (tvid.getSource() == GIndex::ITSTPC)
      {
        auto lab = itstpcLabArr->at(tvid.getIndex());
        lab.get(trackID, evID, srcID, fake);
        // if (evID < 0 || evID > nev) continue;
        // if (evId != evID) std::cout<<"Event Id: "<< evID<<std::endl;
        if (evId != evID)
        {
          diffID++;
          continue;
        }
        sameID++;
        info[evID][trackID].ITSTPCPVInd = it;
        if (std::abs(info[evID][trackID].pdg) == 211)
          counterPion++;

        // continue;
      }
      else if (tvid.getSource() == GIndex::ITSAB)
      {
        std::cout << tvid.getIndex() << std::endl;
        auto lab = itstpcabLabArr->at(tvid.getIndex());
        lab.get(trackID, evID, srcID, fake);
        // if (evID < 0 || evID > nev) continue;
        // if (evId != evID) std::cout<<"Event Id: "<< evID<<std::endl;
        if (evId != evID)
        {
          diffID++;
          continue;
        }
        sameID++;
        info[evID][trackID].ITSTPCABPVInd = it;
        if (std::abs(info[evID][trackID].pdg) == 211)
          counterPion++;

        // continue;
      }
    }
  }

  std::cout << "Pions from PV: " << counterPion << std::endl;
  std::cout << "PV with different ID: " << diffID << std::endl; // tvid.getIndex()
  std::cout << "PV with same ID: " << sameID << std::endl;      // tvid.getIndex()
  std::cout << "Primary Vertex done..." << std::endl;

  std::vector<int> failList;

  for (int n = 0; n < 1000; ++n)
  {
    // std::vector<int> wasList;
    mcTree.GetEntry(n);
    for (unsigned int m = 0; m < mcArr->size(); ++m)
    { //&& tempcount < 100000

      if (std::abs(info[n][m].pdg) == 3112 && ((info[n][m].ITSTPCABTrackInd >= 0 && info[n][m].isFakeITSTPCAB <= 0) || (info[n][m].ITSTPCTrackInd >= 0 && info[n][m].isFakeITSTPC <= 0) || (info[n][m].ITSTrackInd >= 0 && info[n][m].isFakeITS <= 0) || (info[n][m].TPCTrackInd >= 0 && info[n][m].isFakeTPC <= 0)))
      {
        std::vector<int> wasList;
        auto track = mcArr->at(m);
        int neutron = 0;
        int pion = 0;
        int tempID;
        float neutronPFull = 0;
        int daughcount = track.getLastDaughterTrackId() - track.getFirstDaughterTrackId() + 1;
        for (int iDaug = track.getFirstDaughterTrackId(); iDaug <= track.getLastDaughterTrackId(); iDaug++)
        {
          if (iDaug < 0)
            break;
          if (mcArr->at(iDaug).GetPdgCode() == 11)
            daughcount--;

          if (std::abs(mcArr->at(iDaug).GetPdgCode()) == 211)
          {
            pion++;
            tempID = iDaug;
          }
          if (std::abs(mcArr->at(iDaug).GetPdgCode()) == 2112)
            neutron++;
          auto neutronP1 = info[n][iDaug].Vmom;
          neutronPFull = info[n][iDaug].P; // sqrt(neutronP1[0]*neutronP1[0]+neutronP1[1]*neutronP1[1]+neutronP1[2]*neutronP1[2]);
        }
        if (daughcount == 2 && neutron == 1 && pion == 1)
        {
          auto pionR = info[n][tempID].Vcoor; // sigma actually but this is to get where the sigma decays

          for (unsigned int k = 0; k < mcArr->size(); ++k)
          {
            // if (mcArr->at(k).isPrimary())
            // continue;
            auto p = std::find(std::begin(wasList), std::end(wasList), k);
            if (p == std::end(wasList))
            { //&& info[n][k].motherId == m  std::abs(info[n][info[n][k].motherId].pdg) != 3112
              if (std::abs(info[n][k].pdg) == 211 && ((info[n][k].ITSTPCABTrackInd >= 0 && info[n][k].isFakeITSTPCAB <= 0 && info[n][k].ITSTPCPVInd < 0) || (info[n][k].ITSTPCTrackInd >= 0 && info[n][k].isFakeITSTPC <= 0 && info[n][k].ITSTPCPVInd < 0) || (info[n][k].ITSTrackInd >= 0 && info[n][k].isFakeITS <= 0 && info[n][k].ITSPVInd < 0) || (info[n][k].TPCTrackInd >= 0 && info[n][k].isFakeTPC <= 0 && info[n][k].TPCPVInd < 0)))
              {
                // wasList.push_back(k);

                o2::track::TrackParCov SigmaTr;
                o2::track::TrackParCov PionTr;
                std::array<float, 3> sigmaP = {0, 0, 0};
                std::array<float, 3> pionP = {0, 0, 0};
                std::array<float, 3> sigmaV = {0, 0, 0};
                std::array<float, 3> pionV = {0, 0, 0};
                float pionPabs = 0;
                float sigmaPabs = 0;
                float etaS = 0;
                float phiS = 0;
                float etaP = 0;
                float phiP = 0;
                if ((info[n][m].ITSTrackInd >= 0 && info[n][m].isFakeITS <= 0))
                {

                  SigmaTr = itsArr->at(info[n][m].ITSTrackInd).getParamOut(); //.getParamOut();
                  itsArr->at(info[n][m].ITSTrackInd).getPxPyPzGlo(sigmaP);
                  sigmaPabs = itsArr->at(info[n][m].ITSTrackInd).getP();
                  etaS = itsArr->at(info[n][m].ITSTrackInd).getEta();
                  phiS = itsArr->at(info[n][m].ITSTrackInd).getPhi();
                }
                else if ((info[n][m].TPCTrackInd >= 0 && info[n][m].isFakeTPC <= 0))
                {
                  SigmaTr = tpcArr->at(info[n][m].TPCTrackInd).getParamOut();
                  tpcArr->at(info[n][m].TPCTrackInd).getPxPyPzGlo(sigmaP);
                  sigmaPabs = tpcArr->at(info[n][m].TPCTrackInd).getP();
                  etaS = tpcArr->at(info[n][m].TPCTrackInd).getEta();
                  phiS = tpcArr->at(info[n][m].TPCTrackInd).getPhi();
                }
                else if ((info[n][m].ITSTPCTrackInd >= 0 && info[n][m].isFakeITSTPC <= 0))
                {
                  SigmaTr = itstpcArr->at(info[n][m].ITSTPCTrackInd);
                  itstpcArr->at(info[n][m].ITSTPCTrackInd).getPxPyPzGlo(sigmaP);
                  sigmaPabs = itstpcArr->at(info[n][m].ITSTPCTrackInd).getP();
                  etaS = itstpcArr->at(info[n][m].ITSTPCTrackInd).getEta();
                  phiS = itstpcArr->at(info[n][m].ITSTPCTrackInd).getPhi();
                }
                else if ((info[n][m].ITSTPCABTrackInd >= 0 && info[n][m].isFakeITSTPCAB <= 0))
                {
                  SigmaTr = itstpcArr->at(info[n][m].ITSTPCABTrackInd);
                  itstpcArr->at(info[n][m].ITSTPCABTrackInd).getPxPyPzGlo(sigmaP);
                  sigmaPabs = itstpcArr->at(info[n][m].ITSTPCABTrackInd).getP();
                  etaS = itstpcArr->at(info[n][m].ITSTPCABTrackInd).getEta();
                  phiS = itstpcArr->at(info[n][m].ITSTPCABTrackInd).getPhi();
                }

                if ((info[n][k].ITSTrackInd >= 0 && info[n][k].isFakeITS <= 0 && info[n][k].ITSPVInd < 0))
                {
                  PionTr = itsArr->at(info[n][k].ITSTrackInd); //.getParamOut(); //.getParamOut();
                  pionPabs = itsArr->at(info[n][k].ITSTrackInd).getP();
                  itsArr->at(info[n][k].ITSTrackInd).getPxPyPzGlo(pionP);
                  etaP = itsArr->at(info[n][k].ITSTrackInd).getEta();
                  phiP = itsArr->at(info[n][k].ITSTrackInd).getPhi();
                }
                else if ((info[n][k].TPCTrackInd >= 0 && info[n][k].isFakeTPC <= 0 && info[n][k].TPCPVInd < 0))
                {
                  PionTr = tpcArr->at(info[n][k].TPCTrackInd); //.getParamOut();
                  pionPabs = tpcArr->at(info[n][k].TPCTrackInd).getP();
                  tpcArr->at(info[n][k].TPCTrackInd).getPxPyPzGlo(pionP);
                  etaP = tpcArr->at(info[n][k].TPCTrackInd).getEta();
                  phiP = tpcArr->at(info[n][k].TPCTrackInd).getPhi();
                }
                else if ((info[n][k].ITSTPCTrackInd >= 0 && info[n][k].isFakeITSTPC <= 0 && info[n][k].ITSTPCPVInd < 0))
                {
                  PionTr = itstpcArr->at(info[n][k].ITSTPCTrackInd); //.getParamOut();
                  pionPabs = itstpcArr->at(info[n][k].ITSTPCTrackInd).getP();
                  itstpcArr->at(info[n][k].ITSTPCTrackInd).getPxPyPzGlo(pionP);
                  etaP = itstpcArr->at(info[n][k].ITSTPCTrackInd).getEta();
                  phiP = itstpcArr->at(info[n][k].ITSTPCTrackInd).getPhi();
                }
                else if ((info[n][k].ITSTPCABTrackInd >= 0 && info[n][k].isFakeITSTPCAB <= 0 && info[n][k].ITSTPCPVInd < 0))
                {
                  PionTr = itstpcArr->at(info[n][k].ITSTPCABTrackInd); //.getParamOut();
                  pionPabs = itstpcArr->at(info[n][k].ITSTPCABTrackInd).getP();
                  itstpcArr->at(info[n][k].ITSTPCABTrackInd).getPxPyPzGlo(pionP);
                  etaP = itstpcArr->at(info[n][k].ITSTPCABTrackInd).getEta();
                  phiP = itstpcArr->at(info[n][k].ITSTPCABTrackInd).getPhi();
                }

                int nCl{0};
                int nClR{0};
                for (unsigned int bit{0}; bit < sizeof(info[n][k].clusters) * 8; ++bit)
                {
                  nCl += bool(info[n][k].clusters & (1 << bit));
                }
                if ((info[n][k].ITSTPCABTrackInd >= 0 && info[n][k].isFakeITSTPCAB <= 0 && info[n][k].ITSTPCPVInd < 0) || (info[n][k].ITSTPCTrackInd >= 0 && info[n][k].isFakeITSTPC <= 0 && info[n][k].ITSTPCPVInd < 0))
                {
                  for (unsigned int bit{0}; bit < sizeof(info[n][k].clusters) * 8; ++bit)
                  {
                    nClR += bool(info[n][k].clusters & (1 << bit));
                  }
                  h4.Fill(nClR);
                }

                // h3.Fill(nCl);

                // try
                //{ // Magnetic Field 5
                // if (FITTEROPTION == "DCA")
                // {

                //}
                //}

                if (FITTEROPTION == "KFParticle")
                {
                  KFParticle::SetField(grp->getNominalL3Field()); //(int)(grp->getNominalL3Field())
                  KFParticle kfSigma = KFParticle();
                  KFParticle kfPion = KFParticle();

                  // double posmomS[6];

                  // kfNeutron.SetField(grp->getNominalL3Field());
                  //  std::array<float, 3> posmomS = {0.,0.,0.};
                  // gpu::gpustd::array<float, 9> posmomStemp;
                  //  gpu::gpustd::array<float, 6> posmomS;
                  double posmomS[6], covS[21];
                  gpu::gpustd::array<float, 21> covStemp; // = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
                  SigmaTr.getCovXYZPxPyPzGlo(covStemp);
                  int chargeS = SigmaTr.getCharge();
                  SigmaTr.getXYZGlo(sigmaV);
                  // SigmaTr.getPxPyPzGlo(posmomS);
                  // SigmaTr.getPosDirGlo(posmomStemp);
                  for (int h = 0; h < 21; h++)
                  {
                    covS[h] = covStemp[h];
                  }
                  for (int h = 0; h < 6; h++)
                  {
                    if (h >= 3)
                    {
                      posmomS[h] = sigmaP[h - 3];
                    }
                    else
                      posmomS[h] = sigmaV[h];
                  }

                  kfSigma.Create(posmomS, covS, chargeS, 1.1974);

                  // double posmomP[6], covP[15];
                  // gpu::gpustd::array<float, 9> posmomPtemp;
                  // gpu::gpustd::array<float, 6> posmomP;
                  double posmomP[6], covP[21];
                  gpu::gpustd::array<float, 21> covPtemp; // = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
                  PionTr.getCovXYZPxPyPzGlo(covPtemp);
                  int chargeP = PionTr.getCharge();
                  // PionTr.GetXYZ(posmomP);
                  PionTr.getXYZGlo(pionV);
                  // PionTr.getPosDirGlo(posmomPtemp);
                  for (int h = 0; h < 21; h++)
                  {
                    covP[h] = covPtemp[h];
                  }
                  for (int h = 0; h < 6; h++)
                  {
                    if (h >= 3)
                    {
                      posmomP[h] = pionP[h - 3];
                    }
                    else
                      posmomP[h] = pionV[h];
                  }
                  kfPion.Create(posmomP, covP, chargeP, 0.1396);

                  // KFPVertex kfPVertex;
                  // kfPVertex.SetXYZ(SigmaTr.getX(), SigmaTr.getY(), SigmaTr.getZ());
                  // kfPVertex.SetCovarianceMatrix(covS[0], covS[1], covS[2], covS[3], covS[4], covS[5]);
                  // // // kfPVertex.SetChi2(fEventCut.GetPrimaryVertex()->GetChi2());
                  // // // kfPVertex.SetNDF(fEventCut.GetPrimaryVertex()->GetNDF());
                  // kfPVertex.SetNContributors(1);
                  // KFParticle prodVertex{kfPVertex};

                  // kfSigma.SetProductionVertex(prodVertex);

                  // KFParticle kfSigma0{kfSigma};
                  // KFParticle kfPion0{kfPion};
                  KFParticle kfNeutron;

                  kfNeutron.AddDaughter(kfSigma);
                  kfNeutron.NDF() = -1;
                  kfNeutron.Chi2() = 0.f;
                  kfNeutron.SubtractDaughter(kfPion);
                  kfNeutron.SetNonlinearMassConstraint(0.93957);
                  KFParticle Mother{kfNeutron, kfPion};

                  // kfNeutron.TransportToDecayVertex();
                  auto RResol = (sqrt((pionR[0] - Mother.X()) * (pionR[0] - Mother.X()) + (pionR[1] - Mother.Y()) * (pionR[1] - Mother.Y()))) / sqrt(pionR[0] * pionR[0] + pionR[1] * pionR[1]); //- sqrt(kfNeutron.X() * kfNeutron.X() + kfNeutron.Y() * kfNeutron.Y())) / sqrt(pionR[0] * pionR[0] + pionR[1] * pionR[1]);
                  // auto RResol1 = (sqrt(pionR[0] * pionR[0] + pionR[1] * pionR[1]) - sqrt(kfNeutron.X() * kfNeutron.X() + kfNeutron.Y() * kfNeutron.Y())) / sqrt(pionR[0] * pionR[0] + pionR[1] * pionR[1]);
                  //  auto RResol = sqrt(kfNeutron.X() * kfNeutron.X() + kfNeutron.Y() * kfNeutron.Y());
                  //   std::cout<< "--------------------------------"<<std::endl;
                  //   std::cout<<"After Neutron: " <<kfNeutron.GetP() << std::endl;
                  //   std::cout<<"Before Sigma " << sigmaPabs << std::endl;
                  //   std::cout<<"After Sigma " << kfSigma.GetP() << std::endl;
                  //   std::cout<<"Before Pion " << pionPabs << std::endl;
                  //   std::cout<<"After Pion " << kfPion.GetP() << std::endl;
                  //   if (Mother.GetChi2() < 0) //|| Mother.GetChi2() > 10 )
                  //     continue;
                  if (sqrt(kfNeutron.X() * kfNeutron.X() + kfNeutron.Y() * kfNeutron.Y()) < 17)
                    continue;
                  if (std::abs(etaS-etaP) > 0.3 || std::abs(phiS-phiP) > 0.3) continue;
                  if (std::abs(kfNeutron.GetMass()-0.939) > 0.1) continue;

                  h.Fill(Mother.GetChi2());
                  hmass.Fill(Mother.GetMass());
                  h1.Fill(RResol);
                  hmassResol.Fill(Mother.GetMass() - info[n][m].M);
                  h3.Fill(sqrt(kfPion.X() * kfPion.X() + kfPion.Y() * kfPion.Y()));
                  h6.Fill(RResol, sqrt(pionR[0] * pionR[0] + pionR[1] * pionR[1]));
                  h7.Fill(etaS - etaP, phiS - phiP);
                  h12.Fill(Mother.GetChi2(), Mother.GetMass() - info[n][m].M);
                  h10.Fill(Mother.GetP() - sigmaPabs);
                  h11.Fill(Mother.GetChi2(), Mother.GetP() - sigmaPabs);
                }
                if (FITTEROPTION == "DCA")
                {
                  try
                  {
                    DCAFitter2 ft2;
                    PionTr.checkCovariance();
                    SigmaTr.checkCovariance();
                    ft2.setUseAbsDCA(true);
                    ft2.setMaxChi2(5);
                    ft2.setBz(grp->getNominalL3Field());
                    ft2.process(PionTr, SigmaTr);
                    ft2.propagateTracksToVertex();
                    if (ft2.isPropagateTracksToVertexDone() == true)
                    {
                      auto SigmaTrDCA = ft2.getTrack(1);
                      auto PionTrDCA = ft2.getTrack(0);

                      SigmaTrDCA.getPxPyPzGlo(sigmaP);
                      sigmaPabs = SigmaTrDCA.getP();
                      etaS = SigmaTr.getEta();
                      phiS = SigmaTr.getPhi();

                      PionTrDCA.getPxPyPzGlo(pionP);
                      pionPabs = PionTrDCA.getP();
                      etaP = PionTr.getEta();
                      phiP = PionTr.getPhi();

                       if (ft2.getChi2AtPCACandidate() < 0)
                         continue;

                      std::array<float, 3> R = ft2.getPCACandidatePos();
                      auto RResol = (sqrt((pionR[0] - R[0]) * (pionR[0] - R[0]) + (pionR[1] - R[1]) * (pionR[1] - R[1]))) / sqrt(pionR[0] * pionR[0] + pionR[1] * pionR[1]);
                      if (sqrt(R[0] * R[0] + R[1] * R[1]) < 17)
                        continue;
                      if (std::abs(etaS-etaP) > 0.3 || std::abs(phiS-phiP) > 0.3) continue;
                      h.Fill(ft2.getChi2AtPCACandidate());
                      // std::cout<<ft2.getChi2AtPCACandidate()<<std::endl;
                      h1.Fill(RResol);
                      h7.Fill(etaS - etaP, phiS - phiP);
                      tempcount++;
                      h2.Fill(ft2.getChi2AtPCACandidate(), RResol); // ft2.getChi2AtPCACandidate()
                      h2R.Fill(sqrt(pionR[0] * pionR[0] + pionR[1] * pionR[1] + pionR[2] * pionR[2]), RResol);
                      float pionE = sqrt(0.1396 * 0.1396 + pionPabs * pionPabs);
                      // auto neutronE = sqrt(0.9396 * 0.9396 + pow((sigmaP[2] - pionP[2]), 2) + pow((sigmaP[1] - pionP[1]), 2) + pow((sigmaP[0] - pionP[0]), 2));
                      //  hpPion.Fill(pionP_gen);
                      //  hpPion_rec.Fill(pionPabs);
                      //  hpPionResol.Fill(pionP_gen - pionPabs);
                      // hmass.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - sigmaPabs * sigmaPabs));
                      float sigmaE = sqrt(1.1974 * 1.1974 + sigmaPabs * sigmaPabs); // 1.1974
                      float neutronPabs = sqrt(pow((sigmaP[2] - pionP[2]), 2) + pow((sigmaP[1] - pionP[1]), 2) + pow((sigmaP[0] - pionP[0]), 2));
                      float neutronM = sqrt((sigmaE - pionE) * (sigmaE - pionE) - neutronPabs * neutronPabs);
                      if (std::abs(neutronM-0.939) > 0.1) continue;
                      auto neutronE = sqrt(0.9396 * 0.9396 + pow((sigmaP[2] - pionP[2]), 2) + pow((sigmaP[1] - pionP[1]), 2) + pow((sigmaP[0] - pionP[0]), 2));
                      hmass.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - sigmaPabs * sigmaPabs));
                      hmassResol.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - sigmaPabs * sigmaPabs) - info[n][m].M);
                      h10.Fill(sigmaPabs - SigmaTrDCA.getP());
                      h11.Fill(ft2.getChi2AtPCACandidate(), sigmaPabs - SigmaTrDCA.getP());
                      h12.Fill(ft2.getChi2AtPCACandidate(), sqrt((neutronE + pionE) * (neutronE + pionE) - sigmaPabs * sigmaPabs) - info[n][m].M);
                    }
                  }
                  catch (std::runtime_error &e)
                  {
                    continue;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // TH1F h5 = h4 / h3;
  // h5.SetTitle("ITSTPC Efficiency");
  // h5.GetXaxis().("Pion in ITS hits");
  // h5.GetYaxis().("ITSTPC Efficiency");
  std::cout << tempcount << std::endl;
  TFile outputFile(("/home/justas_t/DCA_Fitter_Random_" + FITTEROPTION + "_s+b_cut" + ".root").c_str(), "recreate");
  // TFile outputFile("/home/justas_t/kfPlots_Resol.root", "recreate");

  h.Write();
  h1.Write();
  h2.Write();
  h2R.Write();
  h3.Write();
  h4.Write();
  hnResol.Write();
  // h5.Write();
  h6.Write();
  h7.Write();
  hmass.Write();
  hmassResol.Write();
  h10.Write();
  h11.Write();
  h12.Write();

  outputFile.Close();
}
