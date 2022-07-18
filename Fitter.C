
// step0 : set the DCAFitter : cut on chi2<5
// step 1 : take ITS tracks
// step 2 : take ITS-TPC tracks
// step3 : Cut on the minimum V0radius : the reconstructed one has to be> 17(2 cm of tolerance)
// step3 .1 : Cut on the chi2 < 0
// step4 : Cut on invariant mass
// step5 : --->check the results

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "TGeoGlobalMagField.h"

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
    int pdg;
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
};

void Fitter()
{
    string FITTEROPTION = "DCA"; // "DCA" or "KFParticle"

    using namespace o2;
    using namespace vertexing;
    using namespace o2::its;
    using namespace o2::tpc;
    using namespace o2::dataformats;
    using GIndex = o2::dataformats::VtxTrackIndex;
    using namespace o2::itsmft;

    const auto grp = o2::parameters::GRPObject::loadFrom("/data/hyper/007/tf2/o2sim_grp.root");
    o2::base::GeometryManager::loadGeometry("tf2/");
    auto gman = o2::its::GeometryTGeo::Instance();

    TChain TPCTree("tpcrec");
    TChain ITSTPCTree("matchTPCITS");
    TChain ITSTree("o2sim");
    TChain pTree("PInfo");
    TChain mcTree("o2sim");
    TChain pvTree("o2sim");

    auto fileName_pv = "tf2/o2_primary_vertex.root";
    auto fileName_mc = "tf2/sgn_2_Kine.root";
    auto fileName_its = "tf2/o2trac_its.root";
    auto fileName_tpc = "tf2/tpctracks.root";
    auto fileName_itstpc = "tf2/o2match_itstpc.root";
    ITSTree.AddFile(fileName_its);
    TPCTree.AddFile(fileName_tpc);
    ITSTPCTree.AddFile(fileName_itstpc);
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

    TH1F h("DCA", "; #Chi^{2} ; Counts", 8, -2, 6);
    TH1F hmass("mass_sigma", ";Invariant #Sigma^{-} M (Gev/#it{c^{2}}); Counts per 8 Mev/#it{c^{2}}", 25, 1.1, 1.3);

    for (int n = 0; n < 1000; ++n)
    {
        mcTree.GetEntry(n);
        info[n].resize(mcArr->size());
        for (unsigned int m = 0; m < mcArr->size(); ++m)
        {
            auto part = mcArr->at(m);
            info[n][m].pdg = part.GetPdgCode();
        }
    }

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
        if (std::abs(info[evID][trackID].pdg) == 3112)
        {
            info[evID][trackID].isFakeITS += fake;
            info[evID][trackID].ITSTrackInd = n;
        }
    }
    std::cout << "ITS tracks done." << endl;

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

                continue;
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

                continue;
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

                continue;
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

                continue;
            }
        }
    }

    std::cout << "Pions from PV: " << counterPion << std::endl;
    std::cout << "PV with different ID: " << diffID << std::endl; // tvid.getIndex()
    std::cout << "PV with same ID: " << sameID << std::endl;      // tvid.getIndex()
    std::cout << "Primary Vertex done..." << std::endl;

    for (int n = 0; n < 1000; ++n)
    {
        mcTree.GetEntry(n);
        for (unsigned int m = 0; m < mcArr->size(); ++m)
        {
            if (std::abs(info[n][m].pdg) == 3112 && info[n][m].ITSTrackInd >= 0 && info[n][m].isFakeITS <= 0)
            {
                for (unsigned int k = 0; k < mcArr->size(); ++k)
                {
                    if (std::abs(info[n][k].pdg) == 211 && ((info[n][k].ITSTPCABTrackInd >= 0 && info[n][k].isFakeITSTPCAB <= 0 && info[n][k].ITSTPCPVInd < 0) || (info[n][k].ITSTPCTrackInd >= 0 && info[n][k].isFakeITSTPC <= 0 && info[n][k].ITSTPCPVInd < 0)))
                    {
                        if (FITTEROPTION == "DCA"){
                        o2::track::TrackParCov SigmaTr;
                        o2::track::TrackParCov PionTr;
                        DCAFitter2 ft2;
                        ft2.setBz(grp->getNominalL3Field());
                        ft2.setMaxChi2(5);
                        ft2.setUseAbsDCA(true);

                        std::array<float, 3> sigmaP = {0, 0, 0};
                        std::array<float, 3> pionP = {0, 0, 0};
                        float pionPabs = 0;
                        float sigmaPabs = 0;

                        float etaS = 0;
                        float phiS = 0;
                        float etaP = 0;
                        float phiP = 0;

                        SigmaTr = itsArr->at(info[n][m].ITSTrackInd);//.getParamOut();//=
                        itsArr->at(info[n][m].ITSTrackInd).getPxPyPzGlo(sigmaP);
                        sigmaPabs = itsArr->at(info[n][m].ITSTrackInd).getP();
                        etaS = itsArr->at(info[n][m].ITSTrackInd).getEta();
                        phiS = itsArr->at(info[n][m].ITSTrackInd).getPhi();

                        if ((info[n][k].ITSTPCTrackInd >= 0 && info[n][k].isFakeITSTPC <= 0 && info[n][k].ITSTPCPVInd < 0))
                        {
                            PionTr = itstpcArr->at(info[n][k].ITSTPCTrackInd);//.getParamOut();
                            pionPabs = itstpcArr->at(info[n][k].ITSTPCTrackInd).getP();
                            itstpcArr->at(info[n][k].ITSTPCTrackInd).getPxPyPzGlo(pionP);
                            etaP = itstpcArr->at(info[n][k].ITSTPCTrackInd).getEta();
                            phiP = itstpcArr->at(info[n][k].ITSTPCTrackInd).getPhi();
                        }
                        else if ((info[n][k].ITSTPCABTrackInd >= 0 && info[n][k].isFakeITSTPCAB <= 0 && info[n][k].ITSTPCPVInd < 0))
                        {
                            PionTr = itstpcArr->at(info[n][k].ITSTPCABTrackInd);//.getParamOut();
                            pionPabs = itstpcArr->at(info[n][k].ITSTPCABTrackInd).getP();
                            itstpcArr->at(info[n][k].ITSTPCABTrackInd).getPxPyPzGlo(pionP);
                            etaP = itstpcArr->at(info[n][k].ITSTPCABTrackInd).getEta();
                            phiP = itstpcArr->at(info[n][k].ITSTPCABTrackInd).getPhi();
                        }

                        PionTr.checkCovariance();
                        SigmaTr.checkCovariance();
                        std::array<float, 3> R = ft2.getPCACandidatePos();

                        ft2.process(PionTr, SigmaTr);
                        //if (ft2.getChi2AtPCACandidate() <= 0 || ft2.getChi2AtPCACandidate() > 3) // cut on chi2
                            //continue;

                        //if (sqrt(R[0] * R[0] + R[1] * R[1]) < 17) // cut on radius
                            //continue;

                        //if (std::abs(etaS - etaP) > 0.3 || std::abs(phiS - phiP) > 0.3)
                            //continue;

                        h.Fill(ft2.getChi2AtPCACandidate());
                        tempcount++;
                        float pionE = sqrt(0.13957 * 0.13957 + pionPabs * pionPabs);
                        float sigmaE = sqrt(1.1974 * 1.1974 + sigmaPabs * sigmaPabs); // 1.1974
                        float neutronPabs = sqrt(pow((sigmaP[2] - pionP[2]), 2) + pow((sigmaP[1] - pionP[1]), 2) + pow((sigmaP[0] - pionP[0]), 2));
                        float neutronM = sqrt((sigmaE - pionE) * (sigmaE - pionE) - neutronPabs * neutronPabs);
                        // if (neutronM > 0.990 || neutronM < 0.89 ) continue;
                        auto neutronE = sqrt(0.93957 * 0.93957 + pow((sigmaP[2] - pionP[2]), 2) + pow((sigmaP[1] - pionP[1]), 2) + pow((sigmaP[0] - pionP[0]), 2));
                        float sigmaM = sqrt((neutronE + pionE) * (neutronE + pionE) - sigmaPabs * sigmaPabs);
                        // hmass.Fill(neutronM);//sigma and pion hypothesis
                        if (sigmaM > 1.3 || sigmaM < 1.1)
                            continue;
                        hmass.Fill(sigmaM); // pion and neutron hypothesis
                        }
                        if (FITTEROPTION == "KFParticle"){

                        }
                    }
                }
            }
        }
    }

    std::cout << tempcount << std::endl;
    TFile outputFile("/home/justas_t/SigmaP/DCA_Fitter_SigmaP.root", "recreate");

    h.Write();

    hmass.Write();

    outputFile.Close();
}
