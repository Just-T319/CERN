#include <iostream>
#include "ROOT/RVec.hxx"

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

struct FilteredInfo
{
    int pdg;
    int pId;
    int evId;
    int motherId;
    float pt_mother;
    bool isRecoSigma = false;
    int isFakeSigma;
    int ITSTrackInd = -1;
    int TPCTrackInd = -1;
    int ITSTPCTrackInd = -1;
    int ITSTPCABTrackInd = -1;
    unsigned char isFakeITS = 0u;
    unsigned char isFakeTPC = 0u;
    unsigned char isFakeITSTPC = 0u;
    unsigned char isFakeITSTPCAB = 0u;
    float pt;
    float rapid;
    float rapid_mother;
    float Rmother;
};

void Filtering()
{

    ///////////Loading Data//////////////////////////////////
    auto fileName = "~/SigmaP/pion_detect_007.root";
    auto fileName_mc = "tf1/sgn_1_Kine.root";
    ROOT::RDataFrame parTree("PInfo", fileName);
    TChain mcTree("o2sim");
    mcTree.AddFile(fileName_mc);

    ///////////For Invariant Mass//////////////////////////////
    using namespace o2::its;
    using namespace o2::tpc;
    using namespace o2::dataformats;

    TChain TPCTree("tpcrec");
    TChain ITSTPCTree("matchTPCITS");
    TChain ITSTree("o2sim");

    auto fileName_its = "tf1/o2trac_its.root";
    auto fileName_tpc = "tf1/tpctracks.root";
    auto fileName_itstpc = "tf1/o2match_itstpc.root";
    ITSTree.AddFile(fileName_its);
    TPCTree.AddFile(fileName_tpc);
    ITSTPCTree.AddFile(fileName_itstpc);

    std::vector<TrackITS> *itsArr = nullptr;
    std::vector<o2::MCCompLabel> *itsLabArr = nullptr;
    std::vector<TrackTPC> *tpcArr = nullptr;
    std::vector<o2::MCCompLabel> *tpcLabArr = nullptr;
    std::vector<o2::MCCompLabel> *itstpcabLabArr = nullptr;
    std::vector<TrackTPCITS> *itstpcArr = nullptr;
    // std::vector<TrackTPCITS> *itstpcabArr = nullptr;
    std::vector<o2::MCCompLabel> *itstpcLabArr = nullptr;

    ITSTree.SetBranchAddress("ITSTrack", &itsArr);
    ITSTree.SetBranchAddress("ITSTrackMCTruth", &itsLabArr);
    TPCTree.SetBranchAddress("TPCTracks", &tpcArr);
    TPCTree.SetBranchAddress("TPCTracksMCTruth", &tpcLabArr);
    ITSTPCTree.SetBranchAddress("TPCITS", &itstpcArr);
    // ITSTPCTree.SetBranchAddress("TPCITSABRefs", &itstpcabArr);
    ITSTPCTree.SetBranchAddress("MatchMCTruth", &itstpcLabArr);
    ITSTPCTree.SetBranchAddress("MatchABMCTruth", &itstpcabLabArr);

    ITSTPCTree.GetEntry(0);
    TPCTree.GetEntry(0);
    ITSTree.GetEntry(0);

    ///////////////////////////////////////////////////////////

    //////////Some Preliminary Info///////////////////////////
    auto nsigma = parTree.Filter("std::abs(pdg) == 3222");
    auto npion = parTree.Filter("std::abs(pdg) == 211");
    auto nsigmaits = parTree.Filter("std::abs(pdg) == 3222 && ITSTrackInd >=0");
    auto npionits = parTree.Filter("std::abs(pdg) == 211 && ITSTrackInd >=0");
    auto nsigmatpc = parTree.Filter("std::abs(pdg) == 3222 && TPCTrackInd >=0");
    auto npiontpc = parTree.Filter("std::abs(pdg) == 211 && TPCTrackInd >=0");
    auto nsigmaitstpc = parTree.Filter("std::abs(pdg) == 3222 && (ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0)");
    auto npionitstpc = parTree.Filter("std::abs(pdg) == 211 && (ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0)");
    auto nsigmaf = parTree.Filter("std::abs(pdg) == 3222 && (isFakeITS >0 || isFakeTPC > 0 || isFakeITSTPC >0 || isFakeITSTPCAB >0)");
    auto nsigmarec = parTree.Filter("std::abs(pdg) == 3222 && (ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0 || TPCTrackInd >=0 || ITSTrackInd >=0)");
    auto npionrec = parTree.Filter("std::abs(pdg) == 211 && (ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0 || TPCTrackInd >=0 || ITSTrackInd >=0)");
    auto npionf = parTree.Filter("std::abs(pdg) == 211 && (isFakeITS >0 || isFakeTPC > 0 || isFakeITSTPC >0 || isFakeITSTPCAB >0)");

    // auto npiondoubles = parTree.Filter("std::abs(pdg) == 211 && (ITSTPCTrackInd >=0 && ITSTrackInd >=0)");
    // auto npiondoubles1 = parTree.Filter("std::abs(pdg) == 211 && (TPCTrackInd >=0 && ITSTrackInd >=0)");
    // auto npiondoubles2 = parTree.Filter("std::abs(pdg) == 211 && ((ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0) && TPCTrackInd >=0)");
    // auto npiondoubles3 = parTree.Filter("std::abs(pdg) == 211 && ((ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0)&& TPCTrackInd >=0 && ITSTrackInd >=0)");
    // auto nsigmadoubles = parTree.Filter("std::abs(pdg) == 3112 && (ITSTPCTrackInd >=0 && ITSTrackInd >=0)");
    // auto nsigmadoubles1 = parTree.Filter("std::abs(pdg) == 3112 && (TPCTrackInd >=0 && ITSTrackInd >=0)");
    // auto nsigmadoubles2 = parTree.Filter("std::abs(pdg) == 3112 && ((ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0) && TPCTrackInd >=0)");
    // auto nsigmadoubles3 = parTree.Filter("std::abs(pdg) == 3112 && ((ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0) && TPCTrackInd >=0 && ITSTrackInd >=0)");

    std::cout << "Total Sigma: " << *nsigma.Count() << std::endl;
    std::cout << "Total Pion: " << *npion.Count() << std::endl;
    std::cout << "Total Sigma reconstructed in ITS: " << *nsigmaits.Count() << std::endl;
    std::cout << "Total Pion reconstructed in ITS: " << *npionits.Count() << std::endl;
    std::cout << "Total Sigma reconstructed in TPC: " << *nsigmatpc.Count() << std::endl;
    std::cout << "Total Pion reconstructed in TPC: " << *npiontpc.Count() << std::endl;
    std::cout << "Total Sigma reconstructed in ITSTPC (inc AB): " << *nsigmaitstpc.Count() << std::endl;
    std::cout << "Total Pion reconstructed in ITSTPC (inc AB): " << *npionitstpc.Count() << std::endl;
    std::cout << "Total Sigma reconstructed: " << *nsigmarec.Count() << std::endl;
    std::cout << "Total Pion reconstructed: " << *npionrec.Count() << std::endl;
    std::cout << "Total Sigma fake: " << *nsigmaf.Count() << std::endl;
    std::cout << "Total Pion fake: " << *npionf.Count() << std::endl;

    // std::cout << "Pion doubles ITS and ITSPTC: " << *npiondoubles.Count() << std::endl;
    // std::cout << "Pion doubles ITS  and PTC: " << *npiondoubles1.Count() << std::endl;
    // std::cout << "Pion doubles TPC and ITSPTC: " << *npiondoubles2.Count() << std::endl;
    // std::cout << "Sigma doubles ITS and ITSPTC: " << *nsigmadoubles.Count() << std::endl;
    // std::cout << "Sigma doubles ITS  and PTC: " << *nsigmadoubles1.Count() << std::endl;
    // std::cout << "Sigma doubles TPC and ITSPTC: " << *nsigmadoubles2.Count() << std::endl;
    /////////////////////////////////////////////////////////////////////////////////

    /////////////////////Analysis////////////////////////////////////////////////////

    TH1F hfake("fake_info_sigma", "Fake clusters (wrt. Sigma)", 2, 0, 2);
    TH1F hfake1("fake_info_pion", "Fake clusters (wrt. Pion)", 2, 0, 2);
    TH1F hmass("mass_sigma", ";Invariant #Sigma^{-} M (Gev/#it{c^{2}}); Counts per 30 Mev/#it{c^{2}}", 50, 0.5, 2);
    TH1F hmassResol("mass_sigma_resol", ";Generated - Reconstructed #Sigma^{-} M (Gev/#it{c^{2}}); Counts per 50 Mev/#it{c^{2}}", 40, -1, 1);
    TH1F hmass_gen("mass_sigma_rec", "; Generated #Sigma^{-} M (Gev/#it{c^{2}}); Counts", 60, 0.8, 2);
    TH1F hpSigma("P_sigma", ";Generated #Sigma^{-} #p (GeV/#it{c}); Counts", 60, 0, 15);
    TH1F hpSigma_rec("P_sigma_rec", "; Reconstructed #Sigma^{-} #p (GeV/#it{c}); Counts", 60, 0, 15);
    TH1F hpSigmaResol("P_sigma_resol", "; Generated - Reconstructed #Sigma^{-} p (GeV/#it{c}); Counts", 120, -15, 15);
    TH1F hpPion("P_Pion", ";Generated #pi^{-} #p (GeV/#it{c}); Counts", 30, 0, 5);
    TH1F hpPion_rec("P_Pion_rec", "; Reconstructed #pi^{-} #p (GeV/#it{c}); Counts", 30, 0, 5);
    TH1F hpPionResol("P_Pion_resol", "; Generated - Reconstructed #pi^{-} p (GeV/#it{c}); Counts", 100, -3, 3);

    TChain pTree("PInfo");
    pTree.AddFile(fileName);
    int nev = 1000;
    PInfo *pArr = nullptr;
    std::vector<o2::MCTrack> *mcArr = nullptr;

    pTree.SetBranchAddress("particle", &pArr);
    mcTree.SetBranchAddress("MCTrack", &mcArr);

    int nparticles = pTree.GetEntries();
    std::vector<std::vector<int>> motherList(nev);
    std::vector<std::vector<float>> motherpt(nev);
    std::vector<std::vector<float>> motherrapid(nev);
    std::vector<std::vector<float>> motherfake(nev);
    std::vector<std::vector<float>> motherR(nev);
    std::vector<std::vector<float>> motherP(nev);
    std::vector<std::vector<float>> motherE(nev);
    std::vector<std::vector<float>> motherM(nev);
    std::vector<std::vector<float>> motherTempID(nev);
    std::vector<std::vector<std::array<float, 3>>> motherP_recV(nev);
    std::vector<std::vector<float>> motherP_rec(nev);
    std::vector<std::vector<int>> motherListTrackType(nev);
    std::vector<FilteredInfo> info(nparticles);

    auto totalPS = 0;
    auto totalPrecS = 0;
    auto totalPrecSrec = 0;
    auto totalPSrec = 0;

    for (int n = 0; n < nparticles; ++n)
    {
        pTree.GetEntry(n);

        int mId, evId;
        int tracktype = -1;
        int faketype = -1;
        if (std::abs(pArr->pdg) == 3222)
        {
            if (pArr->ITSTrackInd >= 0)
                tracktype = 0;
            if (pArr->TPCTrackInd >= 0)
                tracktype = 1;
            if (pArr->ITSTPCTrackInd >= 0)
                tracktype = 2;
            if (pArr->ITSTPCABTrackInd >= 0)
                tracktype = 3;

            if (pArr->isFakeITS > 0)
                faketype = 0;
            if (pArr->isFakeTPC > 0)
                faketype = 1;
            if (pArr->isFakeITSTPC > 0)
                faketype = 2;
            if (pArr->isFakeITSTPCAB > 0)
                faketype = 3;

            evId = pArr->evId;
            mcTree.GetEntry(evId);
            auto track = mcArr->at(pArr->pId);
            int neutron = 0;
            int pion = 0;
            int daughcount = track.getLastDaughterTrackId() - track.getFirstDaughterTrackId() + 1;

            for (int iDaug = track.getFirstDaughterTrackId(); iDaug <= track.getLastDaughterTrackId(); iDaug++)
            {
                if (iDaug < 0)
                    break;
                if (mcArr->at(iDaug).GetPdgCode() == 11)
                    daughcount--;

                if (std::abs(mcArr->at(iDaug).GetPdgCode()) == 211)
                    pion++;
                if (std::abs(mcArr->at(iDaug).GetPdgCode()) == 2112)
                    neutron++;
            }
            if (daughcount == 2 && neutron == 1 && pion == 1)
            {

                if (faketype >= 0 && tracktype >= 0)
                    hfake.Fill(0);

                motherList[evId].push_back(pArr->pId);
                motherListTrackType[evId].push_back(tracktype);
                motherpt[evId].push_back(pArr->pt);
                motherfake[evId].push_back(faketype);
                motherrapid[evId].push_back(pArr->rapid);
                auto P = pArr->Vmom;
                motherP[evId].push_back(sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]));
                auto R = pArr->Vcoor;
                motherR[evId].push_back(sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]));
                motherM[evId].push_back(pArr->M);
                motherE[evId].push_back(pArr->E);

                std::array<float, 3> sigmaP = {0, 0, 0};

                if (pArr->ITSTrackInd >= 0)
                {
                    motherP_rec[evId].push_back(itsArr->at(pArr->ITSTrackInd).getP());
                    itsArr->at(pArr->ITSTrackInd).getPxPyPzGlo(sigmaP);
                    motherTempID[evId].push_back(pArr->ITSTrackInd);
                }
                else if (pArr->TPCTrackInd >= 0)
                {
                    motherP_rec[evId].push_back(tpcArr->at(pArr->TPCTrackInd).getP());
                    tpcArr->at(pArr->TPCTrackInd).getPxPyPzGlo(sigmaP);
                    motherTempID[evId].push_back(pArr->TPCTrackInd);
                }
                else if (pArr->ITSTPCTrackInd >= 0)
                {
                    motherP_rec[evId].push_back(itstpcArr->at(pArr->ITSTPCTrackInd).getP());
                    itstpcArr->at(pArr->ITSTPCTrackInd).getPxPyPzGlo(sigmaP);
                    motherTempID[evId].push_back(pArr->ITSTPCTrackInd);
                }
                else if (pArr->ITSTPCABTrackInd >= 0)
                {
                    motherP_rec[evId].push_back(itstpcArr->at(pArr->ITSTPCABTrackInd).getP());
                    itstpcArr->at(pArr->ITSTPCABTrackInd).getPxPyPzGlo(sigmaP);
                    motherTempID[evId].push_back(pArr->ITSTPCABTrackInd);
                }
                else
                {
                    motherP_rec[evId].push_back(-1);
                    motherTempID[evId].push_back(-1);
                } // So that index are always aligned
                motherP_recV[evId].push_back(sigmaP);
            }
            else if (faketype >= 0 && tracktype >= 0)
                hfake.Fill(1);
        }
    }
    int counter = 0;
    int counterrec = 0;

    for (int n = 0; n < nev; ++n)
    {
        counter += motherList[n].size();
        for (int k = 0; k < motherListTrackType[n].size(); ++k)
        {
            if (motherListTrackType[n][k] >= 0)
            {
                counterrec++;
            }
        }
    }
    std::cout << "Sigmas resulting in correct decay:  " << counter << std::endl;
    std::cout << "Reconstructed Sigmas resulting in correct decay:  " << counterrec << std::endl;
    for (int n = 0; n < nparticles; ++n)
    {
        pTree.GetEntry(n);
        int mId, evId;
        if (std::abs(pArr->pdg) == 211)
        {
            mId = pArr->motherId;
            evId = pArr->evId;
            int faketype = -1;
            auto p = std::find(std::begin(motherList[evId]), std::end(motherList[evId]), mId);
            int index = std::distance(motherList[evId].begin(), p);
            if (p != std::end(motherList[evId]))
            {
                totalPS++;
                bool Reco = false;
                bool Fake = false;
                auto pionR = pArr->Vcoor;

                info[n].pdg = pArr->pdg;
                info[n].pt = pArr->pt;
                info[n].pId = pArr->pId;
                info[n].evId = pArr->evId;
                info[n].motherId = pArr->motherId;
                info[n].isFakeITS = pArr->isFakeITS;
                info[n].isFakeTPC = pArr->isFakeTPC;
                info[n].isFakeITSTPC = pArr->isFakeITSTPC;
                info[n].isFakeITSTPCAB = pArr->isFakeITSTPCAB;
                info[n].rapid = pArr->rapid;
                info[n].ITSTrackInd = pArr->ITSTrackInd;
                info[n].TPCTrackInd = pArr->TPCTrackInd;
                info[n].ITSTPCTrackInd = pArr->ITSTPCTrackInd;
                info[n].ITSTPCABTrackInd = pArr->ITSTPCABTrackInd;
                info[n].pt_mother = motherpt[evId][index];
                info[n].rapid_mother = motherrapid[evId][index];
                info[n].Rmother = (sqrt(pionR[0] * pionR[0] + pionR[1] * pionR[1] + pionR[2] * pionR[2]) - motherR[evId][index]);

                if (pArr->ITSTrackInd >= 0 || pArr->TPCTrackInd >= 0 || pArr->ITSTPCTrackInd >= 0 || pArr->ITSTPCABTrackInd >= 0)
                {

                    if (pArr->isFakeITS > 0 || pArr->isFakeTPC > 0 || pArr->isFakeITSTPC > 0 || pArr->isFakeITSTPCAB > 0)
                        hfake1.Fill(0);

                    totalPrecS++;
                    if (motherListTrackType[evId][index] == 0 || motherListTrackType[evId][index] == 1 || motherListTrackType[evId][index] == 2 || motherListTrackType[evId][index] == 3)
                    {

                        if ((pArr->isFakeITS <= 0 || pArr->isFakeTPC <= 0 || pArr->isFakeITSTPC <= 0 || pArr->isFakeITSTPCAB <= 0) && motherfake[evId][index] < 0)
                        {
                            hpSigma.Fill(motherP[evId][index]);
                            hpSigma_rec.Fill(motherP_rec[evId][index]);
                            hpSigmaResol.Fill(motherP[evId][index] - motherP_rec[evId][index]);
                            hmass_gen.Fill(motherM[evId][index]);
                            auto P = pArr->Vmom;
                            auto pionP_gen = sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
                            // std::vector<o2::track::TrackParametrization> pionP;
                            std::array<float, 3> pionP = {0, 0, 0};
                            float pionPabs = 0;

                            if (pArr->ITSTrackInd >= 0 && pArr->isFakeITS <= 0)
                            {
                                pionPabs = itsArr->at(pArr->ITSTrackInd).getP();
                                itsArr->at(pArr->ITSTrackInd).getPxPyPzGlo(pionP);
                                float pionE = sqrt(0.1396 * 0.1396 + pionPabs * pionPabs);
                                auto neutronE = sqrt(0.9396 * 0.9396 + pow((motherP_recV[evId][index][2] - pionP[2]), 2) + pow((motherP_recV[evId][index][1] - pionP[1]), 2) + pow((motherP_recV[evId][index][0] - pionP[0]), 2));
                                hpPion.Fill(pionP_gen);
                                hpPion_rec.Fill(pionPabs);
                                hpPionResol.Fill(pionP_gen - pionPabs);
                                hmass.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - motherP_rec[evId][index] * motherP_rec[evId][index]));
                                hmassResol.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - motherP_rec[evId][index] * motherP_rec[evId][index]) - motherM[evId][index]);
                            }
                            else if (pArr->TPCTrackInd >= 0 && pArr->isFakeTPC <= 0)
                            {
                                pionPabs = tpcArr->at(pArr->TPCTrackInd).getP();
                                tpcArr->at(pArr->TPCTrackInd).getPxPyPzGlo(pionP);
                                float pionE = sqrt(0.1396 * 0.1396 + pionPabs * pionPabs);
                                auto neutronE = sqrt(0.9396 * 0.9396 + pow((motherP_recV[evId][index][2] - pionP[2]), 2) + pow((motherP_recV[evId][index][1] - pionP[1]), 2) + pow((motherP_recV[evId][index][0] - pionP[0]), 2));
                                hpPion.Fill(pionP_gen);
                                hpPion_rec.Fill(pionPabs);
                                hpPionResol.Fill(pionP_gen - pionPabs);
                                hmass.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - motherP_rec[evId][index] * motherP_rec[evId][index]));
                                hmassResol.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - motherP_rec[evId][index] * motherP_rec[evId][index]) - motherM[evId][index]);
                            }

                            else if (pArr->ITSTPCTrackInd >= 0 && pArr->isFakeITSTPC <= 0)
                            {
                                pionPabs = itstpcArr->at(pArr->ITSTPCTrackInd).getP();
                                itstpcArr->at(pArr->ITSTPCTrackInd).getPxPyPzGlo(pionP);
                                float pionE = sqrt(0.1396 * 0.1396 + pionPabs * pionPabs);
                                auto neutronE = sqrt(0.9396 * 0.9396 + pow((motherP_recV[evId][index][2] - pionP[2]), 2) + pow((motherP_recV[evId][index][1] - pionP[1]), 2) + pow((motherP_recV[evId][index][0] - pionP[0]), 2));
                                hpPion.Fill(pionP_gen);
                                hpPion_rec.Fill(pionPabs);
                                hpPionResol.Fill(pionP_gen - pionPabs);
                                hmass.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - motherP_rec[evId][index] * motherP_rec[evId][index]));
                                hmassResol.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - motherP_rec[evId][index] * motherP_rec[evId][index]) - motherM[evId][index]);
                            }
                            else if (pArr->ITSTPCABTrackInd >= 0 && pArr->isFakeITSTPCAB <= 0)
                            {
                                pionPabs = itstpcArr->at(pArr->ITSTPCABTrackInd).getP();
                                itstpcArr->at(pArr->ITSTPCABTrackInd).getPxPyPzGlo(pionP);
                                float pionE = sqrt(0.1396 * 0.1396 + pionPabs * pionPabs);
                                auto neutronE = sqrt(0.9396 * 0.9396 + pow((motherP_recV[evId][index][2] - pionP[2]), 2) + pow((motherP_recV[evId][index][1] - pionP[1]), 2) + pow((motherP_recV[evId][index][0] - pionP[0]), 2));
                                hpPion.Fill(pionP_gen);
                                hpPion_rec.Fill(pionPabs);
                                hpPionResol.Fill(pionP_gen - pionPabs);
                                hmass.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - motherP_rec[evId][index] * motherP_rec[evId][index]));
                                hmassResol.Fill(sqrt((neutronE + pionE) * (neutronE + pionE) - motherP_rec[evId][index] * motherP_rec[evId][index]) - motherM[evId][index]);
                            }

                            // if (pArr->ITSTrackInd >= 0) {auto lab0 = itsArr->at(pArr->ITSTrackInd);
                            //}

                            // else if (pArr->TPCTrackInd >= 0) {auto lab1 = tpcArr->at(pArr->TPCTrackInd);
                            // hmass.Fill(sqrt(motherE[evId][index]*motherE[evId][index] - lab1.getP()*lab1.getP()));}
                            // else if (pArr->ITSTPCTrackInd >= 0) {auto lab2 = itstpcArr->at(pArr->ITSTPCTrackInd);
                            // hmass.Fill(sqrt(motherE[evId][index]*motherE[evId][index] - lab2.getP()*lab2.getP()));}
                            // else if (pArr->ITSTPCABTrackInd >= 0){ auto lab3 = itstpcArr->at(pArr->ITSTPCABTrackInd);
                            // hmass.Fill(sqrt(motherE[evId][index]*motherE[evId][index] - lab3.getP()*lab3.getP()));}
                        }
                        totalPrecSrec++;
                    }
                }

                if (motherListTrackType[evId][index] == 0 || motherListTrackType[evId][index] == 1 || motherListTrackType[evId][index] == 2 || motherListTrackType[evId][index] == 3)
                {
                    totalPSrec++;
                    Reco = true;
                }

                info[n].isRecoSigma = Reco;
                info[n].isFakeSigma = motherfake[evId][index];
            }

            else
            {
                if (pArr->isFakeITS > 0 || pArr->isFakeTPC > 0 || pArr->isFakeITSTPC > 0 || pArr->isFakeITSTPCAB > 0)
                    hfake1.Fill(1);
            }
        }
    }

    std::cout << "Total Pions from Sigma decay: " << totalPS << std::endl;
    std::cout << "Total reconstructed Pions from Sigma decay: " << totalPrecS << std::endl;
    std::cout << "Total Pions from reconstructed Sigma decay: " << totalPSrec << std::endl;
    std::cout << "Total reconstructed Pions from reconstructed Sigma decay: " << totalPrecSrec << std::endl;

    ////////////////////////////////////////////////////////////

    hfake.GetXaxis()->SetBinLabel(1, "Pion");
    hfake.GetXaxis()->SetBinLabel(2, "Not Pion");
    hfake1.GetXaxis()->SetBinLabel(1, "Sigma");
    hfake1.GetXaxis()->SetBinLabel(2, "Not Sigma");

    //////////////////Saving Filtered Data///////////////////

    TFile file("/home/justas_t/filtered_pion_007.root", "recreate");
    TTree tree("FilteredInfo", "FilteredInfo");
    FilteredInfo pInfo;
    tree.Branch("particle", &pInfo);
    for (auto &part : info)
    {

        pInfo = part;
        tree.Fill();
    }

    tree.Write();
    hfake.Write();
    hfake1.Write();
    hpSigma.Write();
    hpSigma_rec.Write();
    hpSigmaResol.Write();
    hpPion.Write();
    hpPion_rec.Write();
    hpPionResol.Write();
    hmass.Write();
    hmassResol.Write();
    hmass_gen.Write();
    file.Close();

    // for (int i =4; i<8;++i){

    // auto nnsigma = parTree.Filter(Form("std::abs(pdg) == 3112 && clusters == (pow(2,%i)-1)", i) ).Count();
    // auto crsigma = parTree.Filter(Form("std::abs(pdg) == 3112 && clusters == (pow(2,%i)-1) && isReco > 0 ", i) ).Count();
    // auto recofakesigma = parTree.Filter(Form("std::abs(pdg) == 3112 && clusters == (pow(2,%i)-1) && isFake > 0 ", i)).Count();
    // float temp1 = *crsigma ;

    // float temp2 = *recofakesigma ;
    // float temp3 = *nnsigma;

    //  std::cout<< "Total Sigma that has "<<i<<" hits: "<<*nnsigma<<std::endl;
    //  std::cout<< "Total Reco tracks that has (t+f) "<<i<<" hits: "<<*crsigma<<std::endl;
    // if (temp3!=0){

    // std::cout<< "Efficiency: " << (temp1+temp2)/temp3*100 <<"%"<<std::endl;
    // }
    // for (int q = 0; q<temp3; q++){
    // h1.Fill(i);
    // }

    // for (int q = 0; q<(temp1+temp2); q++){
    // h2.Fill(i);
    // }

    // //gr->AddPoint(i,(temp0+temp1)/temp3*100);

    // }

    // //gr->Draw("AL*");
    // h3 = h2/h1;
    // TFile outputFile("/home/justas_t/efficiency.root", "recreate");

    // //gr->Write();
    // h1.Write();
    // h2.Write();
    // h3.Write();

    // outputFile.Close();
}