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

void EffPion()
{

    ROOT::RDataFrame par0Tree("FilteredInfo", "~/filtered_pion_007.root");
    auto parTree = par0Tree.Filter("std::abs(rapid_mother) < 1");
    auto parTreeSigma = par0Tree.Filter("std::abs(rapid_mother) < 1");//for sigma rapidity

    std::cout << "Total Pions from Sigmas: " << *parTree.Filter("std::abs(pdg) == 211").Count() << std::endl;

    auto pt_gen_temp = parTree.Filter("std::abs(pdg) == 211");
    auto pt_gen = pt_gen_temp.Take<float>("pt");

    auto pt_gen_temp_sigma = parTreeSigma.Filter("std::abs(pdg) == 211");
    auto pt_gen_mother = pt_gen_temp_sigma.Take<float>("pt_mother");
    auto R_gen_mother = pt_gen_temp_sigma.Take<float>("Rmother");

    TH1F pth("Eff_pt_gen", "Gener. pt Distribution; GeV/c", 40, 0, 5);

     for (auto &myItem : pt_gen)
    {
        pth.Fill(myItem);
    }

    TH1F pth_mother("Eff_pt_gen_mother", "Gen. Pion as a fn of mother pt; GeV/c", 60, 0, 15);

     for (auto &myItem : pt_gen_mother)
    {
        pth_mother.Fill(myItem);
    }

    TH1F r_mother("Eff_R_gen_mother", "Gen. Pion as a fn of mother R; cm", 120, 0, 60);

     for (auto &myItem : R_gen_mother)
    {
        r_mother.Fill(myItem);
    }


    ////////////////ITS DETECTOR///////////////////////////////////
    TH1F ITSh1("ITS_Eff_pt_rec", "Recons. pt Distribution; GeV/c", 40, 0, 5);
    TH1F ITSh2("ITS_pt_notFake", "Not Fake pt Distribution; GeV/c", 40, 0, 5);
    TH1F ITSh5("ITS_pt_decay_Eff", "Topology Eff; GeV/c", 60, 0, 15);
    TH1F ITSh7("ITS_pt_decay_Eff_NotFake", "Topology Eff (Not Fake); GeV/c", 60, 0, 15);
    TH1F ITSh10("ITS_R_decay_Eff_Radius_NotF", "Topology Eff (Not Fake) (Sigma R); cm", 120, 0, 60);
    TH1F ITSh11("ITS_R_decay_Eff_Radius", "Topology Eff (Sigma R); cm", 120, 0, 60);


    auto ITS_pt_rec_temp = parTree.Filter("std::abs(pdg) == 211 && ITSTrackInd >=0");
    auto ITS_pt_rec = ITS_pt_rec_temp.Take<float>("pt");
    auto ITSpurNotf = parTree.Filter("std::abs(pdg) == 211 && (ITSTrackInd >=0 && isFakeITS <= 0)");
    auto ITSptNotf = ITSpurNotf.Take<float>("pt");

    auto ITS_pt_decay_temp = parTreeSigma.Filter("std::abs(pdg) == 211 && ITSTrackInd >=0 && isRecoSigma == true");
    auto ITS_pt_decay = ITS_pt_decay_temp.Take<float>("pt_mother");

    auto ITS_pt_decay_tempNotFake = parTreeSigma.Filter("std::abs(pdg) == 211 && ITSTrackInd >=0 && isRecoSigma == true && isFakeSigma < 0 && isFakeITS <= 0 ");
    auto ITS_pt_decay_NotFake = ITS_pt_decay_tempNotFake.Take<float>("pt_mother");

    auto ITS_R_decay_tempNotFake = parTreeSigma.Filter("std::abs(pdg) == 211 && ITSTrackInd >=0 && isRecoSigma == true && isFakeSigma < 0 && isFakeITS <= 0 ");
    auto ITS_R_decay_NotFake = ITS_R_decay_tempNotFake.Take<float>("Rmother");
    auto ITS_R_decay_temp = parTreeSigma.Filter("std::abs(pdg) == 211 && ITSTrackInd >=0 && isRecoSigma == true");
    auto ITS_R_decay = ITS_R_decay_temp.Take<float>("Rmother");

    for (auto &myItem : ITS_R_decay_NotFake)
    {
        ITSh10.Fill(myItem);
    }

    for (auto &myItem : ITS_R_decay)
    {
        ITSh11.Fill(myItem);
    }

    for (auto &myItem : ITS_pt_decay_NotFake)
    {
        ITSh7.Fill(myItem);
    }

    for (auto &myItem : ITS_pt_decay)
    {
        ITSh5.Fill(myItem);
    }

    for (auto &myItem : ITSptNotf)
    {
        ITSh2.Fill(myItem);
    }
    for (auto &myItem : ITS_pt_rec)
    {
        ITSh1.Fill(myItem);
    }

    TH1F ITSh3 = ITSh1 / pth;
    TH1F ITSh4 = ITSh2 / ITSh1;

    TH1F ITSh6 = ITSh5 / pth_mother;
    TH1F ITSh8 = ITSh7 / pth_mother;
    TH1F ITSh9 = ITSh7 / ITSh5;

    TH1F ITSh12 = ITSh11 / r_mother;
    TH1F ITSh13 = ITSh10 / ITSh11;
    TH1F ITSh14 = ITSh10 / r_mother;
    ITSh14.SetTitle("Efficiency of topology reconstruction (sigma R) (Not Fake) (ITS)");
    ITSh13.SetTitle("Total Topology Purity (sigma R) (ITS)");
    ITSh12.SetTitle("Efficiency of Topology reconstruction (sigma R) (ITS)");

    ITSh9.SetTitle("Total Topology Purity (ITS)");
    ITSh6.SetTitle("Efficiency of Topology reconstruction (sigma pt) (ITS)");
    ITSh8.SetTitle("Efficiency of topology reconstruction (sigma pt) (Not Fake) (ITS)");

    ITSh3.SetTitle("Efficiency of Pion reconstruction (ITS)");
    ITSh4.SetTitle("Purity of Pion Track (ITS)");

    ////////////////TPC DETECTOR///////////////////////////////////
    TH1F TPCh1("TPC_Eff_pt_rec", "Recons. pt Distribution; GeV/c", 40, 0, 5);
    TH1F TPCh2("TPC_pt_notFake", "Not Fake pt Distribution; GeV/c", 40, 0, 5);
    TH1F TPCh5("TPC_pt_decay_Eff", "Topology Eff (TPC); GeV/c", 60, 0, 15);
    TH1F TPCh7("TPC_pt_decay_Eff_NotFake", "Topology Eff (Not Fake) (TPC); GeV/c", 60, 0, 15);
    TH1F TPCh10("TPC_R_decay_Eff_Radius_NotF", "Topology Eff (Not Fake) (Sigma R) ; cm", 120, 0, 60);
    TH1F TPCh11("TPC_R_decay_Eff_Radius", "Topology Eff (Sigma R); cm", 120, 0, 60);

    auto TPC_pt_rec_temp = parTree.Filter("std::abs(pdg) == 211 && TPCTrackInd >=0");
    auto TPC_pt_rec = TPC_pt_rec_temp.Take<float>("pt");
    auto TPCpurNotf = parTree.Filter("std::abs(pdg) == 211 && (TPCTrackInd >=0 && isFakeTPC <= 0)");
    auto TPCptNotf = TPCpurNotf.Take<float>("pt");

    auto TPC_pt_decay_temp = parTreeSigma.Filter("std::abs(pdg) == 211 && TPCTrackInd >=0 && isRecoSigma == true");
    auto TPC_pt_decay = TPC_pt_decay_temp.Take<float>("pt_mother");

    auto TPC_pt_decay_tempNotFake = parTreeSigma.Filter("std::abs(pdg) == 211 && TPCTrackInd >=0 && isRecoSigma == true && isFakeSigma < 0 && isFakeTPC <= 0 ");
    auto TPC_pt_decay_NotFake = TPC_pt_decay_tempNotFake.Take<float>("pt_mother");

    auto TPC_R_decay_tempNotFake = parTreeSigma.Filter("std::abs(pdg) == 211 && TPCTrackInd >=0 && isRecoSigma == true && isFakeSigma < 0 && isFakeTPC <= 0 ");
    auto TPC_R_decay_NotFake = TPC_R_decay_tempNotFake.Take<float>("Rmother");
    auto TPC_R_decay_temp = parTreeSigma.Filter("std::abs(pdg) == 211 && TPCTrackInd >=0 && isRecoSigma == true");
    auto TPC_R_decay = TPC_R_decay_temp.Take<float>("Rmother");

    for (auto &myItem : TPC_R_decay_NotFake)
    {
        TPCh10.Fill(myItem);
    }

    for (auto &myItem : TPC_R_decay)
    {
        TPCh11.Fill(myItem);
    }

    for (auto &myItem : TPC_pt_decay_NotFake)
    {
        TPCh7.Fill(myItem);
    }

    for (auto &myItem : TPC_pt_decay)
    {
        TPCh5.Fill(myItem);
    }

    for (auto &myItem : TPCptNotf)
    {
        TPCh2.Fill(myItem);
    }
    for (auto &myItem : TPC_pt_rec)
    {
        TPCh1.Fill(myItem);
    }

    TH1F TPCh3 = TPCh1 / pth;
    TH1F TPCh4 = TPCh2 / TPCh1;

    TH1F TPCh6 = TPCh5 / pth_mother;
    TH1F TPCh8 = TPCh7 / pth_mother;
    TH1F TPCh9 = TPCh7 / TPCh5;

    TH1F TPCh12 = TPCh11 / r_mother;
    TH1F TPCh13 = TPCh10 / TPCh11;
    TH1F TPCh14 = TPCh10 / r_mother;
    TPCh14.SetTitle("Efficiency of topology reconstruction (sigma R) (Not Fake) (TPC)");
    TPCh13.SetTitle("Total Topology Purity (sigma R) (TPC)");
    TPCh12.SetTitle("Efficiency of Topology reconstruction (sigma R) (TPC)");

    TPCh9.SetTitle("Total Topology Purity (TPC)");
    TPCh6.SetTitle("Efficiency of topology reconstruction (sigma pt) (TPC)");
    TPCh8.SetTitle("Efficiency of topology reconstruction (sigma pt) (Not Fake) (TPC)");

    TPCh3.SetTitle("Efficiency of Pion reconstruction (TPC)");
    TPCh4.SetTitle("Purity of Pion Track (TPC)");

    ////////////////ITSTPC DETECTOR///////////////////////////////////
    TH1F ITSTPCh1("ITSTPC_Eff_pt_rec", "Recons. pt Distribution; GeV/c", 40, 0, 5);
    TH1F ITSTPCh2("ITSTPC_pt_notFake", "Not Fake pt Distribution; GeV/c", 40, 0, 5);
    TH1F ITSTPCh5("ITSTPC_pt_decay_Eff", "Topology Eff (ITSTPC); GeV/c", 60, 0, 15);
    TH1F ITSTPCh7("ITSTPC_pt_decay_Eff_NotFake", "Topology Eff (Not Fake) (ITSTPC); GeV/c", 60, 0, 15);
    TH1F ITSTPCh10("ITSTPC_R_decay_Eff_Radius_NotF", "Topology Eff (Not Fake) (Sigma R) ; cm", 120, 0, 60);
    TH1F ITSTPCh11("ITSTPC_R_decay_Eff_Radius", "Topology Eff (Sigma R); cm", 120, 0, 60);

    auto ITSTPC_pt_rec_temp = parTree.Filter("std::abs(pdg) == 211 && (ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0)");
    auto tempcount =  parTree.Filter("std::abs(pdg) == 211 && ITSTPCTrackInd >=0 && ITSTPCABTrackInd >= 0").Count();
    std::cout<<"Bad: "<<*tempcount<<std::endl;
    auto ITSTPC_pt_rec = ITSTPC_pt_rec_temp.Take<float>("pt");
    auto ITSTPCpurNotf = parTree.Filter("std::abs(pdg) == 211 && ((ITSTPCTrackInd >=0 && isFakeITSTPC <= 0 )|| (ITSTPCABTrackInd >=0 && isFakeITSTPCAB <= 0))");
    auto ITSTPCptNotf = ITSTPCpurNotf.Take<float>("pt");

    auto ITSTPC_pt_decay_temp = parTreeSigma.Filter("std::abs(pdg) == 211 && (ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0 )&& isRecoSigma == true");
    auto ITSTPC_pt_decay = ITSTPC_pt_decay_temp.Take<float>("pt_mother");

    auto ITSTPC_pt_decay_tempNotFake = parTreeSigma.Filter("std::abs(pdg) == 211 && ((ITSTPCTrackInd >=0 && isFakeITSTPC <= 0 )|| (ITSTPCABTrackInd >=0 && isFakeITSTPCAB <= 0)) && isRecoSigma == true && isFakeSigma < 0 ");
    auto ITSTPC_pt_decay_NotFake = ITSTPC_pt_decay_tempNotFake.Take<float>("pt_mother");

    auto ITSTPC_R_decay_tempNotFake = parTreeSigma.Filter("std::abs(pdg) == 211 && ((ITSTPCTrackInd >=0 && isFakeITSTPC <= 0 )|| (ITSTPCABTrackInd >=0 && isFakeITSTPCAB <= 0)) && isRecoSigma == true && isFakeSigma < 0");
    auto ITSTPC_R_decay_NotFake = ITSTPC_R_decay_tempNotFake.Take<float>("Rmother");
    auto ITSTPC_R_decay_temp = parTreeSigma.Filter("std::abs(pdg) == 211 && (ITSTPCTrackInd >=0 || ITSTPCABTrackInd >=0 ) && isRecoSigma == true");
    auto ITSTPC_R_decay = ITSTPC_R_decay_temp.Take<float>("Rmother");

    for (auto &myItem : ITSTPC_R_decay_NotFake)
    {
        ITSTPCh10.Fill(myItem);
    }

    for (auto &myItem : ITSTPC_R_decay)
    {
        ITSTPCh11.Fill(myItem);
    }

    for (auto &myItem : ITSTPC_pt_decay_NotFake)
    {
        ITSTPCh7.Fill(myItem);
    }

    for (auto &myItem : ITSTPC_pt_decay)
    {
        ITSTPCh5.Fill(myItem);
    }

    for (auto &myItem : ITSTPCptNotf)
    {
        ITSTPCh2.Fill(myItem);
    }
    for (auto &myItem : ITSTPC_pt_rec)
    {
        ITSTPCh1.Fill(myItem);
    }

    TH1F ITSTPCh3 = ITSTPCh1 / pth;
    TH1F ITSTPCh4 = ITSTPCh2 / ITSTPCh1;

    TH1F ITSTPCh6 = ITSTPCh5 / pth_mother;
    TH1F ITSTPCh8 = ITSTPCh7 / pth_mother;
    TH1F ITSTPCh9 = ITSTPCh7 / ITSTPCh5;

    TH1F ITSTPCh12 = ITSTPCh11 / r_mother;
    TH1F ITSTPCh13 = ITSTPCh10 / ITSTPCh11;
    TH1F ITSTPCh14 = ITSTPCh10 / r_mother;
    ITSTPCh14.SetTitle("Efficiency of topology reconstruction (sigma R) (Not Fake) (ITSTPC)");
    ITSTPCh13.SetTitle("Total Topology Purity (sigma R) (ITSTPC)");
    ITSTPCh12.SetTitle("Efficiency of Topology reconstruction (sigma R) (ITSTPC)");


    ITSTPCh9.SetTitle("Total Topology Purity (ITSTPC)");
    ITSTPCh6.SetTitle("Efficiency of Topology reconstruction (sigma pt) (ITSTPC)");
    ITSTPCh8.SetTitle("Efficiency of topology reconstruction (sigma pt) (Not Fake) (ITSTPC)");

    ITSTPCh3.SetTitle("Efficiency of Pion reconstruction (ITSTPC)");
    ITSTPCh4.SetTitle("Purity of Pion Track (ITSTPC)");

    ///////////////////////////////////////////////



    //////////////////////////////////////////////



    TFile outputFile("/home/justas_t/efficiency_007_pion_central.root", "recreate");

    pth.Write();
    ITSh1.Write();
    ITSh2.Write();
    ITSh3.Write();
    ITSh4.Write();
    TPCh1.Write();
    TPCh2.Write();
    TPCh3.Write();
    TPCh4.Write();
    ITSTPCh1.Write();
    ITSTPCh2.Write();
    ITSTPCh3.Write();
    ITSTPCh4.Write();
    TPCh5.Write();
    ITSTPCh5.Write();
    ITSh5.Write();
    TPCh6.Write();
    ITSTPCh6.Write();
    ITSh6.Write();
    TPCh7.Write();
    ITSTPCh7.Write();
    ITSh7.Write();
    TPCh8.Write();
    ITSTPCh8.Write();
    ITSh8.Write();
    TPCh9.Write();
    ITSTPCh9.Write();
    ITSh9.Write();
    TPCh10.Write();
    ITSTPCh10.Write();
    ITSh10.Write();
    TPCh11.Write();
    ITSTPCh11.Write();
    ITSh11.Write();
    TPCh12.Write();
    ITSTPCh12.Write();
    ITSh12.Write();
    TPCh13.Write();
    ITSTPCh13.Write();
    ITSh13.Write();
    TPCh14.Write();
    ITSTPCh14.Write();
    ITSh14.Write();

    pth_mother.Write();
    r_mother.Write();

    outputFile.Close();
}
