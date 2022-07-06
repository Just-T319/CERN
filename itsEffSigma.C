#include <iostream>
#include "ROOT/RVec.hxx"

struct ParticleInfo
{
  int event;
  int pdg;
  float pt;
  float eta;
  float phi;
  float rapid;
  int mother;
  int first;
  unsigned short clusters = 0u;
  unsigned char isReco = 0u;
  unsigned char isFake = 0u;
  bool isPrimary = 0u;
  unsigned char storedStatus = 2; /// not stored = 2, fake = 1, good = 0
  o2::its::TrackITS track;
};

void itsEffSigma()
{

  ROOT::RDataFrame par0Tree("ParticleInfo", "~/CheckTracksCA_007.root");
  auto parTree = par0Tree.Filter("std::abs(rapid) < 1");
  TH1F h1("itsEff_gen", "Eff_gen; hits", 4, 3.5, 7.5);
  TH1F h2("itsEff_rec", "Eff_rec; hits", 4, 3.5, 7.5);
  TH1F h3("itsEff", "Efficiency; hits", 4, 3.5, 7.5);
  TH1F h4("itsEff_mom_rec", "Efficiency; GeV/c", 60, 0, 15);
  TH1F h5("itsEff_mom_gen", "Efficiency; GeV/c", 60, 0, 15);
  TH1F h8("its_mom_notFake", "MomentumNotFake; GeV/c", 60, 0, 15);
  // TH1F h5("itsEff_gen_all", "Eff_gen; hits", 4, 3.5, 7.5);
  // TH1F h3("Eff_rec", "Eff_rec; hits", 4, 4, 8);
  // auto ncluster = parTree.Take<unsigned short>("clusters");
  // auto nsigma0 = parTree.Take<int>("pdg");
  // int i = 0;
  // auto smt = parTree.Filter("std::abs(pdg) == 3112 && isReco > 0");
  // auto smt1 = smt.Take<unsigned short>("clusters");
  // for (auto& myItem : smt1) {

  //   auto temp = myItem;
  //   std::cout << temp << std::endl;
  //   ++i;

  // }
  //   std::cout<<"Total: " << i << std::endl;
  //   //auto temp = nsigma0->GetValue();
  //  // std::cout << *temp << std::endl;

  auto nsigma = parTree.Filter("std::abs(pdg) == 3112 && clusters >= (pow(2,4)-1)").Count();
  float temp0 = *nsigma;

  // auto gr = new TGraph();
  // gr->SetTitle("Efficiency;Clusters;Efficiency");

  std::cout << "Total Sigma in ITS (more than 4 hits): " << *nsigma << std::endl;

  for (int i = 4; i < 8; ++i)
  {

    auto nnsigma = parTree.Filter(Form("std::abs(pdg) == 3112 && clusters == (pow(2,%i)-1)", i)).Count();
    auto crsigma = parTree.Filter(Form("std::abs(pdg) == 3112 && clusters == (pow(2,%i)-1) && isReco > 0 ", i)).Count();
    auto recofakesigma = parTree.Filter(Form("std::abs(pdg) == 3112 && clusters == (pow(2,%i)-1) && isFake > 0 ", i)).Count();
    float temp1 = *crsigma;
    float temp2 = *recofakesigma;
    float temp3 = *nnsigma;

    std::cout << "Total Sigma that has " << i << " hits: " << *nnsigma << std::endl;
    std::cout << "Total Reco tracks that has (t+f) " << i << " hits: " << (*crsigma + *recofakesigma) << std::endl;
    if (temp3 != 0)
    {

      std::cout << "Efficiency: " << (temp1 + temp2) / temp3 * 100 << "%" << std::endl;
    }
    for (int q = 0; q < temp3; q++)
    {
      h1.Fill(i);
    }

    for (int q = 0; q < (temp1 + temp2); q++)
    {
      h2.Fill(i);
    }

    // //gr->AddPoint(i,(temp0+temp1)/temp3*100);
  }

  // gr->Draw("AL*");
  auto totalgensigma = parTree.Filter("std::abs(pdg) == 3112").Count();
  float tempgen = *totalgensigma;
  std::cout << "  " << tempgen << std::endl;
  h3 = h2 / h1;
  TH1F h6 = h2 * (1 / tempgen);

  auto mom_rec = parTree.Filter("std::abs(pdg) == 3112 && (isReco > 0 || isFake > 0)");
  auto mom_gen = parTree.Filter("std::abs(pdg) == 3112");
  auto h4c = mom_rec.Take<float>("pt"); // Histo1D({"itsEff_mom_rec", "histTitle", 20, 0., 5.},
  auto h5c = mom_gen.Take<float>("pt"); // Histo1D({"itsEff_mom_rec", "histTitle", 20, 0., 5.},

  // auto pt = parTree.Take<float>("pt"); //Histo1D({"itsEff_mom_rec", "histTitle", 20, 0., 5.},
  // auto h5c = mom_gen.Take<float>("pt"); //Histo1D({"itsEff_mom_rec", "histTitle", 20, 0., 5.},
  auto purNotf = parTree.Filter("std::abs(pdg) == 3112 && isReco > 0 && isFake <= 0");
  auto purAll = parTree.Filter("std::abs(pdg) == 3112");
  auto ptNotf = purNotf.Take<float>("pt");

  for (auto &myItem : ptNotf)
  {
    h8.Fill(myItem);
  }

  for (auto &myItem : h4c)
  {
    h4.Fill(myItem);
  }

  for (auto &myItem : h5c)
  {
    h5.Fill(myItem);
  }

  TH1F h7 = h4 / h5;
  TH1F h9 = h8 / h4;

  h9.SetTitle("Purity of Sigma tracks");
  h7.SetTitle("Efficiency of Sigma tracks");

  TFile outputFile("/home/justas_t/efficiency_007_sigma_central.root", "recreate");

  // gr->Write();
  h1.Write();
  h2.Write();
  h3.Write();
  h4.Write();
  h5.Write();
  h6.Write();
  h7.Write();
  h8.Write();
  h9.Write();

  outputFile.Close();
}