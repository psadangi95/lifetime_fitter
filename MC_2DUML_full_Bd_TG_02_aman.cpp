
#include <TLatex.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <TString.h>
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooHist.h"
#include "RooGenericPdf.h"
#include "RooTruthModel.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooEffProd.h"
using namespace std;
using namespace RooFit;
void MC_2DUML_full_Bd_TG_02_aman()
{
    TFile *fin = new TFile("BdMC_2018_bdt_merge.root");
    TFile *fin1 = new TFile("BdMC_2017_bdt_merge.root");
    TFile *fin2 = new TFile("BdMC_2016_bdt_merge.root");
    TTree *tree = (TTree *)fin->Get("tree");
    TTree *tree1 = (TTree *)fin1->Get("tree");
    TTree *tree2 = (TTree *)fin2->Get("tree");
    int entries = tree->GetEntries();
    int entries1 = tree1->GetEntries();
    int entries2 = tree2->GetEntries();
    int nentries_ = entries+entries1+entries2;
    cout << "\n=> total entries in data tree = " << nentries_ << endl;

    Float_t m_t,  mhsb_t, fpuw8, Bprob_t,treco_t, trecoe_t, bdt_t = 0., fpu_w8_t, BKs0mass, BJmass, lambdamass, genb_ct;
    Float_t m_t1, mhsb_t1,fpuw81, Bprob_t1,treco_t1, trecoe_t1, bdt_t1 = 0.;
    Float_t m_t2, mhsb_t2,fpuw82, Bprob_t2,treco_t2, trecoe_t2, bdt_t2 = 0.;
    bool gmuid_t, istruebs;
    Char_t truebs, truebs1 , truebs2 ;

    tree->SetBranchAddress("Bmass", &m_t);
    tree->SetBranchAddress("BProb", &Bprob_t);
    tree->SetBranchAddress("Btrecoxy_PDG", &treco_t);
    tree->SetBranchAddress("Btrecoxye_J", &trecoe_t);
    tree->SetBranchAddress("Bdt", &bdt_t);
    tree->SetBranchAddress("is_truebs", &truebs);
    tree->SetBranchAddress("fpu_w8", &fpuw8);

    tree1->SetBranchAddress("Bmass", &m_t1);
    tree1->SetBranchAddress("BProb", &Bprob_t1);
    tree1->SetBranchAddress("Btrecoxy_PDG", &treco_t1);
    tree1->SetBranchAddress("Btrecoxye_J", &trecoe_t1);
    tree1->SetBranchAddress("Bdt", &bdt_t1);
    tree1->SetBranchAddress("is_truebs", &truebs1);
    tree1->SetBranchAddress("fpu_w8", &fpuw81);

    tree2->SetBranchAddress("Bmass", &m_t2);
    tree2->SetBranchAddress("BProb", &Bprob_t2);
    tree2->SetBranchAddress("Btrecoxy_PDG", &treco_t2);
    tree2->SetBranchAddress("Btrecoxye_J", &trecoe_t2);
    tree2->SetBranchAddress("Bdt", &bdt_t2);
    tree2->SetBranchAddress("is_truebs", &truebs2);
    tree2->SetBranchAddress("fpu_w8", &fpuw82);

    RooRealVar mass("mass", "Invariant mass(B_{s}^{0})[GeV/c^{2}]", 5.0, 5.6);
    RooRealVar masshsb("masshsb", "Invariant mass(B_{s}^{0})[GeV/c^{2}]", 5.45, 5.6);
    RooRealVar treco("treco", "DecayTime[ps]", 0.2, 10.);
    RooRealVar trecoe("trecoe", "DecayTime[ps]", 0.01, 0.3);
    RooRealVar lt("lt", "lt", 1.525, 0.35, 1.8);
    RooRealVar fpu_w8("fpu_w8", "s_wt", -10, 10);
    RooArgList arglist(mass, treco, trecoe, fpu_w8);
  
    RooDataSet *reduc_tree = new RooDataSet("reduc_tree", "", arglist,"fpu_w8");
    RooDataSet *reduc_tree1 = new RooDataSet("reduc_tree1", "", arglist,"fpu_w8");
    RooDataSet *reduc_tree2 = new RooDataSet("reduc_tree2", "", arglist,"fpu_w8");

    for (int j = 0; j < entries; j++)
    {
        tree->GetEntry(j);
        if (bdt_t <= 0.61 )
            continue;
        if(treco_t*1e12<0.2 || treco_t*1e12>10 ) continue;
        if(truebs==0) continue;
            mass.setVal(m_t);
            treco.setVal(treco_t * 1e12);
            trecoe.setVal(trecoe_t * 1e12);
            reduc_tree->add(arglist,fpuw8);
    }

    for (int j = 0; j < entries1; j++)
    {
        tree1->GetEntry(j);
        if (bdt_t1 <= 0.64  )
            continue;
        if(truebs1==0) continue;
        if(treco_t1*1e12<0.2 || treco_t1*1e12>10) continue;
            mass.setVal(m_t1);
            treco.setVal(treco_t1 * 1e12);
            trecoe.setVal(trecoe_t1 * 1e12);
            reduc_tree1->add(arglist,fpuw81);
    } 

    for (int j = 0; j < entries2; j++)
    {
        tree2->GetEntry(j);
        if (bdt_t2 <= 0.66)
        continue;
        if(truebs2==0) continue;
        if(treco_t2*1e12<0.2 || treco_t2*1e12>10) continue;
            mass.setVal(m_t2);
            treco.setVal(treco_t2 * 1e12);
            trecoe.setVal(trecoe_t2 * 1e12);
            reduc_tree2->add(arglist,fpuw82);
    }

    reduc_tree->Print();

    TCut c1 = "mass>5.0 && mass<5.6";
    RooDataSet *reduce_tree = (RooDataSet*)reduc_tree->reduce(c1);
    RooDataSet *reduce_tree1 =(RooDataSet*)reduc_tree1->reduce(c1);
    RooDataSet *reduce_tree2 =(RooDataSet*)reduc_tree2->reduce(c1);
    reduce_tree->Print();
    
    // using triple gaussians
    RooRealVar fg3_16("fg1","",7.31341e-01,0.1,0.99);
    RooRealVar fg4_16("fg2","",9.91368e-01,0.1,0.99);
    RooRealVar mean2_16("mean2_16", "mean",5.27973e+00,5.27,5.29);
    RooRealVar sigma4_16("sigma4_16", "sigma1 of Gaussian",1.13345e-02,0.0001,0.2);//
    RooRealVar sigma5_16("sigma5_16", "sigma2 of Gaussian",2.61697e-02,0.0001,0.2);//
    RooRealVar sigma6_16("sigma6_16", "sigma3 of Gaussian",1.26803e-01,0.0001,0.2);//

    RooRealVar fg3_17("fg1","",2.38446e-01,0.1,0.99);
    RooRealVar fg4_17("fg2","",9.93663e-01,0.1,0.99);
    RooRealVar mean2_17("mean2_17", "mean",5.28001e+00,5.27,5.29);
    RooRealVar sigma4_17("sigma4_17", "sigma1 of Gaussian",2.56142e-02,0.0001,0.2);
    RooRealVar sigma5_17("sigma5_17", "sigma2 of Gaussian",1.05222e-02,0.0001,0.2);
    RooRealVar sigma6_17("sigma6_17", "sigma3 of Gaussian",2.10193e-01,0.0001,0.2);

    RooRealVar fg3_18 ("fg3_18","",2.44072e-01,0.1,0.99);
    RooRealVar fg4_18 ("fg4_18","",9.88339e-01,0.1,0.99);
    RooRealVar mean2_18 ("mean2_18","",5.27978e+00,5.27,5.29);
    RooRealVar sigma4_18 ("sigma4_18","",2.49649e-02,0.0001,0.2);
    RooRealVar sigma5_18 ("sigma5_18","",1.06898e-02,0.0001,0.2);
    RooRealVar sigma6_18 ("sigma6_18","",1.09579e-01,0.0001,0.2);

    // gaussian fit from Bd

    RooGaussian gausian4_16("gausian4_16", "gaus4", mass, mean2_16, sigma4_16);
    RooGaussian gausian5_16("gausian5_16", "gaus5", mass, mean2_16, sigma5_16);
    RooGaussian gausian6_16("gausian6_16", "gaus6", mass, mean2_16, sigma6_16);

    RooGaussian gausian4_17("gausian4_17", "gaus4", mass, mean2_17, sigma4_17);
    RooGaussian gausian5_17("gausian5_17", "gaus5", mass, mean2_17, sigma5_17);
    RooGaussian gausian6_17("gausian6_17", "gaus6", mass, mean2_17, sigma6_17);

    RooGaussian gausian4_18("gausian4_18", "gaus4", mass, mean2_18, sigma4_18);
    RooGaussian gausian5_18("gausian5_18", "gaus5", mass, mean2_18, sigma5_18);
    RooGaussian gausian6_18("gausian6_18", "gaus6", mass, mean2_18, sigma6_18);
    
    RooAddPdf mass216("mass216", "mass2", RooArgList(gausian4_16, gausian5_16), fg3_16);
    RooAddPdf mass217("mass217", "mass2", RooArgList(gausian4_17, gausian5_17), fg3_17);
    RooAddPdf mass218("mass218", "mass2", RooArgList(gausian4_18, gausian5_18), fg3_18);

    RooAddPdf massBd16("massBd16", "mass5", RooArgList(mass216, gausian6_16), fg4_16);
    RooAddPdf massBd17("massBd17", "mass5", RooArgList(mass217, gausian6_17), fg4_17);
    RooAddPdf massBd18("massBd18", "mass5", RooArgList(mass218, gausian6_18), fg4_18);

    
    RooTruthModel gm("gm", "truth model", treco);
    RooRealVar mean("mean","",0.0);
    RooRealVar scale("scale","",1.0);

    //RooGaussModel gm("gm", "truth model", treco, mean, scale, trecoe);
    
    RooRealVar ltBd("ltBd", "lt", 1.519, 1.35, 4.4);

    RooDecay *decay_Bd16 = new RooDecay("decay_Bd16", "decay", treco, ltBd, gm, RooDecay::SingleSided);
    RooDecay *decay_Bd17 = new RooDecay("decay_Bd17", "decay", treco, ltBd, gm, RooDecay::SingleSided);
    RooDecay *decay_Bd18 = new RooDecay("decay_Bd18", "decay", treco, ltBd, gm, RooDecay::SingleSided);

    RooRealVar dl1bd66("dl1bd66", "effpdf meanbd", -0.01021); // define variables for the fit
    RooRealVar dl2bd66("dl2bd66", "effpdf sigbd", -4.022e-11);
    RooRealVar powlbd66("powlbd66", "powlbd", 0.03343);
    RooRealVar powl2bd66("powl2bd66", "powl2bd", 3.443);
    RooRealVar powl3bd66("powl3bd66", "powl3bd", -0.0001924);
    RooFormulaVar effbd16("effbd16", "(dl1bd66+dl2bd66*treco  + powlbd66 /(1+exp( -treco*powl2bd66)) +powl3bd66*treco*treco)", RooArgSet(treco, dl1bd66, dl2bd66, powlbd66, powl2bd66, powl3bd66));
    RooEffProd decaytime_bd16 ("decaytime_bd16", "", *decay_Bd16, effbd16);

    RooRealVar dl1bd77("dl1bd77","effpdf meanbd77",-0.01034); //define variables for the fit
    RooRealVar dl2bd77("dl2bd77","effpdf sigbd",-0.0002096);
    RooRealVar powlbd77("powlbd77","powlbd",0.02853);
    RooRealVar powl2bd77("powl2bd77","powl2bd",4.228);
    RooRealVar powl3bd77("powl3bd77","powl3bd",-0.00007697);
    RooFormulaVar effbd77("effbd77", "(dl1bd77+dl2bd77*treco  + powlbd77 /(1+exp( -treco*powl2bd77)) +powl3bd77*treco*treco)", RooArgSet(treco, dl1bd77,dl2bd77, powlbd77, powl2bd77, powl3bd77));
    RooEffProd decaytime_bd17("decaytime_bd17", "model with efficiencybd77", *decay_Bd17, effbd77);
    
    RooRealVar dl1bd88("dl1bd88","effpdf meanbd",-0.01221); //define variables for the fit
    RooRealVar dl2bd88("dl2bd88","effpdf sigbd",-4.263e-13);
    RooRealVar powlbd88("powlbd88","powlbd",0.03354);
    RooRealVar powl2bd88("powl2bd88","powl2bd",4.508);
    RooRealVar powl3bd88("powl3bd88","powl3bd",-0.000144);
    RooFormulaVar effbd88("effbd88", "(dl1bd88+dl2bd88*treco  + powlbd88 /(1+exp( -treco*powl2bd88)) +powl3bd88*treco*treco)", RooArgSet(treco, dl1bd88,dl2bd88, powlbd88, powl2bd88, powl3bd88));
    RooEffProd decaytime_bd18("decaytime_bd18", "model with efficiency", *decay_Bd18, effbd88);

    RooProdPdf *Bd16 = new RooProdPdf("Bd16", "", RooArgSet(massBd16, decaytime_bd16));
    RooProdPdf *Bd17 = new RooProdPdf("Bd17", "", RooArgSet(massBd17, decaytime_bd17));
    RooProdPdf *Bd18 = new RooProdPdf("Bd18", "", RooArgSet(massBd18, decaytime_bd18));
  
    RooCategory era("era","era") ;
    era.defineType("2018") ;
    era.defineType("2017") ;
    era.defineType("2016") ;
    RooDataSet combData("combData","combined data",RooArgSet(treco,trecoe, mass),Index(era),Import("2018",*reduce_tree),Import("2017",*reduce_tree1),Import("2016",*reduce_tree2)) ;
    RooSimultaneous simPdf("simPdf","simultaneous pdf",era) ;

    simPdf.addPdf(*Bd18,"2018") ;
    simPdf.addPdf(*Bd17,"2017") ;
    simPdf.addPdf(*Bd16,"2016") ;

    RooFitResult *fit = simPdf.fitTo(combData,  ConditionalObservables(trecoe),Save(true));

}
