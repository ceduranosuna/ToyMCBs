#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <RooGaussian.h>
#include <TCanvas.h>
#include <RooPlot.h>
#include <TAxis.h>
#include <RooChebychev.h>
#include <RooMCStudy.h>
#include <RooAddPdf.h>
#include <RooGenericPdf.h>
#include <RooEffProd.h>
#include <RooFormulaVar.h>
#include <RooEffProd.h>
#include <RooProdPdf.h>
#include <RooPolynomial.h>
#include <RooArgList.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <RooFitResult.h>
#include <RooClassFactory.h>
#include <RooAddition.h>
#include <RooMinuit.h>
#include <RooRandom.h>
#include <string.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <TChain.h>
#include "TMatrixDSym.h"

#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"

#include "Datos_LamBout.C"

using namespace RooFit ;
using namespace TMath;
using namespace std;

int cpu = 1;

//# ifndef __CINT__
//int main()
void toyMC_signal()
{

    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().getStream(1).removeTopic(Integration) ;
    RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
    RooMsgService::instance().getStream(1).removeTopic(Fitting) ;
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
    RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
    RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
    RooMsgService::instance().getStream(1).removeTopic(Eval) ;   
    RooMsgService::instance().Print() ;
    
    Double_t Msigmamin = 5.6194-3.5*0.0165;//5.3217
    Double_t Msigmamax=  5.6194+3.5*0.0165;//5.4165
    Double_t infMass = 5.6194-10*0.0165; //Double_t infMass = mass.getMin();
    Double_t supMass = 5.6194+10*0.0165; //Double_t supMass = mass.getMax();
    RooRealVar cos_l("cos_l","cos#theta_{#Lambda}",-1.,1.);
    RooRealVar cos_p("cos_p","cos#theta_{p}",-1.,1.);
    RooRealVar cos_m("cos_m","cos#theta_{#mu}",-1.,1.);
    RooRealVar mass("mass","mass",5.4,5.84,"GeV");  
    
    //mass models
    RooRealVar mean("#mu","mean of  gaussian",5.40,5.84);
    RooRealVar s1("s1","sigma 1",0.,1);
    RooRealVar s2("s2","sigma 2",0.,1.);
    RooRealVar f("f","fraction of gaussian",0.,1.);
   
    RooRealVar s1_a("s1_a","asigma 1",0.,1);
    RooRealVar s2_a("s2_a","asigma 2",0.,1.);
    RooRealVar f_a("f_a","fraction of gaussian",0.,1.);
    
    RooRealVar s1_X("s1_X","sigma 1 X",0.,1);
    RooRealVar s2_X("s2_X","sigma 2 X",0.,1.);
    RooRealVar f_X("f_X","fraction of gaussian X",0.,1.);
    
    RooRealVar s1_a_X("s1_a_X","a sigma 1 X",0.,1);
    RooRealVar s2_a_X("s2_a_X","a sigma 2 X",0.,1.);
    RooRealVar f_a_X("f_a_X","fraction of gaussian X",0.,1.);
    
    RooRealVar p_("p_","pendiente poloinomio bkg",-1.,1.);
    RooRealVar p_a("p_a","pendiente poloinomio bkg",-1.,1.);
    RooRealVar p_X("p_X","pendiente poloinomio bkg X",-1.,1.);
    RooRealVar p_a_X("p_a_X","pendiente poloinomio bkg X",-1.,1.);
    
    RooGaussian G_one("G_one","Gaussian 1",mass,mean,s1);
    RooGaussian G_two("G_two","Gaussian 2",mass,mean,s2);
    
    RooGaussian G_one_a("G_one_a","Gaussian 1",mass,mean,s1_a);
    RooGaussian G_two_a("G_two_a","Gaussian 2",mass,mean,s2_a);
    
    RooGaussian G_one_X("G_one_X","Gaussian 1",mass,mean,s1_X);
    RooGaussian G_two_X("G_two_X","Gaussian 2",mass,mean,s2_X);
    
    RooGaussian G_one_a_X("G_one_a_X","a Gaussian 1 X",mass,mean,s1_a_X);
    RooGaussian G_two_a_X("G_two_a_X","a Gaussian 2 X",mass,mean,s2_a_X);
    
    RooAddPdf DGauss("DGauss","g_one+g_two",RooArgList(G_one,G_two),f);
    RooAddPdf DGauss_a("DGauss_a","ag_one+ag_two",RooArgList(G_one_a,G_two_a),f_a);
    RooAddPdf DGauss_X("DGauss_X","g_oneX+g_twoX",RooArgList(G_one_X,G_two_X),f_X);
    RooAddPdf DGauss_a_X("DGauss_a_X","ag_oneX+ag_twXo",RooArgList(G_one_a_X,G_two_a_X),f_a_X);
    
    RooPolynomial poly("poly","Background_mass",mass,p_);
    RooPolynomial poly_a("poly_a","Background_mass",mass,p_a);
    RooPolynomial poly_X("poly_X","Background_mass",mass,p_X);
    RooPolynomial poly_a_X("poly_a_X","Background_mass",mass,p_a_X);
    
    RooRealVar N_sig("N_sig","",800,0.,10000.);
    RooRealVar N_bkg("N_bkg","",1000,0.,10000.);
    RooRealVar N_sig_a("N_sig_a","",800,0.,10000.);
    RooRealVar N_bkg_a("N_bkg_a","",1000,0.,10000.);
    RooRealVar N_sig_X("N_sig_X","",800,0.,10000.);
    RooRealVar N_bkg_X("N_bkg_X","",1000,0.,10000.);
    RooRealVar N_sig_a_X("N_sig_a_X","",800,0.,10000.);
    RooRealVar N_bkg_a_X("N_bkg_a_X","",1000,0.,10000.);
    
    RooArgList NList;
    NList.add(N_sig);
    NList.add(N_bkg);
    
    RooArgList NList_a;
    NList_a.add(N_sig_a);
    NList_a.add(N_bkg_a);
    
    RooArgList NList_X;
    NList_X.add(N_sig_X);
    NList_X.add(N_bkg_X);
    
    RooArgList NList_a_X;
    NList_a_X.add(N_sig_a_X);
    NList_a_X.add(N_bkg_a_X);
    
    RooAddPdf PDFT_mass("MassPro","mass Projection",RooArgList(DGauss,poly),NList);
    RooAddPdf PDFT_mass_a("MassPro_a","mass Projection",RooArgList(DGauss_a,poly_a),NList_a);
    RooAddPdf PDFT_mass_X("MassPro_X","mass Projection",RooArgList(DGauss_X,poly_X),NList_X);
    RooAddPdf PDFT_mass_a_X("MassPro_a_X","mass Projection",RooArgList(DGauss_a_X,poly_a_X),NList_a_X);
    
    //Angular models
    
    /* ------------- Eficciency ------------------- */
    //Lambda
    RooRealVar A1("A1","A1",0);
    RooRealVar A2("A2","A2",0);
    RooRealVar A3("A3","A3",0);
    RooRealVar A4("A4","A4",0);
    RooRealVar A5("A5","A5",0);
    //RooRealVar A6("A6","A6",0);
    
    RooRealVar B1("B1","B1",0);
    RooRealVar B2("B2","B2",0);
    RooRealVar B3("B3","B3",0);
    RooRealVar B4("B4","B4",0);
    RooRealVar B5("B5","B5",0);
    
    RooRealVar C1("C1","C1",0);
    RooRealVar C2("C2","C2",0);
    RooRealVar C3("C3","C3",0);
    RooRealVar C4("C4","C4",0);
    RooRealVar C5("C5","C5",0);
    RooRealVar C6("C6","C6",0);
    
    // LambdaBar
    RooRealVar aA1("aA1","aA1",0);
    RooRealVar aA2("aA2","aA2",0);
    RooRealVar aA3("aA3","aA3",0);
    RooRealVar aA4("aA4","aA4",0);
    RooRealVar aA5("aA5","aA5",0);
    //RooRealVar aA6("aA6","aA6",0);
    
    RooRealVar aB1("aB1","aB1",0);
    RooRealVar aB2("aB2","aB2",0);
    RooRealVar aB3("aB3","aB3",0);
    RooRealVar aB4("aB4","aB4",0);
    RooRealVar aB5("aB5","aB5",0);
    
    RooRealVar aC1("aC1","aC1",0);
    RooRealVar aC2("aC2","aC2",0);
    RooRealVar aC3("aC3","aC3",0);
    RooRealVar aC4("aC4","aC4",0);
    RooRealVar aC5("aC5","aC5",0);
    RooRealVar aC6("aC6","aC6",0);
    
    //LambdaX
    RooRealVar A1X("A1X","A1X",0);
    RooRealVar A2X("A2X","A2X",0);
    RooRealVar A3X("A3X","A3X",0);
    RooRealVar A4X("A4X","A4X",0);
    RooRealVar A5X("A5X","A5X",0);
    //RooRealVar A6("A6","A6",0);
    
    RooRealVar B1X("B1X","B1X",0);
    RooRealVar B2X("B2X","B2X",0);
    RooRealVar B3X("B3X","B3X",0);
    RooRealVar B4X("B4X","B4X",0);
    RooRealVar B5X("B5X","B5X",0);
    
    RooRealVar C1X("C1X","C1X",0);
    RooRealVar C2X("C2X","C2X",0);
    RooRealVar C3X("C3X","C3X",0);
    RooRealVar C4X("C4X","C4X",0);
    RooRealVar C5X("C5X","C5X",0);
    RooRealVar C6X("C6X","C6X",0);
    
    // LambdaBarX
    RooRealVar aA1X("aA1X","aA1X",0);
    RooRealVar aA2X("aA2X","aA2X",0);
    RooRealVar aA3X("aA3X","aA3X",0);
    RooRealVar aA4X("aA4X","aA4X",0);
    RooRealVar aA5X("aA5X","aA5X",0);
    //RooRealVar aA6X("aA6X","aA6X",0);
    
    RooRealVar aB1X("aB1X","aB1X",0);
    RooRealVar aB2X("aB2X","aB2X",0);
    RooRealVar aB3X("aB3X","aB3X",0);
    RooRealVar aB4X("aB4X","aB4X",0);
    RooRealVar aB5X("aB5X","aB5X",0);
    
    RooRealVar aC1X("aC1X","aC1X",0);
    RooRealVar aC2X("aC2X","aC2X",0);
    RooRealVar aC3X("aC3X","aC3X",0);
    RooRealVar aC4X("aC4X","aC4X",0);
    RooRealVar aC5X("aC5X","aC5X",0);
    RooRealVar aC6X("aC6X","aC6X",0);
    
    RooArgList eps_lCoef;
    eps_lCoef.add(A1);
    eps_lCoef.add(cos_l);
    eps_lCoef.add(A2);
    eps_lCoef.add(A3);
    eps_lCoef.add(A4);
    eps_lCoef.add(A5);
    //eps_lCoef.add(A6);
    
    RooArgList aLeps_lCoef;
    aLeps_lCoef.add(aA1);
    aLeps_lCoef.add(cos_l);
    aLeps_lCoef.add(aA2);
    aLeps_lCoef.add(aA3);
    aLeps_lCoef.add(aA4);
    aLeps_lCoef.add(aA5);
    //aLeps_lCoef.add(aA6);
    
    RooArgList eps_lCoefX;
    eps_lCoefX.add(A1X);
    eps_lCoefX.add(cos_l);
    eps_lCoefX.add(A2X);
    eps_lCoefX.add(A3X);
    eps_lCoefX.add(A4X);
    eps_lCoefX.add(A5X);
    //eps_lCoef.add(A6);
    
    RooArgList aLeps_lCoefX;
    aLeps_lCoefX.add(aA1X);
    aLeps_lCoefX.add(cos_l);
    aLeps_lCoefX.add(aA2X);
    aLeps_lCoefX.add(aA3X);
    aLeps_lCoefX.add(aA4X);
    aLeps_lCoefX.add(aA5X);
    //aLeps_lCoef.add(aA6);
    
    RooFormulaVar eps_l("eps_l","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) )",eps_lCoef);
    RooFormulaVar aLeps_l("aLeps_l","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) )",aLeps_lCoef);
    RooFormulaVar eps_lX("eps_lX","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) )",eps_lCoefX);
    RooFormulaVar aLeps_lX("aLeps_lX","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) )",aLeps_lCoefX);
    
    RooArgList eps_pCoef;
    eps_pCoef.add(B1);
    eps_pCoef.add(cos_p);
    eps_pCoef.add(B2);
    eps_pCoef.add(B3);
    eps_pCoef.add(B4);
    eps_pCoef.add(B5);
    
    RooArgList aLeps_pCoef;
    aLeps_pCoef.add(aB1);
    aLeps_pCoef.add(cos_p);
    aLeps_pCoef.add(aB2);
    aLeps_pCoef.add(aB3);
    aLeps_pCoef.add(aB4);
    aLeps_pCoef.add(aB5);
    
    RooArgList eps_pCoefX;
    eps_pCoefX.add(B1X);
    eps_pCoefX.add(cos_p);
    eps_pCoefX.add(B2X);
    eps_pCoefX.add(B3X);
    eps_pCoefX.add(B4X);
    eps_pCoefX.add(B5X);
    
    RooArgList aLeps_pCoefX;
    aLeps_pCoefX.add(aB1X);
    aLeps_pCoefX.add(cos_p);
    aLeps_pCoefX.add(aB2X);
    aLeps_pCoefX.add(aB3X);
    aLeps_pCoefX.add(aB4X);
    aLeps_pCoefX.add(aB5X);
    
    RooFormulaVar eps_p("eps_p","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) )",eps_pCoef);
    RooFormulaVar aLeps_p("aLeps_p","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) )",aLeps_pCoef);
    RooFormulaVar eps_pX("eps_pX","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) )",eps_pCoefX);
    RooFormulaVar aLeps_pX("aLeps_pX","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) )",aLeps_pCoefX);
    
    RooArgList eps_mCoef;
    eps_mCoef.add(C1);
    eps_mCoef.add(cos_m);
    eps_mCoef.add(C2);
    eps_mCoef.add(C3);
    eps_mCoef.add(C4);
    eps_mCoef.add(C5);
    eps_mCoef.add(C6);
    
    RooArgList aLeps_mCoef;
    aLeps_mCoef.add(aC1);
    aLeps_mCoef.add(cos_m);
    aLeps_mCoef.add(aC2);
    aLeps_mCoef.add(aC3);
    aLeps_mCoef.add(aC4);
    aLeps_mCoef.add(aC5);
    aLeps_mCoef.add(aC6);
    
    RooArgList eps_mCoefX;
    eps_mCoefX.add(C1X);
    eps_mCoefX.add(cos_m);
    eps_mCoefX.add(C2X);
    eps_mCoefX.add(C3X);
    eps_mCoefX.add(C4X);
    eps_mCoefX.add(C5X);
    eps_mCoefX.add(C6X);
    
    RooArgList aLeps_mCoefX;
    aLeps_mCoefX.add(aC1X);
    aLeps_mCoefX.add(cos_m);
    aLeps_mCoefX.add(aC2X);
    aLeps_mCoefX.add(aC3X);
    aLeps_mCoefX.add(aC4X);
    aLeps_mCoefX.add(aC5X);
    aLeps_mCoefX.add(aC6X);
    
    
    RooFormulaVar eps_m("eps_m","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) + @6*(32.0*@1*@1*@1*@1*@1*@1-48.0*@1*@1*@1*@1+18.0*@1*@1-1.0) )",eps_mCoef);
    RooFormulaVar aLeps_m("aLeps_m","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) + @6*(32.0*@1*@1*@1*@1*@1*@1-48.0*@1*@1*@1*@1+18.0*@1*@1-1.0) )",aLeps_mCoef);
    RooFormulaVar eps_mX("eps_mX","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) + @6*(32.0*@1*@1*@1*@1*@1*@1-48.0*@1*@1*@1*@1+18.0*@1*@1-1.0) )",eps_mCoefX);
    RooFormulaVar aLeps_mX("aLeps_mX","( 1.0 + @0*@1 + @2*(2.0*@1*@1-1.0) + @3*(4.0*@1*@1*@1-3.0*@1) + @4*(8.0*@1*@1*@1*@1-8.0*@1*@1+1.0) + @5*(16.0*@1*@1*@1*@1*@1-20.0*@1*@1*@1+5.0*@1) + @6*(32.0*@1*@1*@1*@1*@1*@1-48.0*@1*@1*@1*@1+18.0*@1*@1-1.0) )",aLeps_mCoefX);
    
    
    RooArgList eps_Coef;
    eps_Coef.add(eps_l);
    eps_Coef.add(eps_p);
    eps_Coef.add(eps_m);
    
    RooArgList aLeps_Coef;
    aLeps_Coef.add(aLeps_l);
    aLeps_Coef.add(aLeps_p);
    aLeps_Coef.add(aLeps_m);
    
    RooArgList eps_CoefX;
    eps_CoefX.add(eps_lX);
    eps_CoefX.add(eps_pX);
    eps_CoefX.add(eps_mX);
    
    RooArgList aLeps_CoefX;
    aLeps_CoefX.add(aLeps_lX);
    aLeps_CoefX.add(aLeps_pX);
    aLeps_CoefX.add(aLeps_mX);
    
    RooFormulaVar Eps("Eps","@0*@1*@2",eps_Coef);///3.08414
    RooFormulaVar aLEps("aLEps","@0*@1*@2",aLeps_Coef);///3.05274
    RooFormulaVar EpsX("EpsX","@0*@1*@2",eps_CoefX);///3.08414
    RooFormulaVar aLEpsX("aLEpsX","@0*@1*@2",aLeps_CoefX);///3.05274
    
    RooAbsReal* I_Eps = Eps.createIntegral(RooArgSet(cos_l,cos_p,cos_m));
    RooAbsReal* I_aLEps = aLEps.createIntegral(RooArgSet(cos_l,cos_p,cos_m));
    RooAbsReal* I_EpsX = EpsX.createIntegral(RooArgSet(cos_l,cos_p,cos_m));
    RooAbsReal* I_aLEpsX = aLEpsX.createIntegral(RooArgSet(cos_l,cos_p,cos_m));
    
    RooArgSet *epsParam = Eps.getParameters(RooArgSet(cos_l,cos_p,cos_m));
    RooArgSet *aLepsParam = aLEps.getParameters(RooArgSet(cos_l,cos_p,cos_m));
    RooArgSet *epsParamX = EpsX.getParameters(RooArgSet(cos_l,cos_p,cos_m));
    RooArgSet *aLepsParamX = aLEpsX.getParameters(RooArgSet(cos_l,cos_p,cos_m));
    
    RooFormulaVar eps("eps","@0/@1",RooArgSet(Eps,*I_Eps));///3.08414
    RooFormulaVar aLeps("aLeps","@0/@1",RooArgSet(aLEps,*I_aLEps));///3.05274
    RooFormulaVar epsX("epsX","@0/@1",RooArgSet(EpsX,*I_EpsX));///3.08414
    RooFormulaVar aLepsX("aLepsX","@0/@1",RooArgSet(aLEpsX,*I_aLEpsX));///3.05274
    
    epsParam->readFromFile("textFiles/eff_lambda_7TeV.txt","READ");
    aLepsParam->readFromFile("textFiles/eff_antilambda_7TeV.txt","READ");
    epsParamX->readFromFile("textFiles/eff_lambda_8TeV.txt","READ") ;
    aLepsParamX->readFromFile("textFiles/eff_antilambda_8TeV.txt","READ") ;
    /*epsParam->Print("v");
    aLepsParam->Print("v");
    epsParamX->Print("v");
    aLepsParamX->Print("v");*/
    
    
    A1.setConstant(kTRUE);    aA1.setConstant(kTRUE);
    A2.setConstant(kTRUE);    aA2.setConstant(kTRUE);
    A3.setConstant(kTRUE);    aA3.setConstant(kTRUE);
    A4.setConstant(kTRUE);    aA4.setConstant(kTRUE);
    A5.setConstant(kTRUE);    aA5.setConstant(kTRUE);
    //A6.setConstant(kTRUE);    aA6.setConstant(kTRUE);
    
    B1.setConstant(kTRUE);    aB1.setConstant(kTRUE);
    B2.setConstant(kTRUE);    aB2.setConstant(kTRUE);
    B3.setConstant(kTRUE);    aB3.setConstant(kTRUE);
    B4.setConstant(kTRUE);    aB4.setConstant(kTRUE);
    B5.setConstant(kTRUE);    aB5.setConstant(kTRUE);
    
    C1.setConstant(kTRUE);    aC1.setConstant(kTRUE);
    C2.setConstant(kTRUE);    aC2.setConstant(kTRUE);
    C3.setConstant(kTRUE);    aC3.setConstant(kTRUE);
    C4.setConstant(kTRUE);    aC4.setConstant(kTRUE);
    C5.setConstant(kTRUE);    aC5.setConstant(kTRUE);
    C6.setConstant(kTRUE);    aC6.setConstant(kTRUE);
    //std::cout << aA_one.getVal() << std::endl; return;
    
    A1X.setConstant(kTRUE);    aA1X.setConstant(kTRUE);
    A2X.setConstant(kTRUE);    aA2X.setConstant(kTRUE);
    A3X.setConstant(kTRUE);    aA3X.setConstant(kTRUE);
    A4X.setConstant(kTRUE);    aA4X.setConstant(kTRUE);
    A5X.setConstant(kTRUE);    aA5X.setConstant(kTRUE);
    //A6.setConstant(kTRUE);    aA6.setConstant(kTRUE);
    
    B1X.setConstant(kTRUE);    aB1X.setConstant(kTRUE);
    B2X.setConstant(kTRUE);    aB2X.setConstant(kTRUE);
    B3X.setConstant(kTRUE);    aB3X.setConstant(kTRUE);
    B4X.setConstant(kTRUE);    aB4X.setConstant(kTRUE);
    B5X.setConstant(kTRUE);    aB5X.setConstant(kTRUE);
    
    C1X.setConstant(kTRUE);    aC1X.setConstant(kTRUE);
    C2X.setConstant(kTRUE);    aC2X.setConstant(kTRUE);
    C3X.setConstant(kTRUE);    aC3X.setConstant(kTRUE);
    C4X.setConstant(kTRUE);    aC4X.setConstant(kTRUE);
    C5X.setConstant(kTRUE);    aC5X.setConstant(kTRUE);
    C6X.setConstant(kTRUE);    aC6X.setConstant(kTRUE);
    
    /* Angular Distribution */
    
    RooRealVar a_l("a_l","Alpha_lambda",0.642);
    //a_l.setVal(0.629);
    //a_l.setVal(0.655);
    //if (sample <0) a_l.setVal(-1.*a_l.getVal());
    a_l.setConstant(kTRUE);
    RooRealVar P("P","Polarizacion"/*,0.05*/,-1.5,1.5);
    RooRealVar alpha_one("alpha1","alpha1",/*,0.07*/ -1.5,1.5);
    RooRealVar alpha_two("alpha2","alpha2",-1.2,1.2);
    RooRealVar gamma("gamma","gamma",/*,-0.59*/ -0.5, -2.1,1.1);
     
    RooArgList FSigCoefList;
    FSigCoefList.add(alpha_two);
    FSigCoefList.add(a_l);
    FSigCoefList.add(cos_p);
    FSigCoefList.add(alpha_one);
    FSigCoefList.add(P);
    FSigCoefList.add(cos_l);
    FSigCoefList.add(gamma);
    FSigCoefList.add(cos_m); 

    RooGenericPdf FSig("FSig","Angular Distribution","(1.+(@0*@1*@2)-(@3*@4*@5)-((1.+2.*@6)*@1*@4*@5*@2)/3.+(@6*(3.*@7*@7-1.))/4.+(@1*(3.*@3-@0)*@2*(3.*@7*@7-1.))/8.+((@3-3.*@0)*@4*@5*(3.*@7*@7-1.))/8.+((@6-4.)*@1*@4*@5*@2*(3.*@7*@7-1.))/12.)",FSigCoefList);  
    // Background....
    
    RooRealVar D1("D1","D1",-10,10);
    RooRealVar D2("D2","D2",-10,10);
    RooRealVar D3("D3","D3",-10,10);
    RooRealVar D4("D4","D4",-10,10);
    RooRealVar D5("D5","D5",-10,10);
    RooRealVar D6("D6","D6",-10,10);
    RooRealVar D7("D7","D7",-10,10);
    
    RooRealVar aD1("aD1","aD1",-10,10);
    RooRealVar aD2("aD2","aD2",-10,10);
    RooRealVar aD3("aD3","aD3",-10,10);
    RooRealVar aD4("aD4","aD4",-10,10);
    RooRealVar aD5("aD5","aD5",-10,10);
    RooRealVar aD6("aD6","aD6",-10,10);
    RooRealVar aD7("aD7","aD7",-10,10);
    
    RooRealVar D1X("D1X","D1X",-10,10);
    RooRealVar D2X("D2X","D2X",-10,10);
    RooRealVar D3X("D3X","D3X",-10,10);
    RooRealVar D4X("D4X","D4X",-10,10);
    RooRealVar D5X("D5X","D5X",-10,10);
    RooRealVar D6X("D6X","D6X",-10,10);
    RooRealVar D7X("D7X","D7X",-10,10);
    
    RooRealVar aD1X("aD1X","aD1X",-10,10);
    RooRealVar aD2X("aD2X","aD2X",-10,10);
    RooRealVar aD3X("aD3X","aD3X",-10,10);
    RooRealVar aD4X("aD4X","aD4X",-10,10);
    RooRealVar aD5X("aD5X","aD5X",-10,10);
    RooRealVar aD6X("aD6X","aD6X",-10,10);
    RooRealVar aD7X("aD7X","aD7X",-10,10);
    
    RooRealVar E1("E1","E1",-10.,10.);
    RooRealVar E2("E2","E2",-10.,10.);//Coeficientes de la distrib bkg para tetha_p
    RooRealVar E3("E3","E3",-10.,10.);
    
    RooRealVar aE1("aE1","aE1",-10.,10.);
    RooRealVar aE2("aE2","aE2",-10.,10.);//Coeficientes de la distrib bkg para tetha_p
    RooRealVar aE3("aE3","aE3",-10.,10.);
    RooRealVar aE4("aE4","aE4",-10.,10.);
    
    RooRealVar E1X("E1X","E1",-10.,10.);
    RooRealVar E2X("E2X","E2",-10.,10.);//Coeficientes de la distrib bkg para tetha_p
    RooRealVar E3X("E3X","E3",-10.,10.);
    
    RooRealVar aE1X("aE1X","aE1",-10.,10.);
    RooRealVar aE2X("aE2X","aE2",-10.,10.);//Coeficientes de la distrib bkg para tetha_p
    RooRealVar aE3X("aE3X","aE3",-10.,10.);
    RooRealVar aE4X("aE4X","aE4X",-10.,10.);
    
    RooRealVar k1("k1","k1",-10,10);
    RooRealVar k2("k2","k2",-10,10);
    RooRealVar k3("k3","k3",-10,10);
    RooRealVar k4("k4","k4",-10,10);
    
    RooRealVar ak1("ak1","ak1",-10,10);
    RooRealVar ak2("ak2","ak2",-10,10);
    RooRealVar ak3("ak3","ak3",-10,10);
    RooRealVar ak4("ak4","ak4",-10,10);
    
    RooRealVar k1X("k1X","k1X",-10,10);
    RooRealVar k2X("k2X","k2X",-10,10);
    RooRealVar k3X("k3X","k3X",-10,10);
    RooRealVar k4X("k4X","k4X",-10,10);
    
    RooRealVar ak1X("ak1X","ak1X",-10,10);
    RooRealVar ak2X("ak2X","ak2X",-10,10);
    RooRealVar ak3X("ak3X","ak3X",-10,10);
    RooRealVar ak4X("ak4X","ak4X",-10,10);
    
    RooArgList Fb_lCoef;
    Fb_lCoef.add(D1);
    Fb_lCoef.add(D2);
    Fb_lCoef.add(D3);
    Fb_lCoef.add(D4);
    Fb_lCoef.add(D5);
    Fb_lCoef.add(D6);
    Fb_lCoef.add(D7);
    
    RooArgList aFb_lCoef;
    aFb_lCoef.add(aD1);
    aFb_lCoef.add(aD2);
    aFb_lCoef.add(aD3);
    aFb_lCoef.add(aD4);
    aFb_lCoef.add(aD5);
    aFb_lCoef.add(aD6);
    aFb_lCoef.add(aD7);
    
    RooArgList Fb_lXCoef;
    Fb_lXCoef.add(D1X);
    Fb_lXCoef.add(D2X);
    Fb_lXCoef.add(D3X);
    Fb_lXCoef.add(D4X);
    Fb_lXCoef.add(D5X);
    Fb_lXCoef.add(D6X);
    Fb_lXCoef.add(D7X);
    
    RooArgList aFb_lXCoef;
    aFb_lXCoef.add(aD1X);
    aFb_lXCoef.add(aD2X);
    aFb_lXCoef.add(aD3X);
    aFb_lXCoef.add(aD4X);
    aFb_lXCoef.add(aD5X);
    aFb_lXCoef.add(aD6X);
    aFb_lXCoef.add(aD7X);
    
    RooArgList Fb_pCoef;
    Fb_pCoef.add(E1);
    Fb_pCoef.add(E2);
    Fb_pCoef.add(E3);
    
    RooArgList aLFb_pCoef;
    aLFb_pCoef.add(aE1);
    aLFb_pCoef.add(aE2);
    aLFb_pCoef.add(aE3);
    aLFb_pCoef.add(aE4);
    
    RooArgList Fb_pCoefX;
    Fb_pCoefX.add(E1X);
    Fb_pCoefX.add(E2X);
    Fb_pCoefX.add(E3X);
    
    RooArgList aLFb_pCoefX;
    aLFb_pCoefX.add(aE1X);
    aLFb_pCoefX.add(aE2X);
    aLFb_pCoefX.add(aE3X);
    aLFb_pCoefX.add(aE4X);
    
    RooArgList Fb_mCoef;
    Fb_mCoef.add(cos_m);
    Fb_mCoef.add(k1);
    Fb_mCoef.add(k2);
    Fb_mCoef.add(k3);
    Fb_mCoef.add(k4);
    
    RooArgList aFb_mCoef;
    aFb_mCoef.add(cos_m);
    aFb_mCoef.add(ak1);
    aFb_mCoef.add(ak2);
    aFb_mCoef.add(ak3);
    aFb_mCoef.add(ak4);
    
    RooArgList Fb_mCoefX;
    Fb_mCoefX.add(cos_m);
    Fb_mCoefX.add(k1X);
    Fb_mCoefX.add(k2X);
    Fb_mCoefX.add(k3X);
    Fb_mCoefX.add(k4X);
    
    RooArgList aFb_mCoefX;
    aFb_mCoefX.add(cos_m);
    aFb_mCoefX.add(ak1X);
    aFb_mCoefX.add(ak2X);
    aFb_mCoefX.add(ak3X);
    aFb_mCoefX.add(ak4X);
    
    
    RooChebychev Fb_l("Fb_l","Background_l",cos_l,Fb_lCoef);
    RooChebychev aLFb_l("aLFb_l","Background_l",cos_l,aFb_lCoef);
    RooChebychev Fb_lX("Fb_lX","Background_l",cos_l,Fb_lXCoef);
    RooChebychev aLFb_lX("aLFb_lX","Background_l",cos_l,aFb_lXCoef);
    
    RooChebychev Fb_p("Fb_p","Background_p",cos_p,Fb_pCoef);
    RooChebychev aLFb_p("aLFb_p","Background_p",cos_p,aLFb_pCoef);
    RooChebychev Fb_pX("Fb_pX","Background_p",cos_p,Fb_pCoefX);
    RooChebychev aLFb_pX("aLFb_pX","Background_p",cos_p,aLFb_pCoefX);
    
    RooGenericPdf Fb_mu("Fb_mu","","(1+TMath::Erf((@0-@1)/@2))*(1+TMath::Erf((-@0-@3)/@4))",Fb_mCoef);
    RooGenericPdf aLFb_mu("aLFb_mu","","(1+TMath::Erf((@0-@1)/@2))*(1+TMath::Erf((-@0-@3)/@4))",aFb_mCoef);
    RooGenericPdf Fb_muX("Fb_muX","","(1+TMath::Erf((@0-@1)/@2))*(1+TMath::Erf((-@0-@3)/@4))",Fb_mCoefX);
    RooGenericPdf aLFb_muX("aLFb_muX","","(1+TMath::Erf((@0-@1)/@2))*(1+TMath::Erf((-@0-@3)/@4))",aFb_mCoefX);
    
    RooProdPdf PDFB_ang("PDFB_ang","PDF Bkg",RooArgList(Fb_l,Fb_p,Fb_mu));
    RooProdPdf aLPDFB_ang("aLPDFB_ang","PDF Bkg antiLambda",RooArgList(aLFb_l,aLFb_p,aLFb_mu));
    RooProdPdf PDFBX_ang("PDFBX_ang","PDF Bkg 2012",RooArgList(Fb_lX,Fb_pX,Fb_muX));
    RooProdPdf aLPDFBX_ang("aLPDFBX_ang","PDF Bkg antiLambda 2012",RooArgList(aLFb_lX,aLFb_pX,aLFb_muX));
    
    RooArgSet *lambBkgParam = PDFB_ang.getParameters(RooArgSet(cos_l,cos_p,cos_m));
    RooArgSet *alambBkgParam = aLPDFB_ang.getParameters(RooArgSet(cos_l,cos_p,cos_m));
    RooArgSet *lambBkgParamX = PDFBX_ang.getParameters(RooArgSet(cos_l,cos_p,cos_m));
    RooArgSet *alambBkgParamX = aLPDFBX_ang.getParameters(RooArgSet(cos_l,cos_p,cos_m));
    
    lambBkgParam->readFromFile("textFiles/LambdaBkgOnce_w.txt");
    alambBkgParam->readFromFile("textFiles/antiLambdaBkgOnce_w.txt");
    lambBkgParamX->readFromFile("textFiles/LambdaBkgDoce_w.txt");
    alambBkgParamX->readFromFile("textFiles/antiLambdaBkgDoce_w.txt");
    
    /*lambBkgParam->Print("v");
    alambBkgParam->Print("v");
    lambBkgParamX->Print("v");
    alambBkgParamX->Print("v");*/
    
    
    RooProdPdf PDFB("PDFB","PDF Bkg",RooArgList(PDFB_ang, poly));
    RooProdPdf aLPDFB("aLPDFB","PDF Bkg antiLambda",RooArgList(aLPDFB_ang,poly_a));
    RooProdPdf PDFBX("PDFBX","PDF Bkg 2012",RooArgList(PDFBX_ang,poly_X));
    RooProdPdf aLPDFBX("aLPDFBX","PDF Bkg antiLambda 2012",RooArgList(aLPDFBX_ang,poly_a_X));
    
    // Total PDF's
    
    RooProdPdf PDFSnoeff("PDFSnoeff","Mass and angular distrib",RooArgList(FSig,DGauss));
    RooProdPdf aLPDFSnoeff("aLPDFSnoeff","Mass and angular distrib antiLambda",RooArgList(FSig,DGauss_a));
    RooProdPdf PDFSnoeffX("PDFSnoeffX","Mass and angular distrib",RooArgList(FSig,DGauss_X));
    RooProdPdf aLPDFSnoeffX("aLPDFSnoeffX","Mass and angular distrib antiLambda",RooArgList(FSig,DGauss_a_X));
    
    RooEffProd Fsigeff("Fsigeff","angular distrib with eff",FSig,eps);
    RooEffProd aLFsigeff("aLFsigeff","angular distrib with eff antiLambda",FSig,aLeps);//PDF señal antiLambda
    RooEffProd FsigeffX("FsigeffX","angular distrib with eff",FSig,epsX);
    RooEffProd aLFsigeffX("aLFsigeffX","angular distrib with eff antiLambda",FSig,aLepsX);//PDF señal antiLambda
    
    //opción 1 con eficiencia
    //RooEffProd PDFS("PDFS","Mass and angular distrib and eff",PDFSnoeff,eps);
    //opción 2 con eficiencia
    RooProdPdf PDFS("PDFS","Mass and angular distrib and eff",RooArgList(Fsigeff,DGauss));
    RooProdPdf aLPDFS("aLPDFS","Mass and angular distrib and eff antiLambda",RooArgList(aLFsigeff,DGauss_a));
    RooProdPdf PDFSX("PDFSX","Mass and angular distrib and eff",RooArgList(FsigeffX,DGauss_X));
    RooProdPdf aLPDFSX("aLPDFSX","Mass and angular distrib and eff antiLambda",RooArgList(aLFsigeffX,DGauss_a_X));
    
    RooArgList PDFList;
    PDFList.add(PDFS);
    PDFList.add(PDFB);
    
    RooArgList aLPDFList;
    aLPDFList.add(aLPDFS);
    aLPDFList.add(aLPDFB);
    
    RooArgList PDFListX;
    PDFListX.add(PDFSX);
    PDFListX.add(PDFBX);
    
    RooArgList aLPDFListX;
    aLPDFListX.add(aLPDFSX);
    aLPDFListX.add(aLPDFBX);
    
    
    RooAddPdf PDFT("PDFT","Total PDF lambda",PDFList,NList);
    RooAddPdf aPDFT("aPDFT","Total PDF antiLambda",aLPDFList,NList_a);
    RooAddPdf PDFTX("PDFTX","Total PDF Lambda 2012",PDFListX,NList_X);
    RooAddPdf aPDFTX("aPDFTX","Total PDF antiLambda 2012",aLPDFListX,NList_a_X);
    
    /* --------- Simultaneous -------------*/
    
    // Define category to distinguish physics and control samples events
    RooCategory sample("sample","sample") ;
    sample.defineType("Lambda") ;
    sample.defineType("AntiLambda") ;
    sample.defineType("LambdaX") ;
    sample.defineType("AntiLambdaX") ;
    
    // C o n s t r u c t  s i m u l t a n e o u s   p d f s   i n   ( P() , s a m p l e )
    // -----------------------------------------------------------------------------------
    RooSimultaneous simPDFT_mass("simPDFT_mass","simPDFT_mass",sample);
    
    simPDFT_mass.addPdf(PDFT_mass,"Lambda");
    simPDFT_mass.addPdf(PDFT_mass_a,"AntiLambda");
    simPDFT_mass.addPdf(PDFT_mass_X,"LambdaX");
    simPDFT_mass.addPdf(PDFT_mass_a_X,"AntiLambdaX");
    
    // Simultaneous PDF bkg
    RooSimultaneous simPDFbkg("simPDFbkg","simPDFbkg",sample);
    
    simPDFbkg.addPdf(PDFB,"Lambda");
    simPDFbkg.addPdf(aLPDFB,"AntiLambda");
    simPDFbkg.addPdf(PDFBX,"LambdaX");
    simPDFbkg.addPdf(aLPDFBX,"AntiLambdaX");
    
    RooSimultaneous simPDFT("simPDFT","simPDFT",sample);
    
    simPDFT.addPdf(PDFT,"Lambda");
    simPDFT.addPdf(aPDFT,"AntiLambda");
    simPDFT.addPdf(PDFTX,"LambdaX");
    simPDFT.addPdf(aPDFTX,"AntiLambdaX");
    
    RooSimultaneous simPDFS("simPDFS","simPDFS",sample);
    
    simPDFS.addPdf(PDFS,"Lambda");
    simPDFS.addPdf(aLPDFS,"AntiLambda");
    simPDFS.addPdf(PDFS,"LambdaX");
    simPDFS.addPdf(aLPDFS,"AntiLambdaX");
    
    RooArgSet* massParam = simPDFT_mass.getParameters(RooArgSet(mass));
    
    
    massParam->readFromFile("textFiles/massParameters.txt");
    massParam->Print("v");
 
    mean.setConstant(kTRUE);
    N_bkg.setConstant(kTRUE);   N_bkg_a.setConstant(kTRUE);
    N_bkg_X.setConstant(kTRUE); N_bkg_a_X.setConstant(kTRUE);
    N_sig.setConstant(kTRUE);   N_sig_a.setConstant(kTRUE);
    N_sig_X.setConstant(kTRUE); N_sig_a_X.setConstant(kTRUE);
    s1.setConstant(kTRUE);      s1_a.setConstant(kTRUE);
    s1_X.setConstant(kTRUE);    s1_a_X.setConstant(kTRUE);
    s2.setConstant(kTRUE);      s2_a.setConstant(kTRUE);
    s2_X.setConstant(kTRUE);    s2_a_X.setConstant(kTRUE);
    f.setConstant(kTRUE);       f_a.setConstant(kTRUE);
    f_X.setConstant(kTRUE);     f_a_X.setConstant(kTRUE);
    p_.setConstant(kTRUE);      p_a.setConstant(kTRUE);
    p_X.setConstant(kTRUE);     p_a_X.setConstant(kTRUE);
    
    //massParam->Print("v"); return;
    
    mass.setRange("R3",Msigmamin,Msigmamax);
    
    RooAbsReal* fracInt = poly.createIntegral(mass,NormSet(mass),Range("R3"));
    //std::cout << fracInt->getVal() << std::endl;
    RooAbsReal* fracInt_a = poly_a.createIntegral(mass,NormSet(mass),Range("R3"));
    //std::cout << fracInt_a->getVal() << std::endl;
    RooAbsReal* fracInt_X = poly_X.createIntegral(mass,NormSet(mass),Range("R3"));
    //std::cout << fracInt_X->getVal() << std::endl;
    RooAbsReal* fracInt_a_X = poly_a_X.createIntegral(mass,NormSet(mass),Range("R3"));
    //std::cout << fracInt_a_X->getVal() << std::endl;
    
    N_bkg.setVal(N_bkg.getVal()*(fracInt->getVal()));
    N_bkg_a.setVal(N_bkg_a.getVal()*(fracInt_a->getVal()));
    N_bkg_X.setVal(N_bkg_X.getVal()*(fracInt_X->getVal()));
    N_bkg_a_X.setVal(N_bkg_a_X.getVal()*(fracInt_a_X->getVal()));
   
    Double_t minNll, edmFit, pol, polError, polPull;
    Double_t alp1,alp1Error, alp1Pull;
    Double_t alp2, alp2Error, alp2Pull, gam, gamError, gamPull;
    Double_t polErrorLo, polErrorHi, alp1ErrorLo, alp1ErrorHi;
    Double_t alp2ErrorLo, alp2ErrorHi, gamErrorLo, gamErrorHi;
    Int_t stat,covqual,badNll;

    Double_t initP, initAlp1, initAlp2, initGamma;
    initP = 0.00;
    initAlp1 = 0.12;
    initAlp2 = -0.93;
    initGamma = -0.46;

    TString masspeakCut = "mass>"; masspeakCut += Msigmamin; masspeakCut += " &&"; masspeakCut += " mass<"; masspeakCut += Msigmamax;
    mass.setRange("senial",Msigmamin,Msigmamax);
    
    P.setVal(initP);
    alpha_one.setVal(initAlp1);
    alpha_two.setVal(initAlp2);
    gamma.setVal(initGamma);

    RooRandom::randomGenerator()->SetSeed(0);

    RooDataSet* dataLb=PDFT.generate(RooArgSet(mass,cos_l,cos_p,cos_m),N_bkg.getVal()+N_sig.getVal());
    dataLb->Print();
    RooDataSet* dataALb=aPDFT.generate(RooArgSet(mass,cos_l,cos_p,cos_m),N_bkg_a.getVal()+N_sig_a.getVal());
    dataALb->Print();
    RooDataSet* dataLbX=PDFTX.generate(RooArgSet(mass,cos_l,cos_p,cos_m),N_bkg_X.getVal()+N_sig_X.getVal());
    dataLbX->Print();
    RooDataSet* dataALbX=aPDFTX.generate(RooArgSet(mass,cos_l,cos_p,cos_m),N_bkg_a_X.getVal()+N_sig_a_X.getVal());
    dataALbX->Print();

    RooDataSet combData("combData","combined data",RooArgSet(cos_l,cos_m,cos_p,mass),Index(sample),Import("Lambda",*dataLb),Import("AntiLambda",*dataALb),Import("LambdaX",*dataLbX),Import("AntiLambdaX",*dataALbX)) ;

    P.setVal(0);
    alpha_one.setVal(0);
    alpha_two.setVal(0);
    gamma.setVal(-0.5);

    RooAbsReal* nll = simPDFT.createNLL(combData,Extended(kTRUE)/*,ExternalConstraints(extconstraints)*/,NumCPU(cpu)) ;
    RooMinuit m(*nll);
    m.migrad();
    RooFitResult* r_migrad = m.save(); if((Int_t)r_migrad->status()!=0){return;}
    m.hesse();
    RooFitResult* r_hesse = m.save(); if((Int_t)r_hesse->status()!=0){return;}
    m.minos();
    RooFitResult* r = m.save();
    //r->Print("v") ;
    if((Int_t)r->status()!=0) {return;}
    if((Int_t)r->covQual()!=3) {return;}
    if((Double_t)r->edm() > 0.001) {return;}

    std::cout << "P: " << P.getValV() << " +/- " << P.getError() << std::endl;
    std::cout << "alpha1: " << alpha_one.getValV() << " +/- " << alpha_one.getError() << std::endl;
    std::cout << "alpha2: " << alpha_two.getValV() << " +/- " << alpha_two.getError() << std::endl;
    std::cout << "gamma: " << gamma.getValV() << " +/- " << gamma.getError() << std::endl;
    
    TFile *f1 = new TFile("ToyMC.root","RECREATE");
    TTree *t = new TTree("toyMC","Tree");
    f1->cd();
    
    t->Branch("P",&pol);
    t->Branch("error_P",&polError);
    t->Branch("error_P_Lo",&polErrorLo);
    t->Branch("error_P_Hi",&polErrorHi);
    t->Branch("alp1",&alp1);
    t->Branch("alp1_Error",&alp1Error);
    t->Branch("alp1_ErrorLo",&alp1ErrorLo);
    t->Branch("alp1_ErrorHi",&alp1ErrorHi);
    t->Branch("alp2",&alp2);
    t->Branch("alp2_Error",&alp2Error);
    t->Branch("alp2_ErrorLo",&alp2ErrorLo);
    t->Branch("alp2_ErrorHi",&alp2ErrorHi);
    t->Branch("gamma",&gam);
    t->Branch("gamma_Error",&gamError);
    t->Branch("gamma_ErrorLo",&gamErrorLo);
    t->Branch("gamma_ErrorHi",&gamErrorHi);
    t->Branch("pull_P",&polPull);
    t->Branch("pull_alp1",&alp1Pull);
    t->Branch("pull_alp2",&alp2Pull);
    t->Branch("pull_gam",&gamPull);
    t->Branch("nll",&minNll);
    t->Branch("edm",&edmFit);
    t->Branch("badNll",&badNll);
    t->Branch("status",&stat);
    t->Branch("covqual",&covqual);

    pol = P.getVal();                    polError = P.getError();
    alp1 = alpha_one.getVal();           alp1Error = alpha_one.getError();
    alp2 = alpha_two.getVal();           alp2Error = alpha_two.getError();
    gam = gamma.getVal();                gamError = gamma.getError();
    
    polPull =  ( pol - initP)/polError;
    alp1Pull = ( alp1 - initAlp1)/alp1Error;
    alp2Pull = ( alp2 - initAlp2)/alp2Error;
    gamPull = (gam - initGamma)/gamError;
    
    polErrorLo = P.getAsymErrorLo();      
    polErrorHi = P.getAsymErrorHi();
    alp1ErrorLo = alpha_one.getAsymErrorLo();
    alp1ErrorHi = alpha_one.getAsymErrorHi();
    alp2ErrorLo = alpha_two.getAsymErrorLo();
    alp2ErrorHi = alpha_two.getAsymErrorHi();
    gamErrorLo = gamma.getAsymErrorLo();
    gamErrorHi = gamma.getAsymErrorHi();

    minNll = r->minNll();
    edmFit = r->edm();
    badNll = r->numInvalidNLL();

    stat = r->status();
    covqual = r->covQual();

    t->Fill();
    
    //delete combData;
    delete dataLb;
    delete dataLbX;
    delete dataALb;
    delete dataALbX;
    delete r_migrad;
    delete r_hesse;
    delete r;
    
    f1->Write("",TObject::kOverwrite);

}
