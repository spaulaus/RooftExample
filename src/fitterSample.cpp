/** \name fitterSample.cpp 
 *  
 *  A sample file for using Roofit to do maximum likelyhood fitting.
 *  This file is accompanied by Makefile and ML.pdf
 *  The code is distributed under the GNU GPL 3.0
 *  
 *  \author G. Cerizza
 *  \date 30 Nov. 2012
 *  \modified S. V. Paulauskas 
 *  \last modified 02 June 2013
 */
#include <fstream>
#include <iostream>

#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooProdPdf.h"

#include "TAxis.h"
#include "TCanvas.h"

using namespace RooFit;
using namespace std;

void fitter_sample();

int main() {
    fitter_sample();
}

void fitter_sample()
{
  // list of variables 
  RooRealVar genergy("genergy","genergy",1000,2000);

  // loading datafile from root file 
  //  RooDataSet *data;
  //  TFile *f = new TFile("file.root");
  //  TTree* tree = (TTree*)f->Get("tree");
  //  data = new RooDataSet("data","data",tree,RooArgList(genergy));
  // or txt file  
  //  RooDataSet* data = RooDataSet::read("file.txt",RooArgList(genergy));

  ///////////////////////////////////////////
  // for a parameterized variable
  ///////////////////////////////////////////
  //RooRealVar mu("mu","", 30., -10., 10.);
  //RooFormulaVar sigma("sigma", "0.0264412*mu0+0.0432495", mu0);
 
  ///////////////////////////////////////////
  // model for gamma peaks and background
  ///////////////////////////////////////////
  // PDFs parameters for peak 1
  RooRealVar sig_m1("sig_m1","",1350,1300,1400);
  RooRealVar sig_w1("sig_w1","",10,0,20);
  RooGaussian peak_1("peak_1","peak distribution",genergy,sig_m1,sig_w1);  
  // PDFs parameters for peak 2
  RooRealVar sig_m2("sig_m2","",1420,1380,1450);
  RooRealVar sig_w2("sig_w2","",15,0,20);
  RooGaussian peak_2("peak_2","peak distribution",genergy,sig_m2,sig_w2);  
  // PDFs parameters for peak 3
  RooRealVar sig_m3("sig_m3","",1850,1750,1950);
  RooRealVar sig_w3("sig_w3","",15,0,20);
  RooGaussian peak_3("peak_3","peak distribution",genergy,sig_m3,sig_w3);
  // PDFs parameters for bkg
  RooRealVar a0("a0","a0",-0.3,-1,1);
  RooRealVar a1("a1","a1",0.1,-3,3);
  RooChebychev bkg("bkg","flat background",genergy,RooArgList(a0,a1));

  /////////////////////////////////
  // yields
  /////////////////////////////////
  RooRealVar nsig1("nsig1","number of events in peak 1",15000,0.,100000) ; // corresponding to peak1
  RooRealVar nsig2("nsig2","number of events in peak 2",10000,0.,100000) ; // corresponding to peak2
  RooRealVar nsig3("nsig3","number of events in peak 3",5000,0.,100000) ; // corresponding to peak3
  RooRealVar nbkg("nbkg","number of background events",70000,0,100000) ; // corresponding to bgk

  /////////////////////////////////
  // add all the pdfs together
  /////////////////////////////////
  RooAddPdf  model("model","model",RooArgList(peak_1,peak_2,peak_3,bkg),RooArgList(nsig1,nsig2,nsig3,nbkg)) ;
  
  // Generation of 100k MonteCarlo events
  RooDataSet* dataMC = model.generate(genergy,100000);
  
  ////////////////////////////////////////
  // fit the added pdf to the data set
  // parameters are fitTo(<dataSet>, number of cpus for calc, 
  //                      save the file, fit range)
  ////////////////////////////////////////
  RooFitResult* fitResult = model.fitTo(*dataMC, NumCPU(3), Save(), 
                                        Range(1000., 2000.));
  ofstream resultsParam("results/fitResults.fit");
  fitResult->printMultiline(resultsParam, 0, false, "");
  resultsParam.close();
  
  // gamma energy plot
  RooPlot* frame = genergy.frame();
  frame = genergy.frame(100);
  frame->SetTitle("Gamma spectra");
  frame->SetXTitle("energy (keV)");
  frame->SetYTitle("Events/(100 keV)");
  frame->GetYaxis()->SetTitleOffset(1.2);
  dataMC->plotOn(frame,Name("data"));
  model.plotOn(frame,Name("model"));
  model.plotOn(frame,RooFit::Components("peak_1,nsig1"),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));  
  model.plotOn(frame,RooFit::Components("peak_2,nsig2"),RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));  
  model.plotOn(frame,RooFit::Components("peak_3,nsig3"),RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));  
  model.plotOn(frame,RooFit::Components("bkg,nbkg"),RooFit::LineColor(29), RooFit::LineStyle(kDashed));  

  TCanvas* c = new TCanvas("c","",0,0,700,500) ;
  c->cd();
  c->SetFillColor(kWhite);
  frame->Draw();
  c->SaveAs("results/results.eps");
}
