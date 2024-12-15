#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TH1.h>

void SinalBackground() {
    
    TFile *file = TFile::Open("DataSet_lowstat.root");
    TH1F* hist = (TH1F*)file->Get("data");
    RooDataSet *data = (RooDataSet*)file->Get("data");
    
    RooRealVar mass("mass", "mass", 2, 4);  

    //sinal
    RooRealVar mean("mean", "mean", 3.1, 3.05, 3.15);
    RooRealVar sigma("sigma", "sigma", 0.1, 0.01, 0.15);
    RooRealVar alpha("alpha", "alpha", 1.5, 0.5, 5);
    RooRealVar N("N", "N", 2, -10, 10);
    RooCBShape signal("signal", "signal", mass, mean, sigma, alpha, N);

    // background
    RooRealVar lambda("lambda", "lambda", 5, 2, 7);
    RooExponential background("background", "background", mass, lambda);
    RooPolynomial pol("pol","pol", mass);

   
    RooRealVar frac("frac","frac",0.9,0.8,1);
    RooAddPdf model("model", "model", RooArgList(signal, pol), RooArgList(frac));
    //RooDataHist *data("data", "data", RooArgList(mass), hist);

    model.fitTo(*data);

    RooPlot* frame = mass.frame();
    data->plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, RooFit::Components("pol"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));  
    model.plotOn(frame, RooFit::Components("signal"), RooFit::LineStyle(kDotted), RooFit::LineColor(kGreen));   

    model.paramOn(frame, RooFit::Layout(0.1, 0.9, 0.9));

    double chi2 = frame->chiSquare();
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    frame->Draw();
    TLegend *leg = new TLegend(-0.1, 0.2, 0.9, 0.9);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f", chi2), "");
    leg->Draw();
    std::cout << "Chi^2 / ndf = " << chi2 << std::endl;
    c->SaveAs("ex3.pdf");

}
