#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooFit.h"
#include "TObject.h"
#include "TMath.h"


void generateCrystalBall() {

// Definir as variaveis x, mean e width com seus intervalos

RooRealVar cbmean("cbmean", "cbmean" , 90.0, 20, 180.0) ;
RooRealVar cbsigma("cbsigma", "cbsigma" , 2.0, 1.0, 40.0) ;
RooRealVar cbsig("cbsig", "cbsignal", 10, 0, 1000000);
RooRealVar n("n","", 5.1);
RooRealVar alpha("alpha","", 1.3);
RooRealVar x("x","x",0,80,100);
RooCBShape cball("cball", "crystal ball", x, cbmean, cbsigma, alpha, n);


RooDataSet* data = cball.generate(RooArgSet(x), 10000);
RooFitResult* fit_result = cball.fitTo(*data,RooFit::Save());
fit_result->Print("v");

TCanvas* c1 = new TCanvas("c1", "Crystal Ball Event Generation", 800,600);
    
//Criar um objeto Rooplot que serve como um frame para o grafico
RooPlot* frame = x.frame();

//Plotar a PDF no frame (a função gaussiana g)
data->plotOn(frame);
cball.plotOn(frame);

frame->Draw();
    

double chi2 = frame->chiSquare();
TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
leg->SetTextSize(0.03);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry((TObject*)0, Form("Mean = %.3f #pm %.3f",cbmean.getVal(), cbmean.getError()), "");
leg->AddEntry((TObject*)0, Form("Sigma = %.3f #pm %.3f",cbsigma.getVal(), cbsigma.getError()), "");
leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f", chi2), "");    
leg->Draw();


c1->Draw();

c1->SaveAs("CrystalBallEvents.pdf");
}
