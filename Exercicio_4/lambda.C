#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "TH1.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooExponential.h"



void lambda(){
    RooRealVar x("x","x",0,10);
    RooRealVar lambda("lambda","lambda",-1,-2,-0.1);
    RooExponential expo("expo", "Decaimento de Lambda",x,lambda);
    RooRealVar N("N", "Numero de eventos", 1500, 0, 2000);

    RooDataSet* data = expo.generate(RooArgSet(x), 1500);

    RooExtendPdf extExp("extExp", "Extended", expo, N);

    extExp.fitTo(*data ,RooFit::Save());

    RooPlot* frame = x.frame();
    data->plotOn(frame);
    extExp.plotOn(frame);

    RooFitResult* fitResult = extExp.fitTo(*data, RooFit::Save());

    
    TCanvas *c =new TCanvas("c","Ajuste Exponencial",800,600);

    frame->Draw();


    double chi2 = frame->chiSquare();
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry((TObject*)0, Form("Lambda = %.3f #pm %.3f",lambda.getVal(), lambda.getError()), "");
    leg->AddEntry((TObject*)0, Form("Eventos ajustados = %.1f #pm %.1f",N.getVal(),N.getError()), "");
    leg->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f", chi2), "");    
    leg->Draw();
    
    
    c->Draw();
    c->SaveAs("Lambda.png");
    
}
