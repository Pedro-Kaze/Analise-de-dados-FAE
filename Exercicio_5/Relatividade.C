#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <filesystem>
#include <algorithm> 
#include <string>    

double calcular_massa_invariante(const std::vector<float>& pt, const std::vector<float>& eta, const std::vector<float>& phi) {
    if (pt.size() >= 2) {
        return sqrt(2 * pt[0] * pt[1] * (TMath::CosH(eta[0] - eta[1]) - TMath::Cos(phi[0] - phi[1])));
    }
    return -1.0; 
}

void Relatividade() {
    std::string diretorio = "/opendata/eos/opendata/cms/Run2016G/DoubleEG/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/250000";

    std::vector<double> e_massas_invariantes;

    // Histogramas
    TH1F* hElectronPt = new TH1F("hElectronPt", "Electron p_{T} Distribution", 50, 0, 200);
    TH1F* hMuonPt = new TH1F("hMuonPt", "Muon p_{T} Distribution", 50, 0, 200);
    TH1F* hTauPt = new TH1F("hTauPt", "Tau p_{T} Distribution", 50, 0, 200);
    TH1F* hJetPt = new TH1F("hJetPt", "Jet p_{T} Distribution", 50, 0, 200);

    TH1F* hElectronEta = new TH1F("hElectronEta", "Electron #eta Distribution", 50, -5, 5);
    TH1F* hMuonEta = new TH1F("hMuonEta", "Muon #eta Distribution", 50, -5, 5);
    TH1F* hTauEta = new TH1F("hTauEta", "Tau #eta Distribution", 50, -5, 5);
    TH1F* hJetEta = new TH1F("hJetEta", "Jet #eta Distribution", 50, -5, 5);

    TH1F* hElectronPhi = new TH1F("hElectronPhi", "Electron #phi Distribution", 50, -3.14, 3.14);
    TH1F* hMuonPhi = new TH1F("hMuonPhi", "Muon #phi Distribution", 50, -3.14, 3.14);
    TH1F* hTauPhi = new TH1F("hTauPhi", "Tau #phi Distribution", 50, -3.14, 3.14);
    TH1F* hJetPhi = new TH1F("hJetPhi", "Jet #phi Distribution", 50, -3.14, 3.14);

    for (const auto& entry : std::filesystem::directory_iterator(diretorio)) {
        std::string file_path = entry.path();
        TFile file(file_path.c_str(), "READ");
        if (!file.IsOpen()) continue;

        TTreeReader reader("Events", &file);
        TTreeReaderArray<float> Electron_pt(reader, "Electron_pt");
        TTreeReaderArray<float> Electron_eta(reader, "Electron_eta");
        TTreeReaderArray<float> Electron_phi(reader, "Electron_phi");
        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
        TTreeReaderArray<float> Jet_pt(reader, "Jet_pt");
        TTreeReaderArray<float> Jet_eta(reader, "Jet_eta");
        TTreeReaderArray<float> Jet_phi(reader, "Jet_phi");
        TTreeReaderArray<float> Tau_pt(reader, "Tau_pt");
        TTreeReaderArray<float> Tau_eta(reader, "Tau_eta");
        TTreeReaderArray<float> Tau_phi(reader, "Tau_phi");

        while (reader.Next()) {
            // Preencher histogramas de pT, η e φ
            for (size_t i = 0; i < Electron_pt.GetSize(); ++i) {
                hElectronPt->Fill(Electron_pt[i]);
                hElectronEta->Fill(Electron_eta[i]);
                hElectronPhi->Fill(Electron_phi[i]);
            }
            for (size_t i = 0; i < Muon_pt.GetSize(); ++i) {
                hMuonPt->Fill(Muon_pt[i]);
                hMuonEta->Fill(Muon_eta[i]);
                hMuonPhi->Fill(Muon_phi[i]);
            }
            for (size_t i = 0; i < Jet_pt.GetSize(); ++i) {
                hJetPt->Fill(Jet_pt[i]);
                hJetEta->Fill(Jet_eta[i]);
                hJetPhi->Fill(Jet_phi[i]);
            }
            for (size_t i = 0; i < Tau_pt.GetSize(); ++i) {
                hTauPt->Fill(Tau_pt[i]);
                hTauEta->Fill(Tau_eta[i]);
                hTauPhi->Fill(Tau_phi[i]);
            }

            // Calcular a massa invariante dos dois léptons com maior pT
            std::vector<std::pair<float, int>> leptons; // (pT, índice)
            for (size_t i = 0; i < Electron_pt.GetSize(); ++i) {
                leptons.emplace_back(Electron_pt[i], i); // Elétrons
            }
            for (size_t i = 0; i < Muon_pt.GetSize(); ++i) {
                leptons.emplace_back(Muon_pt[i], i + Electron_pt.GetSize()); // Muons
            }
            for (size_t i = 0; i < Tau_pt.GetSize(); ++i) {
                leptons.emplace_back(Tau_pt[i], i + Electron_pt.GetSize() + Muon_pt.GetSize()); // Taus
            }

            // Ordenar  pT
            std::sort(leptons.rbegin(), leptons.rend());

            if (leptons.size() >= 2) {
                int idx1 = leptons[0].second;
                int idx2 = leptons[1].second;

                float pt1, eta1, phi1, pt2, eta2, phi2;

                if (idx1 < Electron_pt.GetSize()) {
                    pt1 = Electron_pt[idx1];
                    eta1 = Electron_eta[idx1];
                    phi1 = Electron_phi[idx1];
                } else if (idx1 < Electron_pt.GetSize() + Muon_pt.GetSize()) {
                    pt1 = Muon_pt[idx1 - Electron_pt.GetSize()];
                    eta1 = Muon_eta[idx1 - Electron_pt.GetSize()];
                    phi1 = Muon_phi[idx1 - Electron_pt.GetSize()];
                } else {
                    pt1 = Tau_pt[idx1 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                    eta1 = Tau_eta[idx1 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                    phi1 = Tau_phi[idx1 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                }

                if (idx2 < Electron_pt.GetSize()) {
                    pt2 = Electron_pt[idx2];
                    eta2 = Electron_eta[idx2];
                    phi2 = Electron_phi[idx2];
                } else if (idx2 < Electron_pt.GetSize() + Muon_pt.GetSize()) {
                    pt2 = Muon_pt[idx2 - Electron_pt.GetSize()];
                    eta2 = Muon_eta[idx2 - Electron_pt.GetSize()];
                    phi2 = Muon_phi[idx2 - Electron_pt.GetSize()];
                } else {
                    pt2 = Tau_pt[idx2 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                    eta2 = Tau_eta[idx2 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                    phi2 = Tau_phi[idx2 - Electron_pt.GetSize() - Muon_pt.GetSize()];
                }

                // Calcular massa invariante
                std::vector<float> pt_values = {pt1, pt2};
                std::vector<float> eta_values = {eta1, eta2};
                std::vector<float> phi_values = {phi1, phi2};
                double massa_invariante = calcular_massa_invariante(pt_values, eta_values, phi_values);

                if (massa_invariante >= 0) {
                    e_massas_invariantes.push_back(massa_invariante);
                }
            }
        }
    }

    
    TCanvas* canvas;

    // Salvar histogramas para pT
    canvas = new TCanvas("canvasElectronPt", "Electron p_{T} Distribution", 800, 600);
    hElectronPt->SetLineColor(kBlue);
    hElectronPt->Draw();
    hElectronPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hElectronPt->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("electron_pt_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasMuonPt", "Muon p_{T} Distribution", 800, 600);
    hMuonPt->SetLineColor(kRed);
    hMuonPt->Draw();
    hMuonPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hMuonPt->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("muon_pt_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasTauPt", "Tau p_{T} Distribution", 800, 600);
    hTauPt->SetLineColor(kMagenta);
    hTauPt->Draw();
    hTauPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hTauPt->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("tau_pt_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasJetPt", "Jet p_{T} Distribution", 800, 600);
    hJetPt->SetLineColor(kGreen);
    hJetPt->Draw();
    hJetPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hJetPt->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("jet_pt_distribution.png");
    delete canvas;

    // Salvar histogramas para η
    canvas = new TCanvas("canvasElectronEta", "Electron #eta Distribution", 800, 600);
    hElectronEta->SetLineColor(kBlue);
    hElectronEta->Draw();
    hElectronEta->GetXaxis()->SetTitle("#eta");
    hElectronEta->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("electron_eta_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasMuonEta", "Muon #eta Distribution", 800, 600);
    hMuonEta->SetLineColor(kRed);
    hMuonEta->Draw();
    hMuonEta->GetXaxis()->SetTitle("#eta");
    hMuonEta->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("muon_eta_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasTauEta", "Tau #eta Distribution", 800, 600);
    hTauEta->SetLineColor(kMagenta);
    hTauEta->Draw();
    hTauEta->GetXaxis()->SetTitle("#eta");
    hTauEta->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("tau_eta_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasJetEta", "Jet #eta Distribution", 800, 600);
    hJetEta->SetLineColor(kGreen);
    hJetEta->Draw();
    hJetEta->GetXaxis()->SetTitle("#eta");
    hJetEta->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("jet_eta_distribution.png");
    delete canvas;

    // Salvar histogramas para φ
    canvas = new TCanvas("canvasElectronPhi", "Electron #phi Distribution", 800, 600);
    hElectronPhi->SetLineColor(kBlue);
    hElectronPhi->Draw();
    hElectronPhi->GetXaxis()->SetTitle("#phi");
    hElectronPhi->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("electron_phi_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasMuonPhi", "Muon #phi Distribution", 800, 600);
    hMuonPhi->SetLineColor(kRed);
    hMuonPhi->Draw();
    hMuonPhi->GetXaxis()->SetTitle("#phi");
    hMuonPhi->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("muon_phi_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasTauPhi", "Tau #phi Distribution", 800, 600);
    hTauPhi->SetLineColor(kMagenta);
    hTauPhi->Draw();
    hTauPhi->GetXaxis()->SetTitle("#phi");
    hTauPhi->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("tau_phi_distribution.png");
    delete canvas;

    canvas = new TCanvas("canvasJetPhi", "Jet #phi Distribution", 800, 600);
    hJetPhi->SetLineColor(kGreen);
    hJetPhi->Draw();
    hJetPhi->GetXaxis()->SetTitle("#phi");
    hJetPhi->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("jet_phi_distribution.png");
    delete canvas;

    // plot da Massa Invariante
    TH1F* hMassaInvariante = new TH1F("hMassaInvariante", "Invariant Mass Distribution", 50, 0, 200);
    for (const auto& massa : e_massas_invariantes) {
        if (massa >= 0) hMassaInvariante->Fill(massa);
    }

    canvas = new TCanvas("canvasInvariantMass", "Invariant Mass Distribution", 800, 600);
    hMassaInvariante->SetLineColor(kBlack);
    canvas->SetLogy();
    hMassaInvariante->Draw();
    hMassaInvariante->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    hMassaInvariante->GetYaxis()->SetTitle("Events");
    canvas->SaveAs("invariant_mass_distribution.png");
    delete canvas;

   
    delete hElectronPt;
    delete hMuonPt;
    delete hTauPt;
    delete hJetPt;
    delete hElectronEta;
    delete hMuonEta;
    delete hTauEta;
    delete hJetEta;
    delete hElectronPhi;
    delete hMuonPhi;
    delete hTauPhi;
    delete hJetPhi;
    delete hMassaInvariante;
}
