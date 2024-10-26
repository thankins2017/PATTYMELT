// Checks gain match accuracy for both front and back sets of signals.

// C++ includes
#include <fstream>

// ROOT includes
#include "TCanvas.h"
#include "TApplication.h"
#include "TRootCanvas.h"

// Other includes
#include "PathManager.h"
#include "AnalysisTools.h"

int main(int argc, char **argv) {
    // option manager
    BrAppOptionManager::Instance()->SetVersion(1, 0, "check_gain_match");
    BrAppOptionManager::Instance()->SetHelp("Reads in calib data and checks gain match pars.");
    BrAppOptionManager::Instance()->SetCommandLine(argc, argv);
    auto *opt_run_number = new BrAppIntOption('r', "run_number", "run number", CALIB_RUN);
    auto *opt_events     = new BrAppIntOption('e', "events", "events to analyze", std::numeric_limits<int>::max());
    auto *opt_display    = new BrAppBoolOption('d', "display", "display fitted plots", true);

    if(!BrAppOptionManager::Instance()->ProcessCommandLine()) { return 1; }
    if(BrAppOptionManager::Instance()->ShowVersion()) { return 0; }
    if(BrAppOptionManager::Instance()->ShowHelp()) { return 0; }

    TApplication app("app", &argc, argv);



    // Face calibration parameter read-in
    double front_vs_back_slope {}, front_vs_back_offset {};
    double front_vs_front_slope {}, front_vs_front_offset {}, back_vs_back_slope {}, back_vs_back_offset {};
    read_face_parameters(Form("%sfront_vs_back.dat", CALIB_FILE_DIR), front_vs_back_slope, front_vs_back_offset);
    read_face_parameters(Form("%sfront_vs_front.dat", CALIB_FILE_DIR), front_vs_front_slope, front_vs_front_offset);
    read_face_parameters(Form("%sback_vs_back.dat", CALIB_FILE_DIR), back_vs_back_slope, back_vs_back_offset);
    
    // Reduced file setup and tree read-in
    std::cout << std::endl;
    std::cout << "check_gain_match : reduced file read-in" << std::endl;
    auto in_file = std::make_unique<TFile>(Form(Form("%s%s", REDUCED_FILE_DIR, FILE_NAME_FORMAT), EXPT_DATE, opt_run_number->GetValue()));
    gROOT->cd();
    auto tree = dynamic_cast<TTree *>(in_file->FindObjectAny("T"));
    T063123Event *evin {};
    tree->SetBranchAddress("event", &evin);
    auto num_entries {std::min(static_cast<int>(tree->GetEntries()), opt_events->GetValue())};

    // Histograms/graphs needed for calibration
    auto *g_EFvsY_Uncorrected = new TGraph(); g_EFvsY_Uncorrected->SetTitle(";Y;Front (chan)");
    auto *g_EFvsY_Corrected = new TGraph(); g_EFvsY_Corrected->SetTitle(";Y;Front (chan)");
    auto *g_EBvsX_Uncorrected = new TGraph(); g_EBvsX_Uncorrected->SetTitle(";X;Back (chan)");
    auto *g_EBvsX_Corrected = new TGraph(); g_EBvsX_Corrected->SetTitle(";X;Back (chan)");

    // Event loop
    std::cout << "check_gain_match : entering event loop" << std::endl;
    for(auto i = 0; i < num_entries; ++i) {
        if(i % 10000 == 0) printf("Event %d/%d\n", i, num_entries);
        tree->GetEntry(i);

        double f1 {}, f2 {}, b1 {}, b2 {};
        integrator_method(evin, f1, f2, b1, b2);

        double front_sum {f1 + f2}, back_sum {b1 + b2};
        if(front_sum > 0 && front_sum < 1E10 && back_sum > 0 && back_sum < 1E10) {
            if(front_equal_back(front_sum, back_sum, front_vs_back_slope, front_vs_back_offset)) {
                g_EFvsY_Uncorrected->SetPoint(g_EFvsY_Uncorrected->GetN(), (f2 - f1)/(f1 + f2), front_sum);
                g_EBvsX_Uncorrected->SetPoint(g_EBvsX_Uncorrected->GetN(), (b2 - b1)/(b1 + b2), back_sum);

                auto front_corrected {gain_match(f1, f2, front_vs_front_slope, front_vs_front_offset)};
                auto back_corrected {gain_match(b1, b2, back_vs_back_slope, back_vs_back_offset)};

                g_EFvsY_Corrected->SetPoint(g_EFvsY_Corrected->GetN(), (f2 - f1)/(f1 + f2), front_corrected);
                g_EBvsX_Corrected->SetPoint(g_EBvsX_Corrected->GetN(), (b2 - b1)/(b1 + b2), back_corrected);
            }
        }
    }

    auto *can = (TCanvas *)gROOT->FindObjectAny("can");
    if(!can) can = new TCanvas("can", "", 50, 50, 700, 700);
    can->Divide(2, 2);
    can->cd(1);
    g_EFvsY_Uncorrected->Draw("ap");
    can->cd(2);
    g_EFvsY_Corrected->Draw("ap");
    can->cd(3);
    g_EBvsX_Uncorrected->Draw("ap");
    can->cd(4);
    g_EBvsX_Corrected->Draw("ap");

    in_file->Close();

    std::cout << "check_gain_match : finished" << std::endl;
    std::cout << std::endl;

    if(opt_display->GetValue()) {
        TRootCanvas *root_canvas = (TRootCanvas *)can->GetCanvasImp();
        root_canvas->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

        app.Run();
    }

    return 0;
}