// Determines the front vs. back slope and offset from a calibration run. Gain matching on the
// faces is needed before running this procedure.

// C++ includes
#include <fstream>

// ROOT includes
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TROOT.h"

// Other includes
#include "PathManager.h"
#include "AnalysisTools.h"

int main(int argc, char **argv) {
    // option manager
    BrAppOptionManager::Instance()->SetVersion(1, 0, "front_vs_back");
    BrAppOptionManager::Instance()->SetHelp("Calibrates target testing station DADL front vs. back.");
    BrAppOptionManager::Instance()->SetCommandLine(argc, argv);
    auto *opt_run_number = new BrAppIntOption('r', "run_number", "run number", CALIB_RUN);
    auto *opt_events     = new BrAppIntOption('e', "events", "events to analyze", std::numeric_limits<int>::max());
    auto *opt_display    = new BrAppBoolOption('d', "display", "display fitted plots", true);

    if(!BrAppOptionManager::Instance()->ProcessCommandLine()) { return 1; }
    if(BrAppOptionManager::Instance()->ShowVersion()) { return 0; }
    if(BrAppOptionManager::Instance()->ShowHelp()) { return 0; }

    TApplication app("app", &argc, argv);



    // Reduced file setup and tree read-in
    std::cout << std::endl;
    std::cout << "front_vs_back : reduced file read-in" << std::endl;
    auto in_file = std::make_unique<TFile>(Form(Form("%s%s", REDUCED_FILE_DIR, FILE_NAME_FORMAT), EXPT_DATE, opt_run_number->GetValue()));
    gROOT->cd();
    auto tree = dynamic_cast<TTree *>(in_file->FindObjectAny("T"));
    T063123Event *evin {};
    tree->SetBranchAddress("event", &evin);
    auto num_entries {std::min(static_cast<int>(tree->GetEntries()), opt_events->GetValue())};

    // Face calibration parameter read-in
    double front_vs_front_slope {}, front_vs_front_offset {}, back_vs_back_slope {}, back_vs_back_offset {};
    read_face_parameters(Form("%sfront_vs_front.dat", CALIB_FILE_DIR), front_vs_front_slope, front_vs_front_offset);
    read_face_parameters(Form("%sback_vs_back.dat", CALIB_FILE_DIR), back_vs_back_slope, back_vs_back_offset);

    // Histograms/graphs needed for calibration
    auto *g_FrontVsBack = new TGraph(); g_FrontVsBack->SetTitle("FvB;Back Sum (chan);Front Sum (chan)");

    // Front versus back event loop
    std::cout << "front_vs_back : entering front vs. back event loop" << std::endl;
    for(auto i = 0; i < num_entries; ++i) {
        if(i % 10000 == 0) printf("Event %d/%d\n", i, num_entries);
        tree->GetEntry(i);

        double f1 {}, f2 {}, b1 {}, b2 {};
        integrator_method(evin, f1, f2, b1, b2);

        double front_sum {gain_match(f1, f2, front_vs_front_slope, front_vs_front_offset)};
        double back_sum  {gain_match(b1, b2, back_vs_back_slope, back_vs_back_offset)};
        // front_equal_back function not used since slope and intercept for this call is extracted
        // from this procedure.
        if(front_sum > 0 && front_sum < 1E10 && back_sum > 0 && back_sum < 1E10) {
            if(front_sum < 1.05 * back_sum && front_sum > 0.95 * back_sum) {
                g_FrontVsBack->SetPoint(g_FrontVsBack->GetN(), back_sum, front_sum);
            }
        }
    }

    // Create front vs back fit, extract parameters, print to terminal
    auto f_FrontVsBack = new TF1("f_FrontVsBack", "[0]*x + [1]", 0, 4000);
    f_FrontVsBack->SetParameters(1, 0);
    g_FrontVsBack->Fit(f_FrontVsBack, "Q");

    auto *can = (TCanvas *)gROOT->FindObjectAny("can");
    if(!can) can = new TCanvas("can", "", 50, 50, 500, 500);
    g_FrontVsBack->Draw("ap");
    f_FrontVsBack->Draw("l same");

    auto front_vs_back_slope {f_FrontVsBack->GetParameter(0)}, front_vs_back_offset {f_FrontVsBack->GetParameter(1)};

    std::cout << "front_vs_back : finished with front vs. back event loop" << std::endl;
    std::cout << "front_vs_back : FvB (slope, offset) - " << Form("(%.3f, %.3f)", front_vs_back_slope, front_vs_back_offset) << std::endl;

    std::fstream out_file(Form("%s%s", CALIB_FILE_DIR, "front_vs_back.dat"), std::ios_base::out);
    out_file << front_vs_back_slope << "\t" << front_vs_back_offset << std::endl;
    out_file.close();

    in_file->Close();

    std::cout << "front_vs_back : finished" << std::endl;
    std::cout << std::endl;

    if(opt_display->GetValue()) {
        TRootCanvas *root_canvas = (TRootCanvas *)can->GetCanvasImp();
        root_canvas->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

        app.Run();
    }

    return 0;
}