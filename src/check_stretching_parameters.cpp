// Verifies the stretching parameters extracted with stretching_parameters. This is essentially a 
// carbon copy of stretching_parameters, but the parameters are already known. Originally made for
// use with mask data to certify method.

// C++ includes
#include <fstream>
#include <string>
#include <sstream>

// ROOT includes
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TROOT.h"
#include "TLine.h"

// Other includes
#include "PathManager.h"
#include "AnalysisTools.h"

int main(int argc, char **argv) {
    // option manager
    BrAppOptionManager::Instance()->SetVersion(1, 0, "check_stretching_parameters");
    BrAppOptionManager::Instance()->SetHelp("Checks stretching parameters.");
    BrAppOptionManager::Instance()->SetCommandLine(argc, argv);
    auto *opt_run_number = new BrAppIntOption('r', "run_number", "run number", CALIB_RUN);
    auto *opt_events     = new BrAppIntOption('e', "events", "events to analyze", std::numeric_limits<int>::max());
    auto *opt_display    = new BrAppBoolOption('d', "display", "display fitted plots", true);

    if(!BrAppOptionManager::Instance()->ProcessCommandLine()) { return 1; }
    if(BrAppOptionManager::Instance()->ShowVersion()) { return 0; }
    if(BrAppOptionManager::Instance()->ShowHelp()) { return 0; }

    TApplication app("app", &argc, argv);



    // Front vs. back calibration parameter read-in
    double front_vs_back_slope {}, front_vs_back_offset {};
    std::fstream in_file_calib(Form("%sfront_vs_back.dat", CALIB_FILE_DIR), std::ios_base::in);
    std::string in_line {}; std::getline(in_file_calib, in_line);
    std::istringstream buffer(in_line); buffer >> front_vs_back_slope >> front_vs_back_offset;
    in_file_calib.close();

    // Stretching parameter read-in
    double right_limit {}, left_limit {}, top_limit {}, bottom_limit {};
    in_file_calib.open(Form("%sstretching_parameters.dat", CALIB_FILE_DIR), std::ios_base::in);
    std::getline(in_file_calib, in_line);
    buffer.clear(); buffer.str(in_line);
    buffer >> top_limit >> bottom_limit >> right_limit >> left_limit;
    in_file_calib.close();
    
    // Reduced file setup and tree read-in
    std::cout << std::endl;
    std::cout << "check_stretching_parameters : reduced file read-in" << std::endl;
    auto in_file = std::make_unique<TFile>(Form(Form("%s%s", REDUCED_FILE_DIR, FILE_NAME_FORMAT), EXPT_DATE, opt_run_number->GetValue()));
    gROOT->cd();
    auto tree = dynamic_cast<TTree *>(in_file->FindObjectAny("T"));
    T063123Event *evin {};
    tree->SetBranchAddress("event", &evin);
    auto num_entries {std::min(static_cast<int>(tree->GetEntries()), opt_events->GetValue())};

    // Histograms/graphs needed for calibration
    auto *g_XY_Uncorrected = new TGraph(); g_XY_Uncorrected->SetTitle("Position;X;Y");
    auto *g_XY_Corrected = new TGraph(); g_XY_Corrected->SetTitle("Position;X;Y");

    // Check stretching parameter event loop
    std::cout << "check_stretching_parameters : entering event loop" << std::endl;
    for(auto i = 0; i < num_entries; ++i) {
        if(i % 10000 == 0) printf("Event %d/%d\n", i, num_entries);
        tree->GetEntry(i);

        double f1 {}, f2 {}, b1 {}, b2 {};
        integrator_method(evin, f1, f2, b1, b2);

        double front_sum {f1 + f2}, back_sum {b1 + b2};
        if(front_sum > 0 && front_sum < 1E10 && back_sum > 0 && back_sum < 1E10) {
            if(front_equal_back(front_sum, back_sum, front_vs_back_slope, front_vs_back_offset)) {
                double x_pos {(b2 - b1) / (b1 + b2)}, y_pos {(f2 - f1) / (f1 + f2)};
                g_XY_Uncorrected->SetPoint(g_XY_Uncorrected->GetN(), x_pos, y_pos);
                g_XY_Corrected->SetPoint(g_XY_Corrected->GetN(),
                                         stretch_position(b1, b2, right_limit, left_limit),
                                         stretch_position(f1, f2, top_limit, bottom_limit));
            }
        }
    }



    // Draw result and finalize
    auto *can = (TCanvas *)gROOT->FindObjectAny("can");
    if(!can) can = new TCanvas("can", "", 50, 50, 1000, 500);
    can->Divide(2, 1);
    can->cd(1);

    auto *h_Dummy = new TH1C("h_Dummy", "", 1, -1.5, 1.5);
    h_Dummy->SetStats(false); h_Dummy->SetLineColor(kWhite);
    h_Dummy->GetYaxis()->SetRangeUser(-1.5, 1.5);
    h_Dummy->Draw("hist");

    std::vector<TLine *> lines {};
    lines.push_back(new TLine(-1.25, -1, 1.25, -1));
    lines.push_back(new TLine(-1.25, 1, 1.25, 1));
    lines.push_back(new TLine(-1, -1.25, -1, 1.25));
    lines.push_back(new TLine(1, -1.25, 1, 1.25));
    for(auto *line : lines) line->Draw("l same");

    g_XY_Uncorrected->Draw("p same");
    can->cd(2);
    h_Dummy->Draw("hist");
    g_XY_Corrected->Draw("p same");
    for(auto *line : lines) line->Draw("l same");

    in_file->Close();

    std::cout << "check_stretching_parameters : finished" << std::endl;
    std::cout << std::endl;

    if(opt_display->GetValue()) {
        TRootCanvas *root_canvas = (TRootCanvas *)can->GetCanvasImp();
        root_canvas->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

        app.Run();
    }

    return 0;
}