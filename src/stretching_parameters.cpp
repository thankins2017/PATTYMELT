// Determines the DADL stretching parameters from blank (source) data.

// C++ includes
#include <fstream>

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
    BrAppOptionManager::Instance()->SetVersion(1, 0, "stretching_parameters");
    BrAppOptionManager::Instance()->SetHelp("Extracts target testing station DADL stretching parameters.");
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
    read_face_parameters(Form("%sfront_vs_back.dat", CALIB_FILE_DIR), front_vs_back_slope, front_vs_back_offset);
    
    // Reduced file setup and tree read-in
    std::cout << std::endl;
    std::cout << "stretching_parameters : reduced file read-in" << std::endl;
    auto in_file = std::make_unique<TFile>(Form(Form("%s%s", REDUCED_FILE_DIR, FILE_NAME_FORMAT), EXPT_DATE, opt_run_number->GetValue()));
    gROOT->cd();
    auto tree = dynamic_cast<TTree *>(in_file->FindObjectAny("T"));
    T063123Event *evin {};
    tree->SetBranchAddress("event", &evin);
    auto num_entries {std::min(static_cast<int>(tree->GetEntries()), opt_events->GetValue())};

    // Histograms/graphs needed for calibration
    auto *g_XY_Uncorrected = new TGraph(); g_XY_Uncorrected->SetTitle("Position;X;Y");
    auto *g_XY_Corrected = new TGraph(); g_XY_Corrected->SetTitle("Position;X;Y");
    auto *h_X = new TH1F("h_X", "h_X;X;Yield", 400, -1, 1);
    auto *h_Y = new TH1F("h_Y", "h_Y;Y;Yield", 400, -1, 1);

    // Stretching parameter event loop
    std::cout << "stretching_parameters : entering event loop" << std::endl;
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
                h_X->Fill(x_pos);
                h_Y->Fill(y_pos);
            }
        }
    }

    // Create the sigmoid functions that will be used to extract the stretching parameters, then fit
    auto *f_xpos = new TF1("f_xpos", "[0]*(1/(1 + TMath::Exp((x - [1])/[3])) + 1/(1 + TMath::Exp((x - [2])/(-1.0*[3]))) - 1)", -0.9, 0.9);
    f_xpos->SetParameters(0.95 * h_X->GetMaximum(), 0.4, -0.4, 0.025);
    f_xpos->SetParLimits(1, 0.2, 0.8);
    f_xpos->SetParLimits(2, -0.8, -0.2);
    f_xpos->SetParLimits(3, 0.0, 0.05);
    h_X->Fit(f_xpos, "NQ");
    auto *f_ypos = new TF1("f_ypos", "[0]*(1/(1 + TMath::Exp((x - [1])/[3])) + 1/(1 + TMath::Exp((x - [2])/(-1.0*[3]))) - 1)", -0.9, 0.9);
    f_ypos->SetParameters(0.95 * h_Y->GetMaximum(), 0.4, -0.4, 0.025);
    f_ypos->SetParLimits(1, 0.2, 0.8);
    f_ypos->SetParLimits(2, -0.8, -0.2);
    f_ypos->SetParLimits(3, 0.0, 0.05);
    h_Y->Fit(f_ypos, "NQ");

    std::cout << "stretching_parameters : fits complete" << std::endl;
    std::cout << "Top, Bottom, Right, Left Limits:" << std::endl;
    std::cout << f_ypos->GetParameter(1) << "\t" << f_ypos->GetParameter(2) << "\t" 
              << f_xpos->GetParameter(1) << "\t" << f_xpos->GetParameter(2) << std::endl;

    // Stretching parameter event loop, this time to scale
    std::cout << "stretching_parameters : entering event loop 2" << std::endl;
    for(auto i = 0; i < num_entries; ++i) {
        if(i % 10000 == 0) printf("Event %d/%d\n", i, num_entries);
        tree->GetEntry(i);

        double f1 {}, f2 {}, b1 {}, b2 {};
        integrator_method(evin, f1, f2, b1, b2);

        double front_sum {f1 + f2}, back_sum {b1 + b2};
        if(front_sum > 0 && front_sum < 1E10 && back_sum > 0 && back_sum < 1E10) {
            if(front_equal_back(front_sum, back_sum, front_vs_back_slope, front_vs_back_offset)) {
                double x_pos {(b2 - b1) / (b1 + b2)}, y_pos {(f2 - f1) / (f1 + f2)};
                g_XY_Corrected->SetPoint(g_XY_Corrected->GetN(),
                                         stretch_position(b1, b2, f_xpos->GetParameter(1), f_xpos->GetParameter(2)),
                                         stretch_position(f1, f2, f_ypos->GetParameter(1), f_ypos->GetParameter(2)));
            }
        }
    }



    // Draw results and finalize
    auto *can = (TCanvas *)gROOT->FindObjectAny("can");
    if(!can) can = new TCanvas("can", "", 50, 50, 700, 700);
    can->Divide(2, 2);
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
    can->cd(3);
    h_X->SetStats(false);
    h_X->Draw("hist");
    f_xpos->Draw("l same");
    can->cd(4);
    h_Y->SetStats(false);
    h_Y->Draw("hist");
    f_ypos->Draw("l same");

    std::fstream out_file(Form("%s%s", CALIB_FILE_DIR, "stretching_parameters.dat"), std::ios_base::out);
    out_file << f_ypos->GetParameter(1) << "\t" << f_ypos->GetParameter(2) << "\t" 
             << f_xpos->GetParameter(1) << "\t" << f_xpos->GetParameter(2) << std::endl;
    out_file.close();
 
    in_file->Close();

    std::cout << "stretching_parameters : finished" << std::endl;
    std::cout << std::endl;

    if(opt_display->GetValue()) {
        TRootCanvas *root_canvas = (TRootCanvas *)can->GetCanvasImp();
        root_canvas->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

        app.Run();
    }

    return 0;
}