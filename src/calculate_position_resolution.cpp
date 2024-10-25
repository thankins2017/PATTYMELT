// Calculates the position resolution for the DADL by creating a position map of the highest energy
// alpha particles and fitting with a sigmoid function to determine the tailing factor off the side

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
    BrAppOptionManager::Instance()->SetVersion(1, 0, "calculate_position_resolution");
    BrAppOptionManager::Instance()->SetHelp("Calculates target testing station DADL position resolution.");
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

    // Front vs. front calibration parameter read-in
    double front_vs_front_slope {}, front_vs_front_offset {};
    in_file_calib.open(Form("%sfront_vs_front.dat", CALIB_FILE_DIR), std::ios_base::in);
    std::getline(in_file_calib, in_line);
    buffer.clear(); buffer.str(in_line);
    buffer >> front_vs_front_slope >> front_vs_front_offset;
    in_file_calib.close();

    // Back vs. back calibration parameter read-in
    double back_vs_back_slope {1.0}, back_vs_back_offset {0.0};
    in_file_calib.open(Form("%sback_vs_back.dat", CALIB_FILE_DIR), std::ios_base::in);
    std::getline(in_file_calib, in_line);
    buffer.clear(); buffer.str(in_line);
    buffer >> back_vs_back_slope >> back_vs_back_offset;
    in_file_calib.close();

    // Stretching parameter read-in
    double top_limit {}, bottom_limit {}, right_limit {}, left_limit{};
    in_file_calib.open(Form("%sstretching_parameters.dat", CALIB_FILE_DIR), std::ios_base::in);
    std::getline(in_file_calib, in_line);
    buffer.clear(); buffer.str(in_line);
    buffer >> top_limit >> bottom_limit >> right_limit >> left_limit;
    in_file_calib.close();

    // Energy calibration read-in
    double energy_slope {}, energy_intercept {};
    in_file_calib.open(Form("%senergy_calibration.dat", CALIB_FILE_DIR), std::ios_base::in);
    std::getline(in_file_calib, in_line);
    buffer.clear(); buffer.str(in_line);
    buffer >> energy_slope >> energy_intercept;
    in_file_calib.close();
    
    // Reduced file setup and tree read-in
    std::cout << std::endl;
    std::cout << "stretching_parameters : reduced file read-in" << std::endl;
    auto in_file = std::make_unique<TFile>(Form(Form("%s%s", REDUCED_FILE_DIR, FILE_NAME_FORMAT), EXPT_DATE, opt_run_number->GetValue()));
    gROOT->cd();
    auto tree = dynamic_cast<TTree *>(in_file->FindObjectAny("T"));
    T063123Event *evin {};
    tree->SetBranchAddress("event", &evin);
    auto num_entries {std::min(static_cast<int>(tree->GetEntries()), opt_events->GetValue())};

    // CycSrim initialization, calculation of energy points in source spectrum
    std::cout << "perform_thickness_measurement : calculating adjusted source spectrum" << std::endl;
    auto *gold = new CycSrim(CycSrim::SrimMaterialAu, 0.026, CycSrim::kUnitsMicron);
    auto *mylar = new CycSrim(CycSrim::SrimMaterialMylar, 1.4, CycSrim::kUnitsMicron);
    auto *dead_layer = new CycSrim(CycSrim::SrimMaterialSi, 0.5, CycSrim::kUnitsMicron);

    for(auto i = 0; i < source_energies.size(); ++i) {
        source_energies.at(i) = gold->GetResidualEnergy(2, 4, source_energies.at(i)); // Loss through source gold layer
        source_energies.at(i) = mylar->GetResidualEnergy(2, 4, source_energies.at(i)); // Loss through source mylar cover
        source_energies.at(i) = dead_layer->GetResidualEnergy(2, 4, source_energies.at(i)); // Loss through detector dead layer
        std::cout << source_energies.at(i) << ", ";
    }
    std::cout << std::endl;

    // Histograms/graphs needed for calibration
    auto *h_X = new TH1F("h_X", "h_X;X;Yield", 480, -1.2, 1.2);
    auto *h_Y = new TH1F("h_Y", "h_Y;Y;Yield", 480, -1.2, 1.2);
    auto *h_Energy = new TH1F("h_Energy", "h_Energy;Energy (MeV);Yield", 2000, 0, 10);

    // Calculate resolution event loop
    std::cout << "calculate_position_resolution : entering event loop" << std::endl;
    for(auto i = 0; i < num_entries; ++i) {
        if(i % 10000 == 0) printf("Event %d/%d\n", i, num_entries);
        tree->GetEntry(i);

        double f1 {}, f2 {}, b1 {}, b2 {};
        integrator_method(evin, f1, f2, b1, b2);

        double front_sum {f1 + f2}, back_sum {b1 + b2};
        if(front_sum > 0 && front_sum < 1E10 && back_sum > 0 && back_sum < 1E10) {
            if(front_equal_back(front_sum, back_sum, front_vs_back_slope, front_vs_back_offset)) {
                auto front_corrected {gain_match(f1, f2, front_vs_front_slope, front_vs_front_offset)};
                auto back_corrected {gain_match(b1, b2, back_vs_back_slope, back_vs_back_offset)};
                auto energy {(front_corrected + back_corrected) / 2.0};

                auto calib_energy {energy_slope * energy + energy_intercept};

                h_Energy->Fill(calib_energy);

                if(calib_energy > 8.0) {
                    auto local_y {stretch_position(f1, f2, top_limit, bottom_limit)};
                    auto local_x {stretch_position(b1, b2, right_limit, left_limit)}; 

                    h_X->Fill(local_x);
                    h_Y->Fill(local_y);
                } 
            }
        }
    }

    // Create the sigmoid functions that will be used to extract the stretching parameters, then fit
    auto *f_xpos = new TF1("f_xpos", "[0]*(1/(1 + TMath::Exp((x - [1])/[3])) + 1/(1 + TMath::Exp((x - [2])/(-1.0*[3]))) - 1)", -1.2, 1.2);
    f_xpos->SetParameters(0.95 * h_X->GetMaximum(), 1.0, -1.0, 0.025);
    f_xpos->SetParLimits(1, 0.8, 1.2);
    f_xpos->SetParLimits(2, -1.2, -0.8);
    f_xpos->SetParLimits(3, 0.0, 0.05);
    h_X->Fit(f_xpos, "NQ");
    auto *f_ypos = new TF1("f_ypos", "[0]*(1/(1 + TMath::Exp((x - [1])/[3])) + 1/(1 + TMath::Exp((x - [2])/(-1.0*[3]))) - 1)", -1.2, 1.2);
    f_ypos->SetParameters(0.95 * h_Y->GetMaximum(), 1.0, -1.0, 0.025);
    f_ypos->SetParLimits(1, 0.8, 1.2);
    f_ypos->SetParLimits(2, -1.2, -0.8);
    f_ypos->SetParLimits(3, 0.0, 0.05);
    h_Y->Fit(f_ypos, "NQ");

    // sigma_sigmoid = 0.5300(2) * sigma_gaussian (Andy thesis)
    std::cout << std::endl;
    std::cout << "(x, y) sigmoid: (" << f_xpos->GetParameter(3) << ", " << f_ypos->GetParameter(3) << ")" << std::endl; 
    std::cout << "(x, y) FWHM: (" << f_xpos->GetParameter(3) / 0.53 * 2.355 << ", " << f_ypos->GetParameter(3) / 0.53 * 2.355 << ")" << std::endl;
    std::cout << "Avg.: " << (f_xpos->GetParameter(3) / 0.53 * 2.355 + f_ypos->GetParameter(3) / 0.53 * 2.355) / 2.0 << std::endl;
    std::cout << std::endl;



    // Draw results and finalize
    auto *can = (TCanvas *)gROOT->FindObjectAny("can");
    if(!can) can = new TCanvas("can", "", 50, 50, 350, 350);
    can->cd();
    h_Energy->SetStats(false);
    h_Energy->Draw("hist");

    auto *can2 = (TCanvas *)gROOT->FindObjectAny("can2");
    if(!can2) can2 = new TCanvas("can2", "", 400, 50, 700, 700);
    can2->Divide(2, 2);
    can2->cd(1);
    h_X->SetStats(false);
    h_X->Draw("hist");
    f_xpos->Draw("l same");
    can2->cd(2);
    h_Y->SetStats(false);
    h_Y->Draw("hist");
    f_ypos->Draw("l same");
 
    in_file->Close();

    std::cout << "calculate_position_resolution : finished" << std::endl;
    std::cout << std::endl;

    if(opt_display->GetValue()) {
        TRootCanvas *root_canvas = (TRootCanvas *)can->GetCanvasImp();
        root_canvas->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

        app.Run();
    }

    return 0;
}