// Using source information, performs energy loss calculations and determines energy calibration
// of the testing DADL.

// C++ includes
#include <fstream>
#include <string>
#include <sstream>

// ROOT includes
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TROOT.h"
#include "TSpectrum.h" // target_link_libraries(Spectrum)

// Other includes
#include "PathManager.h"
#include "AnalysisTools.h"

int main(int argc, char **argv) {
    // option manager
    BrAppOptionManager::Instance()->SetVersion(1, 0, "energy_calibration");
    BrAppOptionManager::Instance()->SetHelp("Energy calibrates the testing DADL.");
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
    double back_vs_back_slope {}, back_vs_back_offset {};
    in_file_calib.open(Form("%sback_vs_back.dat", CALIB_FILE_DIR), std::ios_base::in);
    std::getline(in_file_calib, in_line);
    buffer.clear(); buffer.str(in_line);
    buffer >> back_vs_back_slope >> back_vs_back_offset;
    in_file_calib.close();
    
    // Reduced file setup and tree read-in
    std::cout << std::endl;
    std::cout << "energy_calibration : reduced file read-in" << std::endl;
    auto in_file = std::make_unique<TFile>(Form(Form("%s%s", REDUCED_FILE_DIR, FILE_NAME_FORMAT), EXPT_DATE, opt_run_number->GetValue()));
    gROOT->cd();
    auto tree = dynamic_cast<TTree *>(in_file->FindObjectAny("T"));
    T063123Event *evin {};
    tree->SetBranchAddress("event", &evin);
    auto num_entries {std::min(static_cast<int>(tree->GetEntries()), opt_events->GetValue())};

    // CycSrim initialization, calculation of energy points in source spectrum
    std::cout << "energy_calibration : calculating adjusted source spectrum" << std::endl;
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
    auto *h_E_chan = new TH1F("h_E_chan", "", 2000, 0, 2000);

    // Energy calibration event loop
    std::cout << "energy_calibration : entering event loop" << std::endl;
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

                h_E_chan->Fill((front_corrected + back_corrected) / 2.0);
            }
        }
    }

    auto *can = (TCanvas *)gROOT->FindObjectAny("can");
    if(!can) can = new TCanvas("can", "", 50, 50, 500, 500);
    h_E_chan->SetStats(false);
    h_E_chan->Draw("hist");

    // Set up the TSpectrum to identify peaks.
    auto *spec = new TSpectrum(2 * source_energies.size());
    auto threshold {0.35}; // determines which peaks to cut off for consideration.
    auto peak_list_length = spec->Search(h_E_chan, 5, "", threshold);
    double *x_peaks = spec->GetPositionX();
    double par[2 * source_energies.size() + 1];
    par[0] = 5; // universal width of the peaks, in channels
    for(auto i = 0; i < peak_list_length; ++i) {
        par[2*i + 1] = h_E_chan->GetBinContent(h_E_chan->GetXaxis()->FindBin(x_peaks[i])); // Set the approx. norm.
        par[2*i + 2] = x_peaks[i]; // Set the mean

    }

    // Fit the peaks that were found.
    auto *f_spectrum = new TF1("f_spectrum", spectrum, 0, 2000, 2*source_energies.size() + 1);
    f_spectrum->SetNpx(2000);
    f_spectrum->SetParameters(par);
    for(auto i = 0; i < 2 * source_energies.size() + 1; ++i) {
        f_spectrum->SetParLimits(i, f_spectrum->GetParameter(i) * 0.75, f_spectrum->GetParameter(i) * 1.25);
    }
    h_E_chan->Fit(f_spectrum, "EBQ");
    auto *c1 = (TCanvas*)gROOT->FindObjectAny("c1");
    // if(c1) delete c1; // Because TSpectrum can't handle the canvas it creates...

    // Extract the peaks, then sort in descending order
    std::vector<double> peak_centroids {};
    for(auto i = 0; i < (sizeof(par)/sizeof(double) - 1)/2; ++i) {
        if(par[2 +2*i] > 5) peak_centroids.push_back(par[2 + 2*i]);
    }
    std::sort(peak_centroids.begin(), peak_centroids.end(), greater<double>());

    // Build the energy conversion graph, plot, and fit.
    auto *g_E_calib = new TGraph(); g_E_calib->SetMarkerStyle(20);
    for(auto i = 0; i < peak_centroids.size(); ++i) {
        g_E_calib->SetPoint(g_E_calib->GetN(), peak_centroids.at(i), source_energies.at(i));
    }

    auto *can2 = (TCanvas *)gROOT->FindObjectAny("can2");
    if(!can2) can2 = new TCanvas("can2", "", 50, 50, 500, 500);
    g_E_calib->Draw("ap");

    auto *f_E_calib = new TF1("f_E_calib", "[0]*x + [1]", 0, 2000);
    g_E_calib->Fit(f_E_calib, "BQ");
    f_E_calib->Draw("l same");

    std::cout << "energy_calibration : conversion result" << std::endl;
    std::cout << Form("E (MeV) = m * E (chan) + b; m = %.3f, b = %.3f", f_E_calib->GetParameter(0), f_E_calib->GetParameter(1)) << std::endl;

    std::fstream out_file(Form("%s%s", CALIB_FILE_DIR, "energy_calibration.dat"), std::ios_base::out);
    out_file << f_E_calib->GetParameter(0) << "\t" << f_E_calib->GetParameter(1) << std::endl;
    out_file.close();

    auto out_file_t = std::make_unique<TFile>(Form("%senergy_calib_%s%03d.root", OUTPUT_FILE_DIR, EXPT_DATE, opt_run_number->GetValue()), "RECREATE");
    h_E_chan->Write();
    out_file_t->Close();

    std::cout << "energy_calibration : finished" << std::endl;
    std::cout << std::endl;

    if(opt_display->GetValue()) {
        TRootCanvas *root_canvas = (TRootCanvas *)can->GetCanvasImp();
        root_canvas->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

        app.Run();
    }

    return 0;
}