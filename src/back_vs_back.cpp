// Determines the back vs. back relationship from HDBSCAN clustering.

// C++ includes
#include <fstream>

// ROOT includes
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TROOT.h"
#include "TMarker.h"
#include "TPrincipal.h"
#include "TF1.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TMatrixDSymEigen.h"

// Other includes
#include "PathManager.h"
#include "AnalysisTools.h"
#include "hdbscan.h"

int get_cluster_color(int index) {
    if(index < 0) return kBlack;

    const int number_colors {12};
    int colors[number_colors] {static_cast<int>(kRed), static_cast<int>(kRed + 2), static_cast<int>(kMagenta), static_cast<int>(kMagenta + 2), 
                               static_cast<int>(kBlue), static_cast<int>(kBlue + 2), static_cast<int>(kCyan), static_cast<int>(kCyan + 2), 
                               static_cast<int>(kGreen), static_cast<int>(kGreen + 2), static_cast<int>(kYellow), static_cast<int>(kYellow + 2)};
    return colors[index % number_colors];
}

int main(int argc, char **argv) {
    // option manager
    BrAppOptionManager::Instance()->SetVersion(1, 0, "back_vs_back");
    BrAppOptionManager::Instance()->SetHelp("Calibrates target testing station DADL back vs. back.");
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
    std::cout << "back_vs_back : reduced file read-in" << std::endl;
    auto in_file = std::make_unique<TFile>(Form(Form("%s%s", REDUCED_FILE_DIR, FILE_NAME_FORMAT), EXPT_DATE, opt_run_number->GetValue()));
    gROOT->cd();
    auto tree = dynamic_cast<TTree *>(in_file->FindObjectAny("T"));
    T063123Event *evin {};
    tree->SetBranchAddress("event", &evin);
    auto num_entries {std::min(static_cast<int>(tree->GetEntries()), opt_events->GetValue())};

    // Histograms/graphs needed for calibration
    auto *g_BackSumDif = new TGraph(); g_BackSumDif->SetTitle("SumDif;Difference (chan);Sum (chan)");

    // Vector of Points to be passed to HDBSCAN clusterer
    std::vector<Point *> points {};

    // Back versus back event loop
    std::cout << "back_vs_back : entering event loop" << std::endl;
    int number_clustered {};
    for(auto i = 0; i < num_entries; ++i) {
        if(i % 10000 == 0) printf("Event %d/%d\n", i, num_entries);
        tree->GetEntry(i);

        double f1 {}, f2 {}, b1 {}, b2 {};
        integrator_method(evin, f1, f2, b1, b2);

        double front_sum {f1 + f2}, back_sum {b1 + b2};
        if(front_sum > 0 && front_sum < 1E10 && back_sum > 0 && back_sum < 1E10) {
            if(front_equal_back(front_sum, back_sum, front_vs_back_slope, front_vs_back_offset)) {
                g_BackSumDif->SetPoint(g_BackSumDif->GetN(), b2 - b1, back_sum);
                points.push_back(new Point({b2 - b1, back_sum}));
                ++number_clustered;
            }
        }
    }

    std::cout << "back_vs_back : number of points passed to clusterer - " << number_clustered << std::endl;

    // Display the results
    auto *can = (TCanvas *)gROOT->FindObjectAny("can");
    if(!can) can = new TCanvas("can", "", 750, 750);
    g_BackSumDif->Draw("ap");

    std::cout << "back_vs_back : initializing HDBSCAN clusterer" << std::endl;
    auto *clusterer {new HDBSCAN(points)};
    clusterer->set_minimum_cluster_size(1000);
    clusterer->set_alpha(0.8);
    clusterer->fit();

    // Determine number of clusters present, create associated number of TPrincipals, then fill with data. At the same
    // time, draw the clusters with their corresponding colors. Keep track of the mean x and y values for each cluster
    // by summing and dividing accordingly.
    int number_of_clusters {};
    for(auto *p : points) {
        if(p->get_cluster_id() > number_of_clusters) number_of_clusters = p->get_cluster_id();
    }
    ++number_of_clusters; // since get_cluster_id() counts from zero

    std::cout << "back_vs_back : initializing PCA objects" << std::endl;
    std::vector<TPrincipal *> pcas {};
    std::vector<int> num_points {}; // number of points for each cluster id
    std::vector<double> mean_x {}, mean_y {}; // average values for each cluster id
    for(int i = 0; i < number_of_clusters; ++i) {
        pcas.push_back(new TPrincipal(2, "D"));
        num_points.push_back(0);
        mean_x.push_back(0);
        mean_y.push_back(0);
    }

    double variables[2];
    for(auto *p : points) {
        TMarker *m {new TMarker(p->get_x(), p->get_y(), 20)};
        m->SetMarkerSize(0.1);
        m->SetMarkerColor(get_cluster_color(p->get_cluster_id()));
        m->Draw("p same");

        auto id {p->get_cluster_id()};
        if(id >= 0) {
            variables[0] = p->get_x();
            variables[1] = p->get_y();
            pcas.at(p->get_cluster_id())->AddRow(variables);

            ++num_points.at(p->get_cluster_id());
            mean_x.at(p->get_cluster_id()) += p->get_x();
            mean_y.at(p->get_cluster_id()) += p->get_y();
        }
    }

    for(int i = 0; i < static_cast<int>(mean_x.size()); ++i) {
        mean_x.at(i) /= static_cast<double>(num_points.at(i));
        mean_y.at(i) /= static_cast<double>(num_points.at(i));
    }

    // Perform the PCA calculations.
    std::cout << "back_vs_back : performing PCA fits" << std::endl;
    std::vector<TF1 *> fits {};
    for(int i = 0; i < static_cast<int>(pcas.size()); ++i) {
        pcas.at(i)->MakePrincipals();

        auto basis_matrix {*pcas.at(i)->GetEigenVectors()};
        auto eigenvalues {*pcas.at(i)->GetEigenValues()};
        auto lambda_1 {eigenvalues[0]}, lambda_2 {eigenvalues[1]};

        auto slope {(basis_matrix[1][0] * lambda_1 + basis_matrix[1][1] * lambda_2) /
                    (basis_matrix[0][0] * lambda_1 + basis_matrix[0][1] * lambda_2)};

        auto y_intercept {-1.0 * mean_x.at(i) * slope + mean_y.at(i)};
        auto x_intercept {-1.0 * mean_y.at(i) / slope + mean_x.at(i)};

        fits.push_back(new TF1(Form("line_%d", i), "[0]*x + [1]", -2000, 2000));
        fits.back()->SetParameters(slope, y_intercept);
        fits.back()->Draw("l same");
    }
    
    // Create final graph to compare intercepts. In SumDif space, evaluate the fits at
    // channels beyond where the data lie.
    auto *g_final = new TGraph(); g_final->SetTitle(";Low Chan. Eval.;High Chan. Eval.");
    g_final->SetMarkerStyle(20);
    for(auto *f : fits) {
        g_final->SetPoint(g_final->GetN(), f->GetParameter(0) * (-1500.0) + f->GetParameter(1),
                                           f->GetParameter(0) * (1500.0) + f->GetParameter(1));
    }
    auto *f_final = new TF1("f_final", "[0]*x + [1]", 0, 4000);
    g_final->Fit(f_final, "Q");
    std::cout << "Final fit: " << f_final->GetParameter(0) << " " << f_final->GetParameter(1) << std::endl;

    std::fstream out_file(Form("%s%s", CALIB_FILE_DIR, "back_vs_back.dat"), std::ios_base::out);
    out_file << f_final->GetParameter(0) << "\t" << f_final->GetParameter(1) << std::endl;
    out_file.close();

    in_file->Close();

    std::cout << "back_vs_back : finished" << std::endl;
    std::cout << std::endl;

    if(opt_display->GetValue()) {
        TRootCanvas *root_canvas = (TRootCanvas *)can->GetCanvasImp();
        root_canvas->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

        app.Run();
    }

    return 0;
}