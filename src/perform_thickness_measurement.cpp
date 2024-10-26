// After fully calibrating the DADL, perform a thickness measurement on a target.

// C++ includes
#include <fstream>
#include <string>
#include <sstream>
#include <numeric> // std::accumulate
#include <cmath>

// ROOT includes
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TROOT.h"
#include "TSpectrum.h" // target_link_libraries(Spectrum)

// Other includes
#include "PathManager.h"
#include "AnalysisTools.h"
#include "CycSrimHandler.h"

int main(int argc, char **argv) {
    //____________________________________________________________________________________________________
    // Configuration settings that generally do not change run to run.
    auto is_damaged_target {true}; // Adjusts parameters to expect significant non-uniformity //! In development
    auto verbosity {0}; // Adjust level of information printed as output.
	// 0 - no extra information; 1 - minimal information; 2 - more information; 3 - debug
    auto is_228Th {static_cast<bool>(source_energies.front() == 8.785)};
    //____________________________________________________________________________________________________

    // option manager
    BrAppOptionManager::Instance()->SetVersion(1, 0, "perform_thickness_measurement");
    BrAppOptionManager::Instance()->SetHelp("Measure target thicknesses.");
    BrAppOptionManager::Instance()->SetCommandLine(argc, argv);
    auto *opt_run_number = new BrAppIntOption('r', "run_number", "run number", CALIB_RUN);
    auto *opt_events     = new BrAppIntOption('e', "events", "events to analyze", std::numeric_limits<int>::max());
    auto *opt_material   = new BrAppStringOption('m', "material", "target material", "CycSrim::SrimMaterialAu");
    // Thickness map coarseness option. Physical dimensions of DADL are 2cm x 2cm, intrinsic alpha
    // position resolution for ~8 MeV deposited is ~0.5 mm. Doesn't make sense to plot more fine
    // than this. A warning is thrown if the fineness is dropped below this threshold.
    auto *opt_coarseness = new BrAppFloatOption('c', "coarseness", "map coarseness, in mm", 1.0); // In mm

    if(!BrAppOptionManager::Instance()->ProcessCommandLine()) { return 1; }
    if(BrAppOptionManager::Instance()->ShowVersion()) { return 0; }
    if(BrAppOptionManager::Instance()->ShowHelp()) { return 0; }



    // Target initialization; disregard the thickness given in the constructor, as it does not matter.
    auto *TARGET = new CycSrim(evaluate_cycsrim_material(opt_material->GetValue()), 1, CycSrim::kUnitsMgCm2);

    // Coarseness warning
    if(opt_coarseness->GetValue() < 0.5) {
        std::cout << "perform_thickness_measurement : coarseness option (" << opt_coarseness->GetValue()
                  << " mm) less than optimal for given conditions." << std::endl;
    }

    // Division warning
    // Histogram indexing relies on the DADL being divided by nice numbers (e.g., 0.1, 0.2, 0.25, etc.).
    if(std::abs(std::remainder(20.0, opt_coarseness->GetValue())) > 1e-6) { //std::fmod unreliable
        std::cout << "perform_thickness_measurement : coarseness option (" << opt_coarseness->GetValue()
                  << " mm) isn't evenly divisible into physical dimension (20 mm) (remainder: " 
                  << std::abs(std::remainder(20.0, opt_coarseness->GetValue())) << ")." << std::endl;
        std::abort();
    }

    // Face calibration parameter read-in
    double front_vs_back_slope {}, front_vs_back_offset {};
    double front_vs_front_slope {}, front_vs_front_offset {}, back_vs_back_slope {}, back_vs_back_offset {};
    read_face_parameters(Form("%sfront_vs_back.dat", CALIB_FILE_DIR), front_vs_back_slope, front_vs_back_offset);
    read_face_parameters(Form("%sfront_vs_front.dat", CALIB_FILE_DIR), front_vs_front_slope, front_vs_front_offset);
    read_face_parameters(Form("%sback_vs_back.dat", CALIB_FILE_DIR), back_vs_back_slope, back_vs_back_offset);

    // Stretching parameter and energy calibration read-in
    double top_limit {}, bottom_limit {}, right_limit {}, left_limit{};
    read_stretching_parameters(Form("%sstretching_parameters.dat", CALIB_FILE_DIR), top_limit, bottom_limit, right_limit, left_limit);
    double energy_slope {}, energy_intercept {};
    read_energy_parameters(Form("%senergy_calibration.dat", CALIB_FILE_DIR), energy_slope, energy_intercept);
    
    // Reduced file setup and tree read-in
    std::cout << std::endl;
    std::cout << "perform_thickness_measurement : reduced file read-in" << std::endl;
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
        // We are not interested in the energies that the detector observes; we are interested in the energies immediately after
        // the target. We do not calculate dead layer losses in this instance.
        // source_energies.at(i) = dead_layer->GetResidualEnergy(2, 4, source_energies.at(i)); // Loss through detector dead layer
        std::cout << source_energies.at(i) << ", ";
    }
    std::cout << std::endl;

    // Histograms/graphs needed for calibration
    auto out_file = std::make_unique<TFile>(Form("%sthickness_%s%03d.root", OUTPUT_FILE_DIR, EXPT_DATE, opt_run_number->GetValue()), "RECREATE");
    auto coarse {opt_coarseness->GetValue()};
    auto num_divisions {static_cast<int>(std::round(20.0 / coarse))};
    TH2F *h_ThicknessMap {}; // This will eventually be a clone of the last TH2F in the thickness_map histogram vector.
    TH2F *h_UpperError = new TH2F("h_UpperError", "", num_divisions, -1, 1, num_divisions, -1, 1);
    TH2F *h_LowerError = new TH2F("h_LowerError", "", num_divisions, -1, 1, num_divisions, -1, 1);

    std::vector<TH2F *> thickness_maps {};  // Thickness maps that will be added during evaluation
    auto *h_HitMap_Uncorrected = new TH2F("h_HitMap_Uncorrected", "", num_divisions, -1, 1, num_divisions, -1, 1);
    auto *h_HitMap = new TH2F("h_HitMap", "", num_divisions, -1, 1, num_divisions, -1, 1);
    auto *h_Energy = new TH1F("h_Energy", "", 2000, 0, 2000);
    std::vector<TH1F *> histograms {};
    for(auto y = -10 + 0.5 * coarse; y < 10; y += coarse) {
        for(auto x = -10 + 0.5 * coarse; x < 10; x += coarse) {
            histograms.push_back(new TH1F(Form("h_yx_%.2f_%.2f", y, x), "", 2000, 0, 2000));
        }
    }

    // Perform thickness measurement event loop
    std::cout << "perform_thickness_measurement : entering event loop" << std::endl;
    for(auto i = 0; i < num_entries; ++i) {
        if(i % 10000 == 0) printf("Event %d/%d\n", i, num_entries);
        tree->GetEntry(i);

        double f1 {}, f2 {}, b1 {}, b2 {};
        integrator_method(evin, f1, f2, b1, b2);

        h_HitMap_Uncorrected->Fill((b2 - b1)/ (b2 + b1), (f2 - f1) / (f2 + f1));

        double front_sum {f1 + f2}, back_sum {b1 + b2};
        if(front_sum > 0 && front_sum < 1E10 && back_sum > 0 && back_sum < 1E10) {
            if(front_equal_back(front_sum, back_sum, front_vs_back_slope, front_vs_back_offset)) {
                auto local_y {stretch_position(f1, f2, top_limit, bottom_limit)};
                auto local_x {stretch_position(b1, b2, right_limit, left_limit)}; 

                // Truncating information off the face of the DADL.
                if(local_y > -1.0 && local_y < 1.0 && local_x > -1.0 && local_x < 1.0) {
                    auto front_corrected {gain_match(f1, f2, front_vs_front_slope, front_vs_front_offset)};
                    auto back_corrected {gain_match(b1, b2, back_vs_back_slope, back_vs_back_offset)};

                    h_HitMap->Fill(local_x, local_y);

                    h_Energy->Fill((front_corrected + back_corrected) / 2.0);

                    // To help with indexing, the assignment is calculated by in-line shifting the local DADL 
                    // coordinates so that (-1, -1) becomes (0, 0). This is then converted to mm, and histogram 
                    // filling can be performed with modulo techiniques.

                    auto col_coord {static_cast<int>(((local_x + 1.0) * 10.0) / coarse)};
                    auto row_coord {static_cast<int>(((local_y + 1.0) * 10.0) / coarse)};

                    histograms.at(num_divisions * row_coord + col_coord)->Fill((front_corrected + back_corrected) / 2.0);
                }                
            } // front equals back
        } // raw sums are reasonable
    } // event loop

    std::cout << "perform_thickness_measurement : finished with event loop" << std::endl;

    // If we're using 228Th, the thoron leakage can potentially cause misidentification problems (e.g., peaks still show up)
    // even with sufficiently thick targets, meaning that misidentification can frequently occur. To mediate this, loop through
    // the energy histograms and remember what the maximum in any one is, then determine if analyses should take place based on
    // the ratio of the current energy histogram to the maximum. 40% is probably a reasonable cutoff.
    int maximum_number_entries {};
    for(auto i : histograms) {
        if(i->GetEntries() > maximum_number_entries) maximum_number_entries = i->GetEntries();
    }

    // Add the first thickness histogram to the vector. This is the one that contains all thickness values, good or bad.
    thickness_maps.push_back(new TH2F("h_ThicknessMap_0", "", num_divisions, -1, 1, num_divisions, -1, 1));
    
    // With many spectra, begin looping through them and performing the analysis.
    std::cout << "perform_thickness_measurement : beginning iterative peak search" << std::endl;
    for(auto i = 0; i < histograms.size(); ++i) {
        if(i % 10 == 0) printf("Fit %d/%d\n", i, histograms.size());

        if(histograms.at(i)->GetEntries() > 500 && histograms.at(i)->GetEntries() > 0.4 * maximum_number_entries) {
            if(histograms.at(i)->GetEntries() < 2000) histograms.at(i)->Rebin(4);
            // Set up the TSpectrum to identify peaks.
            if (verbosity > 2) std::cout << "----- POINT 1 -----" << std::endl;
            auto *spec = new TSpectrum(2 * source_energies.size());
            auto threshold {0.25}; // determines which peaks to cut off for consideration.
                                   // There is a chance that the fit won't converge on the lowest shoulder for
                                   // 228Th, but the thickness will only be gauged with higher peaks, so no concern.
            auto peak_list_length = spec->Search(histograms.at(i), 5, "nobackground goff", threshold);
            double *x_peaks = spec->GetPositionX();
            double par[2 * source_energies.size() + 1];
            par[0] = 5; // universal width of the peaks, in channels
            for(auto j = 0; j < peak_list_length; ++j) {
                par[2*j + 1] = histograms.at(i)->GetBinContent(histograms.at(i)->GetXaxis()->FindBin(x_peaks[j])); // Set the approx. norm.
                par[2*j + 2] = x_peaks[j]; // Set the mean

            }

            if(verbosity > 2) std::cout << "----- POINT 2 -----" << std::endl;
            // Fit the peaks that were found.
            auto *f_spectrum = new TF1("f_spectrum", spectrum, 0, 2000, 2*peak_list_length + 1);
            f_spectrum->SetNpx(8000);
            f_spectrum->SetParameters(par);
            for(auto j = 0; j < 2 * peak_list_length + 1; ++j) {
                f_spectrum->SetParLimits(j, f_spectrum->GetParameter(j) * 0.75, f_spectrum->GetParameter(j) * 1.25);
            }
            histograms.at(i)->Fit(f_spectrum, "EQ+");
            auto *c1 = (TCanvas*)gROOT->FindObjectAny("c1");
            if(c1) delete c1; // Because TSpectrum can't handle the canvas it creates...
            
            if(verbosity > 2) std::cout << "----- POINT 3 -----" << std::endl;
            // Extract the peaks, then sort in descending order
            std::vector<std::pair<double, int>> peak_centroids {};
            for(auto j = 0; j < (sizeof(par)/sizeof(double) - 1)/2; ++j) {
                if(par[2 + 2*j] > 5) {
                    peak_centroids.push_back(std::make_pair(f_spectrum->GetParameter(2 + 2*j), 2 + 2*j));
                }
            }
            std::sort(peak_centroids.begin(), peak_centroids.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
            
            if(verbosity > 2) std::cout << "----- POINT 4 -----" << std::endl;
            // Convert the channel values to energy using the energy calibration parameters. Also save centroids for error calculation
            std::vector<double> upper_peak_centroids {}, lower_peak_centroids {};
            auto error {std::sqrt(std::pow(f_spectrum->GetParError(peak_centroids.at(0).second), 2.0) + 
                        std::pow(f_spectrum->GetParError(peak_centroids.at(1).second), 2.0))};
            for(auto j = 0; j < peak_centroids.size(); ++j) {
                upper_peak_centroids.push_back(peak_centroids.at(j).first + error);
                lower_peak_centroids.push_back(peak_centroids.at(j).first - error);
				
				if(verbosity > 1) {
					std::cout << "peak (channel) index: " << j << "; centroid: " << peak_centroids.at(j).first << "; upper: " << upper_peak_centroids.back() 
							  << "; lower: " << lower_peak_centroids.back() << std::endl;
				}		
	
                peak_centroids.at(j).first = energy_slope * peak_centroids.at(j).first + energy_intercept;
                upper_peak_centroids.back() = energy_slope * upper_peak_centroids.back() + energy_intercept;
                lower_peak_centroids.back() = energy_slope * lower_peak_centroids.back() + energy_intercept;

				if(verbosity > 1) {
					std::cout << "peak (MeV) index: " << j << "; centroid: " << peak_centroids.at(j).first << "; upper: " << upper_peak_centroids.back()
							  << "; lower: " << lower_peak_centroids.back() << std::endl;
				}
            }

            if(verbosity > 2) std::cout << "----- POINT 5 -----" << std::endl;
            // CycSrim -> Get the energy loss from residual, then give to the range function.
            // Using the CycSrim range functionality, manually calculate the energy differences 
            // between the observed spectra and the no-target values.
            std::vector<double> nominal_energy_losses {}, thicknesses {};
            std::vector<double> lower_energy_losses {}, upper_energy_losses {}, lower_errors {}, upper_errors {};
            auto num_peaks_to_consider {is_228Th ? 3 : 5}; // If 228Th, use three peaks; if 229Th, use 5.
            for(auto j = 0; j < std::min(static_cast<int>(peak_centroids.size()), num_peaks_to_consider); ++j) {
                // We must account for the dead layer impact in our observation of the energy so that we may compare it to the
                // source energies prior to the dead layer. This is what matters for the thickness calculation. We are also 
                // considering thickness error by adjusting the energies with the fit uncertainties.
                nominal_energy_losses.push_back(source_energies.at(j) - (peak_centroids.at(j).first + dead_layer->GetEnergyLossFromResidual(2, 4, peak_centroids.at(j).first)));
                upper_energy_losses.push_back(source_energies.at(j) - (lower_peak_centroids.at(j) + dead_layer->GetEnergyLossFromResidual(2, 4, lower_peak_centroids.at(j))));
                lower_energy_losses.push_back(source_energies.at(j) - (upper_peak_centroids.at(j) + dead_layer->GetEnergyLossFromResidual(2, 4, upper_peak_centroids.at(j))));
                
				if(verbosity > 1) {
					std::cout << "peak index: " << j << "; nominal loss (MeV): " << nominal_energy_losses.back() << "; upper loss (MeV): "
							  << upper_energy_losses.back() << "; lower loss (MeV): " << lower_energy_losses.back() << std::endl;
				}

                // If a peak is misidentified and the energy loss is calculated to be a negative number, ignore it and move on.
                if(nominal_energy_losses.back() > 0) {
                    // Calculate thickness as range_originalE - range_reducedE; returns thickness in units specified by TARGET object.
                    // Same as above, need to account for dead layer, but not necessary since handled earlier in script.
                    thicknesses.push_back(TARGET->GetRange(2, 4, source_energies.at(j)) -
                                        TARGET->GetRange(2, 4, source_energies.at(j) - nominal_energy_losses.back()));

                    // We'll account for errors by calculating the thickness using the upper and lower gauges on the energy spectrum fits.
                    upper_errors.push_back(TARGET->GetRange(2, 4, source_energies.at(j)) - TARGET->GetRange(2, 4, source_energies.at(j) - upper_energy_losses.back()));
                    lower_errors.push_back(TARGET->GetRange(2, 4, source_energies.at(j)) - TARGET->GetRange(2, 4, source_energies.at(j) - lower_energy_losses.back()));
			
					if(verbosity > 1) {
						std::cout << "peak index: " << j << "; thickness: " << thicknesses.back() << "; upper error: " << upper_errors.back() 
								  << "; lower error: " << lower_errors.back() << std::endl;
					}
                }
            }
            
            if(verbosity > 2) std::cout << "----- POINT 6 -----" << std::endl;
            // Average the thickness array to get a value
            double avg_thickness {}, avg_upper_error {}, avg_lower_error {};
            if(!thicknesses.empty()) {
                avg_thickness = std::accumulate(thicknesses.begin(), thicknesses.end(), 0.0) / thicknesses.size();
                avg_upper_error = std::accumulate(upper_errors.begin(), upper_errors.end(), 0.0) / upper_errors.size();
                avg_lower_error = std::accumulate(lower_errors.begin(), lower_errors.end(), 0.0) / lower_errors.size();

				if(verbosity > 0) {
					std::cout << "avg. thickness: " << avg_thickness << "; avg. upper error: " << avg_upper_error << "; avg. lower error: " << avg_lower_error << std::endl;
				}
            }

            // Add the thickness value to the 2D representative histogram.
            auto x_bin {static_cast<int>(i % num_divisions) + 1}; // ROOT counts from 1
            auto y_bin {static_cast<int>(i / num_divisions) + 1}; // ROOT counts from 1
            thickness_maps.back()->SetBinContent(x_bin, y_bin, avg_thickness);

	    	h_UpperError->SetBinContent(x_bin, y_bin, avg_upper_error - avg_thickness);
	    	h_LowerError->SetBinContent(x_bin, y_bin, avg_thickness - avg_lower_error);
            
            if(verbosity > 2) std::cout << "----- POINT 7 -----" << std::endl;
            // Reclaim memory?
        }
    }

    // In general, using a source that doesn't have leakage will not produce a thickness map that requires corrective
    // action. In this case, if the source being used is not 228Th, the iterative correction part of the code is bypassed.
    if(is_228Th) {
        // Loop through the thickness histogram and calculate what the average thickness is. This can be a very poor approximation
        // of the true average thickness, especially if the number of bad bins is large. Use this value to filter out the obviously 
        // incorrect bins with a forgiving cut, then re-loop through the histogram to calculate a better approximation of the 
        // thickness. This is what we'll use to more accurately identify bad fit bins.
        int num_bins_summed {};
        double poor_average_value {};
        for(auto x = 1; x <= thickness_maps.back()->GetNbinsX(); ++x) {
            for(auto y = 1; y <= thickness_maps.back()->GetNbinsY(); ++y) {
                auto cont {thickness_maps.back()->GetBinContent(x, y)};
                if(cont > 0.005 && cont < 1000.0) {
                    poor_average_value += cont;
                    ++num_bins_summed;
                }
            }
        }

        poor_average_value /= num_bins_summed;

        double average_value {};
        num_bins_summed = 0;
        for(auto x = 1; x <= thickness_maps.back()->GetNbinsX(); ++x) {
            for(auto y = 1; y <= thickness_maps.back()->GetNbinsY(); ++y) {
                auto cont {thickness_maps.back()->GetBinContent(x, y)};
                if(cont > 0.1 * poor_average_value && cont < 3.0 * poor_average_value) {
                    average_value += cont;
                    ++num_bins_summed;
                }
            }
        }

        average_value /= num_bins_summed;

        // Add a second thickness map to the vector. This is the first filtered thickness map, produced by gating
        // using arbitrary thresholds.
        thickness_maps.push_back(new TH2F("h_ThicknessMap_1", "", num_divisions, -1, 1, num_divisions, -1, 1));

        std::vector<std::pair<int, int>> poor_fits {};

        double minimum_value {1E6}, maximum_value {1E-10};
        for(auto x = 1; x <= thickness_maps.back()->GetNbinsX(); ++x) {
            for(auto y = 1; y <= thickness_maps.back()->GetNbinsY(); ++y) {
                auto cont {thickness_maps.at(thickness_maps.size() - 2)->GetBinContent(x, y)};
                // We frequently misidentify the highest energy peak
                if(cont > 0.25 * average_value && cont < 1.5 * average_value) {
                    thickness_maps.back()->SetBinContent(x, y, cont);
                    
                    if(cont > maximum_value) maximum_value = cont;
                    if(cont < minimum_value) minimum_value = cont;
                } else if(cont != 0) {
                    poor_fits.push_back(std::make_pair(x, y));
                }
            }
        }
        thickness_maps.back()->GetZaxis()->SetRangeUser(0.5 * minimum_value, 2.0 * maximum_value);

        // With a list of bad fits, recalculate the thicknesses iteratively. Recalculate missing bins using the nearest
        // two or more good bins as an estimate for the thickness, fit with Gaussians, fill histogram. Limit the number
        // of iterations to prevent problems if an unanticipated problem arises. This is done by comparing the number of
        // filled bins in the final iteration to that in the uncorrected map - they should have the same number present.

        int number_trials {}, bins_modified_this_trial {};
        while(true) {
            if(number_trials > 5) break;

            thickness_maps.push_back(static_cast<TH2F *>(thickness_maps.back()->Clone()));
            // thickness_maps.push_back(new TH2F(Form("h_ThicknessMaps_%d", number_trials + 2), "", num_divisions, -1, 1, num_divisions, -1, 1));
            thickness_maps.back()->SetName(Form("h_ThicknessMap_%d", number_trials + 2));

            bins_modified_this_trial = 0;

            // Loop through bad bin IDs. If the value is non-zero, it's been evaluated, so pass.
            for(auto i: poor_fits) {
                if(thickness_maps.back()->GetBinContent(i.first, i.second) == 0) {
                    double average_thickness_bins {};
                    int num_added {};

                    // Check number of filled bins around trial bin; if not two or greater, break and evaluate later. Make sure that the bins
                    // being added aren't the underflow or overflow bins.
                    auto value {thickness_maps.back()->GetBinContent(i.first, i.second + 1)}; // Above
                    if(value > 0 && i.second != thickness_maps.back()->GetNbinsY()) {
                        average_thickness_bins += value; ++num_added;
                    }
                    value = thickness_maps.back()->GetBinContent(i.first - 1, i.second); // Left
                    if(value > 0 && i.first != 1) {
                        average_thickness_bins += value; ++num_added;
                    }
                    value = thickness_maps.back()->GetBinContent(i.first + 1, i.second); // Right
                    if(value > 0 && i.first != thickness_maps.back()->GetNbinsX()) {
                        average_thickness_bins += value; ++num_added;
                    }
                    value = thickness_maps.back()->GetBinContent(i.first, i.second - 1); // Below
                    if(value > 0 && i.second != 1) {
                        average_thickness_bins += value; ++num_added;
                    }

                    if(num_added > 1) {
                        average_thickness_bins /= num_added;

                        TARGET->SetThickness(average_thickness_bins);
                        // Perform two Gaussian discrete fit
                        auto mean_chan {(TARGET->GetResidualEnergy(2, 4, source_energies.front()) - energy_intercept) / energy_slope}; // only approximate
                        auto *f_peak_1 {new TF1("f_peak_1", "gaus", mean_chan - 50, mean_chan + 50)};
                        f_peak_1->SetParameters(10, mean_chan, 6);
                        f_peak_1->SetParLimits(1, 0.8 * mean_chan, 1.2 * mean_chan);
                        histograms.at(num_divisions * (i.second - 1) + (i.first - 1))->Fit(f_peak_1, "BQR+"); // TH2 count from 1, so adjust; Fit only specified range
                        mean_chan = (TARGET->GetResidualEnergy(2, 4, source_energies.at(1)) - energy_intercept) / energy_slope;
                        auto *f_peak_2 {new TF1("f_peak_2", "gaus", mean_chan - 50, mean_chan + 50)};
                        f_peak_2->SetParameters(10, mean_chan, 6);
                        f_peak_2->SetParLimits(1, 0.8 * mean_chan, 1.2 * mean_chan);
                        histograms.at(num_divisions * (i.second - 1) + (i.first - 1))->Fit(f_peak_2, "BQR+"); // TH2 count from 1, so adjust; Fit only specified range

                        // With the fits, calculate the thickness as done in the previous full loop
                        auto detected_energy_1 {energy_slope * f_peak_1->GetParameter(1) + energy_intercept};
                        auto detected_energy_2 {energy_slope * f_peak_2->GetParameter(1) + energy_intercept};
                        auto energy_loss_1 {source_energies.front() - (detected_energy_1 + dead_layer->GetEnergyLossFromResidual(2, 4, detected_energy_1))};
                        auto energy_loss_2 {source_energies.at(1) - (detected_energy_2 + dead_layer->GetEnergyLossFromResidual(2, 4, detected_energy_2))};
                        auto thickness_1 {TARGET->GetRange(2, 4, source_energies.front()) - TARGET->GetRange(2, 4, source_energies.front() - energy_loss_1)};
                        auto thickness_2 {TARGET->GetRange(2, 4, source_energies.at(1)) - TARGET->GetRange(2, 4, source_energies.at(1) - energy_loss_2)};
                        auto thickness {(thickness_1 + thickness_2) / 2.0};

                        if(thickness > 1.05 * average_thickness_bins || thickness < 0.95 * average_thickness_bins) {
                            thickness_maps.back()->SetBinContent(i.first, i.second, average_thickness_bins);
                        } else {
                            thickness_maps.back()->SetBinContent(i.first, i.second, thickness);
                        }

                        ++bins_modified_this_trial;
                    }
                }
            } // Poor fits
        
            if(bins_modified_this_trial == 0) break; // Checking number of entries in histograms unreliable for some reason.

            ++number_trials;
        }
    }
    
    h_ThicknessMap = (TH2F *)thickness_maps.back()->Clone();

    // To make viewing of this histogram as efficient as possible
    auto minimum_cont {std::numeric_limits<float>::max()}, maximum_cont {0.0f};
    for(int x = 1; x <= h_ThicknessMap->GetNbinsX(); ++x) {
       for(int y = 1; y <= h_ThicknessMap->GetNbinsY(); ++y) {
	  auto cont {h_ThicknessMap->GetBinContent(x, y)};
	  if(cont < minimum_cont && cont != 0) minimum_cont = cont;
	  if(cont> maximum_cont) maximum_cont = cont;
       }
    }

    h_ThicknessMap->SetName("h_ThicknessMap");
    h_ThicknessMap->GetZaxis()->SetRangeUser(minimum_cont, maximum_cont);

    auto *root_dir {gDirectory};
    auto *diagnostic_dir {gDirectory->mkdir("Diagnostics")};
    auto *source_spectra_dir {gDirectory->mkdir("Source Spectra")};

    h_ThicknessMap->Write();

    diagnostic_dir->cd();
    h_LowerError->Write();
    h_UpperError->Write();
    h_HitMap_Uncorrected->Write();
    h_HitMap->Write();
    h_Energy->Write();
    for(auto i : thickness_maps) i->Write();
    root_dir->cd();
    source_spectra_dir->cd();
    for(auto i : histograms) i->Write();
    root_dir->cd();

    out_file->Close();

    std::cout << "perform_thickness_measurement : finished" << std::endl;
    std::cout << std::endl;
    std::cout << "Please refer to the README for descriptions of each histogram in the file." << std::endl;

    return 0;
}