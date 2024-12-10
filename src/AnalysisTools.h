#ifndef ANALYSIS_TOOLS__H
#define ANALYSIS_TOOLS__H

#include <fstream>
#include <string>
#include <sstream>

#include "PathManager.h"

//___________________________________________________________________________________________________
// Global properties of the testing station
// std::vector<double> source_energies {8.785, 6.778, 6.288, 5.685, 5.423, 5.340}; // Global 228Th values
std::vector<double> source_energies {8.38, 7.07, 6.34, 5.83, 4.85}; // Global 229Th values

float source_detector_distance {12.7}; // Distance between source and det face, measured in cm.
//___________________________________________________________________________________________________


// This integrator method is the old version that should be compatible with all event types. Newer
// versions of the event have streamlined functions and methods.
void integrator_method(T063123Event *event, double &f1, double &f2, double &b1, double &b2) {
    // Decide which time to use; ideally all DADL signals are received
    double f1Time_u {event->GetF1Time(1) / 1000.0}, f2Time_u {event->GetF2Time(1) / 1000.0};
    double b1Time_u {event->GetB1Time(1) / 1000.0}, b2Time_u {event->GetB2Time(1) / 1000.0};

    // Initialize the correct integrator values
    double IdealFrontTime {}, IdealBackTime {};

    double IntegratorValueF1[8];
    double IntegratorStartTimeF1[8];
    double IntegratorValueF2[8];
    double IntegratorStartTimeF2[8];
    double IntegratorValueB1[8];
    double IntegratorStartTimeB1[8];
    double IntegratorValueB2[8];
    double IntegratorStartTimeB2[8];

    double IntegratorDiffF1[8];
    double IntegratorDiffF2[8];
    double IntegratorDiffB1[8];
    double IntegratorDiffB2[8];

    // Choose the best integrator based on ADC timestamps. BEGINNING OF INTEGRATOR ALGORITHM
    // Use integrator 1 and integrator baseline as rough estimates.
    double f1_t {event->F1Fired(1) ? event->GetF1Integrator(1, 1) : 0};
    double f2_t {event->F2Fired(1) ? event->GetF2Integrator(1, 1) : 0};
    double b1_t {event->B1Fired(1) ? event->GetB1Integrator(1, 1) : 0};
    double b2_t {event->B2Fired(1) ? event->GetB2Integrator(1, 1) : 0};

    if(f1_t > 2 || f2_t > 2 || b1_t > 2 || b2_t > 2) {
        Bool_t useF1 {false}, useF2 {false}, useB1 {false}, useB2 {false};

        if(f1Time_u <= f2Time_u) useF1 = true;
        if(f2Time_u < f1Time_u) useF2 = true;
        if(b1Time_u <= b2Time_u) useB1 = true;
        if(b2Time_u < b1Time_u) useB2 = true;

        if(useF1) IdealFrontTime = f1Time_u + 1.6;
        if(useF2) IdealFrontTime = f2Time_u + 1.6;
        if(useB1) IdealBackTime = b1Time_u + 1.6;
        if(useB2) IdealBackTime = b2Time_u + 1.6;

        // Initialize array that has the true timestamps for each integrator start.
        for(Int_t j = 0; j < 7; ++j) {
            IntegratorStartTimeF1[j + 1] = f1Time_u + 1.6 - (0.088) * j;
            IntegratorStartTimeF2[j + 1] = f2Time_u + 1.6 - (0.088) * j;
            IntegratorStartTimeB1[j + 1] = b1Time_u + 1.6 - (0.088) * j;
            IntegratorStartTimeB2[j + 1] = b2Time_u + 1.6 - (0.088) * j;
        }

        // Set baselines, then add the other 7 integrator values to the array
        IntegratorValueF1[0] = event->GetF1Baseline(1);
        IntegratorValueF2[0] = event->GetF2Baseline(1);
        IntegratorValueB1[0] = event->GetB1Baseline(1);
        IntegratorValueB2[0] = event->GetB2Baseline(1);
        for(Int_t j = 1; j < 8; ++j) {
            IntegratorValueF1[j] = event->GetF1Integrator(1, j);
            IntegratorValueF2[j] = event->GetF2Integrator(1, j);
            IntegratorValueB1[j] = event->GetB1Integrator(1, j);
            IntegratorValueB2[j] = event->GetB2Integrator(1, j);
        }

        // Set high value for the baseline integrator difference so it is never picked
        IntegratorDiffF1[0] = std::numeric_limits<double>::max();
        IntegratorDiffF2[0] = std::numeric_limits<double>::max();
        IntegratorDiffB1[0] = std::numeric_limits<double>::max();
        IntegratorDiffB2[0] = std::numeric_limits<double>::max();

        // Calculate the difference between the integrator start time and the ideal front/back time
        for(int j = 1; j <= 7; ++j) {
            IntegratorDiffF1[j] = abs(IntegratorStartTimeF1[j] - IdealFrontTime);
            IntegratorDiffF2[j] = abs(IntegratorStartTimeF2[j] - IdealFrontTime);
            IntegratorDiffB1[j] = abs(IntegratorStartTimeB1[j] - IdealBackTime);
            IntegratorDiffB2[j] = abs(IntegratorStartTimeB2[j] - IdealBackTime);
        }

        // Choose the integrator closest to the best time and use that one
        int f1Integ = std::distance((IntegratorDiffF1), std::min_element(IntegratorDiffF1 + 1, IntegratorDiffF1 + 8));
        int f2Integ = std::distance((IntegratorDiffF2), std::min_element(IntegratorDiffF2 + 1, IntegratorDiffF2 + 8));
        int b1Integ = std::distance((IntegratorDiffB1), std::min_element(IntegratorDiffB1 + 1, IntegratorDiffB1 + 8));
        int b2Integ = std::distance((IntegratorDiffB2), std::min_element(IntegratorDiffB2 + 1, IntegratorDiffB2 + 8));

        f1_t = IntegratorValueF1[f1Integ] - IntegratorValueF1[0];
        f2_t = IntegratorValueF2[f2Integ] - IntegratorValueF2[0];
        b1_t = IntegratorValueB1[b1Integ] - IntegratorValueB1[0];
        b2_t = IntegratorValueB2[b2Integ] - IntegratorValueB2[0];

        // Divide by number of bins integrated
        f1 = static_cast<double>(f1_t) / 500.0;
        f2 = static_cast<double>(f2_t) / 500.0;
        b1 = static_cast<double>(b1_t) / 500.0;
        b2 = static_cast<double>(b2_t) / 500.0;
    }
}

void file_not_open(const char *file) {
    std::cout << "AnalysisTools::file_not_open - missing parameter file: " << std::endl;
    std::cout << file << std::endl;
    std::cout << "AnalysisTools::file_not_open - check calibrations directory or consult README." << std::endl;
    std::cout << std::endl;
    std::abort();
}

void read_face_parameters(const char *file, double &slope, double &offset) {
    std::fstream in_file(file, std::ios_base::in);
    if(!in_file.is_open()) file_not_open(file);
    std::string in_line {}; std::getline(in_file, in_line);
    std::istringstream buffer(in_line); buffer >> slope >> offset;
    in_file.close();
}

void read_energy_parameters(const char *file, double &slope, double &intercept) {
    read_face_parameters(file, slope, intercept);
}

void read_stretching_parameters(const char *file, double &top, double &bottom, double &right, double &left) {
    std::fstream in_file(file, std::ios_base::in);
    if(!in_file.is_open()) file_not_open(file);
    std::string in_line {}; std::getline(in_file, in_line);
    std::istringstream buffer(in_line);
    buffer >> top >> bottom >> right >> left;
    in_file.close();
}

bool front_equal_back(double front, double back, double slope, double offset) {
    // F = s * B + o; assume that the offset doesn't play a significant role for this
    double upper_limit {back * (1.1*slope)}, lower_limit {back * (0.9 * slope)};
    if(front < upper_limit && front > lower_limit) return true;
    return false;
}

double gain_match(double sig_1, double sig_2, double slope, double intercept) {
    double result = (sig_1 + ((intercept / (2.0 * slope)) * (sig_1 / (sig_1 + sig_2)))) * ((slope + 1.0) / 2.0) + 
                    (sig_2 - ((intercept / 2.0) * (sig_2 / (sig_1 + sig_2)))) * ((slope + 1.0) / (2.0 * slope));
    return result;
}

// For X, edge_1 = right, edge_2 = left; for Y, edge_1 = top, edge_2 = bottom
double stretch_position(double sig_1, double sig_2, double edge_1, double edge_2) {
    double local {(sig_2 - sig_1) / (sig_1 + sig_2)};
    local = local - ((edge_1 + edge_2) / 2.0); // Shift position to center it
    local = local * (2.0 / (edge_1 - edge_2)); // Stretch position after centering shift
    return local;
}

double spectrum(double *x, double *p) {
    double result{0};
    double sigm = p[0];
    for(int i = 0; i < static_cast<int>(source_energies.size()); ++i) {
        double norm = p[2*i + 1];
        double mean = p[2*i + 2];

        result += norm*TMath::Gaus(x[0], mean, sigm);
    }
    return result;
}

#endif
