// This executable creates a ROOT histogram file filled with various summary plots.

// Other includes
#include "PathManager.h"
#include "AnalysisTools.h"

int main(int argc, char **argv) {
    // option manager
    BrAppOptionManager::Instance()->SetVersion(1, 0, "generate_raw_summary_file");
    BrAppOptionManager::Instance()->SetHelp("Produces summary histogram file for given run.");
    BrAppOptionManager::Instance()->SetCommandLine(argc, argv);
    auto *opt_run_number = new BrAppIntOption('r', "run_number", "run number", 0);
    auto *opt_events     = new BrAppIntOption('e', "events", "events to analyze", std::numeric_limits<int>::max());

    if(!BrAppOptionManager::Instance()->ProcessCommandLine()) { return 1; }
    if(BrAppOptionManager::Instance()->ShowVersion()) { return 0; }
    if(BrAppOptionManager::Instance()->ShowHelp()) { return 0; }



    // Reduced file setup and tree read-in
    std::cout << std::endl;
    std::cout << "generate_raw_summary_file : reduced file read-in" << std::endl;
    auto in_file = std::make_unique<TFile>(Form(Form("%s%s", REDUCED_FILE_DIR, FILE_NAME_FORMAT), EXPT_DATE, opt_run_number->GetValue()));
    gROOT->cd();
    auto tree = dynamic_cast<TTree *>(in_file->FindObjectAny("T"));
    T063123Event *evin {};
    tree->SetBranchAddress("event", &evin);
    auto num_entries {std::min(static_cast<int>(tree->GetEntries()), opt_events->GetValue())};

    // Output file and histograms
    auto out_file = std::make_unique<TFile>(Form("%ssummary_%s%03d.root", OUTPUT_FILE_DIR, EXPT_DATE, opt_run_number->GetValue()), "RECREATE");
    // Signals
    auto channel_limit {2000};
    auto h_F1 = new TH1F("h_F1", "h_F1", 1000, 0, channel_limit);
    auto h_F2 = new TH1F("h_F2", "h_F2", 1000, 0, channel_limit);
    auto h_B1 = new TH1F("h_B1", "h_B1", 1000, 0, channel_limit);
    auto h_B2 = new TH1F("h_B2", "h_B2", 1000, 0, channel_limit);
    // Energy Spectra
    auto h_EF_ch = new TH1F("h_EF_ch", "h_EF_ch", 2000, 0, channel_limit);
    auto h_EB_ch = new TH1F("h_EB_ch", "h_EB_ch", 2000, 0, channel_limit);
    auto h_ES_ch = new TH1F("h_ES_ch", "h_ES_ch", 2000, 0, channel_limit);
    auto h_EFvsY = new TH2F("h_EFvsY", "h_EFvsY", 100, -1, 1, 1000, 0, channel_limit);
    auto h_EBvsX = new TH2F("h_EBvsX", "h_EBvsX", 100, -1, 1, 1000, 0, channel_limit);
    auto h_FSumDif = new TH2F("h_FSumDif", "h_FSumDif", 1000, -channel_limit/2.0, channel_limit/2.0, 1000, 0, channel_limit);
    auto h_BSumDif = new TH2F("h_BSumDif", "h_BSumDif", 1000, -channel_limit/2.0, channel_limit/2.0, 1000, 0, channel_limit);
    auto h_F1vsF2 = new TH2F("h_F1vsF2", "h_F1vsF2", 500, 0, channel_limit, 500, 0, channel_limit);
    auto h_B1vsB2 = new TH2F("h_B1vsB2", "h_B1vsB2", 500, 0, channel_limit, 500, 0, channel_limit);
    auto h_FvsB = new TH2F("h_FvsB", "h_FvsB", 500, 0, channel_limit, 500, 0, channel_limit);
    // Position Spectra
    auto h_X = new TH1F("h_X", "h_X", 200, -1, 1);
    auto h_Y = new TH1F("h_Y", "h_Y", 200, -1, 1);
    auto h_XY = new TH2F("h_XY", "h_XY", 200, -1, 1, 200, -1, 1);


    // Event loop
    std::cout << "generate_raw_summary_file : entering event loop" << std::endl;
    for(auto i = 0; i < num_entries; ++i) {
        if(i % 10000 == 0) printf("Event %d/%d\n", i, num_entries);
        tree->GetEntry(i);

        double f1 {}, f2 {}, b1 {}, b2 {};
        integrator_method(evin, f1, f2, b1, b2);

        if(f1 > 5) h_F1->Fill(f1);
        if(f2 > 5) h_F2->Fill(f2);
        if(b1 > 5) h_B1->Fill(b1);
        if(b2 > 5) h_B2->Fill(b2);

        double front_sum {}, back_sum {}, x_pos {}, y_pos {};

        if(f1 > 5 && f2 > 5) {
            front_sum = f1 + f2;
            y_pos = (f2 - f1) / (f1 + f2);

            h_EF_ch->Fill(front_sum);
            h_EFvsY->Fill(y_pos, front_sum);
            h_F1vsF2->Fill(f2, f1);
            h_FSumDif->Fill(f2 - f1, front_sum);
            h_Y->Fill(y_pos);
        }

        if(b1 > 5 && b2 > 5) {
            back_sum = b1 + b2;
            x_pos = (b2 - b1) / (b1 + b2);

            h_EB_ch->Fill(back_sum);
            h_EBvsX->Fill(x_pos, back_sum);
            h_B1vsB2->Fill(b2, b1);
            h_BSumDif->Fill(b2 - b1, back_sum);
            h_X->Fill(x_pos);
        }

        if(f1 > 5 && f2 > 5 && b1 > 5 && b2 > 5) {
            h_FvsB->Fill(back_sum, front_sum);
            h_XY->Fill(x_pos, y_pos);
            h_ES_ch->Fill((front_sum + back_sum) / 2.0);
        }
    }

    out_file->Write();
    out_file->Close();
    in_file->Close();

    std::cout << "generate_raw_summary_file : finished" << std::endl;
    std::cout << std::endl;

    return 0;
}