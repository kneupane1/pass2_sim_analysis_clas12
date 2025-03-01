#include <fstream>
#include <string>
#include "TFile.h"
#include "THnSparse.h"
#include <iostream>

static const int W_bins_no = 15;
static const int Q2_bins_no = 8;

float q2_low_values[Q2_bins_no] = {2.0, 2.40, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0};
float q2_up_values[Q2_bins_no] = {2.4, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0};

THnSparseD *h_simu_prot[W_bins_no];
THnSparseD *h_simu_prot_background[W_bins_no];

THnSparseD *h_simu_pip[W_bins_no];
THnSparseD *h_simu_pip_background[W_bins_no];

THnSparseD *h_simu_pim[W_bins_no];
THnSparseD *h_simu_pim_background[W_bins_no];

THnSparseD *h_simu_thrown_prot[W_bins_no];
THnSparseD *h_simu_thrown_pip[W_bins_no];
THnSparseD *h_simu_thrown_pim[W_bins_no];

THnSparseD *h_simu_prot_evt[W_bins_no];
THnSparseD *h_simu_prot_background_evt[W_bins_no];

THnSparseD *h_simu_pip_evt[W_bins_no];
THnSparseD *h_simu_pip_background_evt[W_bins_no];

THnSparseD *h_simu_pim_evt[W_bins_no];
THnSparseD *h_simu_pim_background_evt[W_bins_no];

THnSparseD *h_simu_thrown_prot_evt[W_bins_no];
THnSparseD *h_simu_thrown_pip_evt[W_bins_no];
THnSparseD *h_simu_thrown_pim_evt[W_bins_no];

int main()
{
    // // // File paths
    std::string mcFileName = "/Users/krishnaneupane/Downloads/2025/For_background_subtraction/bkg_files/new/pass2_all_741_files_for_final_cs.root";
    std::string background_mcFileName = "/Users/krishnaneupane/Downloads/2025/For_background_subtraction/bkg_files/new/pass2_final_all_741_files_for_background_subtraction.root";

    // ////File paths
    // std::string background_mcFileName = "/Users/krishnaneupane/Downloads/2025/For_background_subtraction/bkg_files/new/pass2_first_741_files_for_background_subtraction_wt_1.root";
    // std::string mcFileName = "/Users/krishnaneupane/Downloads/2025/Used_for_cs/out_Pass2_sim_twoPi_rga_fall2018_tor-1_sol-1_flagrad_2_all_741_files_with_wt_1.root";

    // Open ROOT files
    TFile *root_mc = new TFile(mcFileName.c_str(), "READ");
    TFile *root_mc_background = new TFile(background_mcFileName.c_str(), "READ");

    if (!root_mc || !root_mc->IsOpen())
    {
        std::cerr << "Error: Could not open MC file!" << std::endl;
        return 1;
    }
    if (!root_mc_background || !root_mc_background->IsOpen())
    {
        std::cerr << "Error: Could not open background MC file!" << std::endl;
        return 1;
    }

    // Output ROOT file
    TFile *outputFile = new TFile("/Users/krishnaneupane/Downloads/2025/For_background_subtraction/bkg_files/new/background_subtracted_741_files_rec_mc_thnsparse_hists_new1.root", "RECREATE");
    // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_thrown_prot");
    outputFile->cd("THnSparse_7D_thrown_prot");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_simu = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", q2_lower_lim, q2_upper_lim, 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_thrown_prot[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_thrown_prot/%s", name_simu));
            // Write to output file
            h_simu_thrown_prot[w]->Write();
        }
    }
    // // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_thrown_pim");
    outputFile->cd("THnSparse_7D_thrown_pim");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_simu = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", q2_lower_lim, q2_upper_lim, 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_thrown_pim[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_thrown_pim/%s", name_simu));
            h_simu_thrown_pim[w]->Write();
        }
    }
    // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_thrown_pip");
    outputFile->cd("THnSparse_7D_thrown_pip");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_simu = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", q2_lower_lim, q2_upper_lim, 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_thrown_pip[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_thrown_pip/%s", name_simu));
            // Write to output file
            h_simu_thrown_pip[w]->Write();
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_prot");
    outputFile->cd("THnSparse_7D_prot");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_simu = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", q2_lower_lim, q2_upper_lim, 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_prot[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_prot/%s", name_simu));
            h_simu_prot_background[w] = (THnSparseD *)root_mc_background->Get(Form("THnSparse_7D_prot/%s", name_simu));
            h_simu_prot[w]->Add(h_simu_prot_background[w], -1.0);

            // Write to output file
            h_simu_prot[w]->Write();
        }
    }
    // // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_pim");
    outputFile->cd("THnSparse_7D_pim");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_simu = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", q2_lower_lim, q2_upper_lim, 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_pim[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_pim/%s", name_simu));
            h_simu_pim_background[w] = (THnSparseD *)root_mc_background->Get(Form("THnSparse_7D_pim/%s", name_simu));
            h_simu_pim[w]->Add(h_simu_pim_background[w], -1.0);

            // Write to output file

            h_simu_pim[w]->Write();
        }
    }
    // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_pip");
    outputFile->cd("THnSparse_7D_pip");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_simu = Form("h_5dim_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", q2_lower_lim, q2_upper_lim, 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_pip[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_pip/%s", name_simu));
            h_simu_pip_background[w] = (THnSparseD *)root_mc_background->Get(Form("THnSparse_7D_pip/%s", name_simu));
            h_simu_pip[w]->Add(h_simu_pip_background[w], -1.0);

            // Write to output file
            h_simu_pip[w]->Write();
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_prot_evt");
    outputFile->cd("THnSparse_7D_prot_evt");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (q2_lower_lim), (q2_upper_lim), 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_prot_evt[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_prot_evt/%s", name_evt));
            h_simu_prot_background_evt[w] = (THnSparseD *)root_mc_background->Get(Form("THnSparse_7D_prot_evt/%s", name_evt));
            h_simu_prot_evt[w]->Add(h_simu_prot_background_evt[w], -1.0);
            h_simu_prot_evt[w]->Write();
        }
    }
    // // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_pim_evt");
    outputFile->cd("THnSparse_7D_pim_evt");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (q2_lower_lim), (q2_upper_lim), 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_pim_evt[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_pim_evt/%s", name_evt));
            h_simu_pim_background_evt[w] = (THnSparseD *)root_mc_background->Get(Form("THnSparse_7D_pim_evt/%s", name_evt));
            h_simu_pim_evt[w]->Add(h_simu_pim_background_evt[w], -1.0);

            // // Write to output file

            h_simu_pim_evt[w]->Write();
        }
    }
    // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_pip_evt");
    outputFile->cd("THnSparse_7D_pip_evt");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (q2_lower_lim), (q2_upper_lim), 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_pip_evt[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_pip_evt/%s", name_evt));
            h_simu_pip_background_evt[w] = (THnSparseD *)root_mc_background->Get(Form("THnSparse_7D_pip_evt/%s", name_evt));
            h_simu_pip_evt[w]->Add(h_simu_pip_background_evt[w], -1.0);

            // Write to output file
            h_simu_pip_evt[w]->Write();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_thrown_prot_evt");
    outputFile->cd("THnSparse_7D_thrown_prot_evt");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (q2_lower_lim), (q2_upper_lim), 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_thrown_prot_evt[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_thrown_prot_evt/%s", name_evt));
            h_simu_thrown_prot_evt[w]->Write();
        }
    }
    // // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_thrown_pim_evt");
    outputFile->cd("THnSparse_7D_thrown_pim_evt");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (q2_lower_lim), (q2_upper_lim), 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_thrown_pim_evt[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_thrown_pim_evt/%s", name_evt));
            h_simu_thrown_pim_evt[w]->Write();
        }
    }
    // Move to the output file
    outputFile->cd();

    // Create the directory for histograms
    outputFile->mkdir("THnSparse_7D_thrown_pip_evt");
    outputFile->cd("THnSparse_7D_thrown_pip_evt");

    for (short q2 = 0; q2 < Q2_bins_no; q2++)
    {
        float q2_lower_lim = q2_low_values[q2];
        float q2_upper_lim = q2_up_values[q2];

        for (short w = 0; w < W_bins_no; w++)
        {
            auto name_evt = Form("h_5dim_evt_%.1f<=Q2<=%.1f GeV2_%.2f<=W<=%.2f GeV", (q2_lower_lim), (q2_upper_lim), 1.40 + 0.05 * w, 1.40 + 0.05 * w + 0.05);

            // Retrieve histograms from ROOT files
            h_simu_thrown_pip_evt[w] = (THnSparseD *)root_mc->Get(Form("THnSparse_7D_thrown_pip_evt/%s", name_evt));
            h_simu_thrown_pip_evt[w]->Write();
        }
    }

    // Close files
    outputFile->Close();
    root_mc->Close();
    root_mc_background->Close();

    std::cout << "Background-subtracted histograms saved successfully!" << std::endl;
    return 0;
}
