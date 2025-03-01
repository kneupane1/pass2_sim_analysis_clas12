#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <vector>
#include <TEllipse.h>

// Set custom palette for better visualization
void setCustomPalette()
{
    const Int_t nRGBs = 10;
    Double_t stops[nRGBs] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    Double_t red[nRGBs] = {0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0};
    Double_t green[nRGBs] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.6, 0.4, 0.0};
    Double_t blue[nRGBs] = {1.0, 0.8, 0.6, 0.4, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0};

    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255);
    gStyle->SetNumberContours(255);
}

// Function to plot and normalize histograms
void plotNormalizedHistograms(const std::string &expFilePath, const std::string &simFilePath,
                              const std::string &folderName, const std::string &histName, const std::string &name)
{
    // Open experimental and simulation ROOT files
    TFile *expFile = TFile::Open(expFilePath.c_str(), "READ");
    TFile *simFile = TFile::Open(simFilePath.c_str(), "READ");

    // if (!expFile || expFile->IsZombie() || !simFile || simFile->IsZombie())
    // {
    //     std::cerr << "Error: Unable to open input files." << std::endl;
    //     return;
    // }

    // Retrieve histograms
    std::string expHistPath = folderName + "/" + histName;
    std::string simHistPath = folderName + "/" + histName;

    TH2D *expHist = dynamic_cast<TH2D *>(expFile->Get(expHistPath.c_str()));
    TH2D *simHist = dynamic_cast<TH2D *>(simFile->Get(simHistPath.c_str()));

    // if (!expHist || !simHist)
    // {
    //     std::cerr << "Error: Unable to find histograms: " << histName << std::endl;
    //     return;
    // }

    double totalExp = expHist->Integral(0, expHist->GetNbinsX() - 0, 0, expHist->GetNbinsY() - 0);
    double totalSim = simHist->Integral(0, simHist->GetNbinsX() - 0, 0, simHist->GetNbinsY() - 0);

    // Normalize histograms
    // double totalExp = expHist->Integral(0, 1);
    // double totalSim = simHist->Integral();
    double normFactor = totalExp / totalSim;

    // // double normFactor = (totalSim != 0) ? (totalExp / totalSim) : 1.0;
    // std::cout << "  normFactor  " << normFactor << std::endl;

    TH2D *effHist = dynamic_cast<TH2D *>(expHist->Clone("effHist"));
    effHist->Divide(simHist);
    effHist->Scale(1.0 / normFactor);

    // TH2D *effHist = dynamic_cast<TH2D *>(simHist->Clone("effHist"));
    // effHist->Divide(expHist);
    // effHist->Scale(1.0 * normFactor);

    // Style and plot
    TCanvas *canvas = new TCanvas("canvas", "Efficiency Analysis", 900, 400);
    canvas->Divide(2, 1, 0.01, 0.01); // 3 rows, with spacing

    gStyle->SetOptStat(0);
    canvas->SetMargin(0.15, 0.2, 0.1, 0.1); // left, right, top, buttom

    canvas->cd(1);
    canvas->SetMargin(0.5, 0.5, 0.2, 0.1); // Wider margins
    expHist->SetTitle((name + " Exp. data").c_str());
    expHist->SetTitleSize(0.06, "t"); // Increase plot title size

    // expHist->GetXaxis()->SetRangeUser(x1, x2);
    // expHist->GetYaxis()->SetRangeUser(y1, y2);
    expHist->GetXaxis()->SetTitleSize(0.05);
    expHist->GetYaxis()->SetTitleSize(0.05);
    expHist->GetXaxis()->SetLabelSize(0.04);
    expHist->GetYaxis()->SetLabelSize(0.04);
    expHist->GetXaxis()->SetTitleOffset(1.);
    expHist->GetYaxis()->SetTitleOffset(1.);
    // expHist->SetMinimum(el);
    // expHist->SetMaximum(eh);
    expHist->Draw("COLZ");

    // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, -50)
    TLine *line = new TLine(1.4, 2, 1.4, 9);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); // Solid line
    line->SetLineWidth(2); // Make it more visible
    line->Draw("same");

    // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, 0)
    TLine *line1 = new TLine(2.15, 2, 2.15, 9);
    line1->SetLineColor(kRed);
    line1->SetLineStyle(1); // Solid line
    line1->SetLineWidth(2); // Make it more visible
    line1->Draw("same");

    // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, -50)
    TLine *line2 = new TLine(1.4, 2.0, 2.15, 2.0);
    line2->SetLineColor(kRed);
    line2->SetLineStyle(1); // Solid line
    line2->SetLineWidth(2); // Make it more visible
    line2->Draw("same");

    // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, -50)
    TLine *line3 = new TLine(1.4, 9, 2.15, 9);
    line3->SetLineColor(kRed);
    line3->SetLineStyle(1); // Solid line
    line3->SetLineWidth(2); // Make it more visible
    line3->Draw("same");
    // Enable grid lines
    canvas->SetGridx();
    canvas->SetGridy();

    // Customize x-axis grid every 0.05
    gPad->SetGridx();
    // Customize y-axis grid at specific points
    double xGrid[] = {1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15};
    int nGridLinesx = sizeof(xGrid) / sizeof(xGrid[0]);

    // Customize y-axis grid at specific points
    double yGrid[] = {2.0, 2.4, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0, 9.0};
    int nGridLines = sizeof(yGrid) / sizeof(yGrid[0]);

    // Draw horizontal grid lines
    for (int i = 0; i < nGridLinesx; i++)
    {
        TLine *gridLine = new TLine(xGrid[i], 1.0, xGrid[i], 10.0);
        gridLine->SetLineColor(kGray + 2); // Light gray
        gridLine->SetLineStyle(2);         // Dashed
        gridLine->Draw("same");
    }
    // Draw horizontal grid lines
    for (int i = 0; i < nGridLines; i++)
    {
        TLine *gridLine = new TLine(1., yGrid[i], 2.5, yGrid[i]);
        gridLine->SetLineColor(kGray + 2); // Light gray
        gridLine->SetLineStyle(2);         // Dashed
        gridLine->Draw("same");
    }
    // Plot simulation histogram
    canvas->cd(2);
    canvas->SetMargin(0.5, 0.2, 0.2, 0.1); // Extra space on the right for color axis
    simHist->SetTitle((name + " MC data").c_str());
    simHist->SetTitleSize(0.06, "t"); // Increase plot title size

    simHist->GetXaxis()->SetTitleSize(0.05);
    simHist->GetYaxis()->SetTitleSize(0.05);
    simHist->GetXaxis()->SetLabelSize(0.04);
    simHist->GetYaxis()->SetLabelSize(0.04);
    simHist->GetXaxis()->SetTitleOffset(1.);
    simHist->GetYaxis()->SetTitleOffset(1.);

    // // Customize the color bar
    // simHist->GetZaxis()->SetLabelSize(0.04);  // Increase color bar label size
    // simHist->GetZaxis()->SetTitleSize(0.05);  // Increase color bar title size
    simHist->GetZaxis()->SetTitleOffset(0.01); // Adjust color bar title position
    simHist->Draw("COLZ");

    line->Draw("same");
    line1->Draw("same");
    line2->Draw("same");
    line3->Draw("same");

    // Enable grid lines
    canvas->SetGridx();
    canvas->SetGridy();
    gPad->SetGridx();

    // Draw horizontal grid lines
    for (int i = 0; i < nGridLinesx; i++)
    {
        TLine *gridLine = new TLine(xGrid[i], 1.0, xGrid[i], 10.0);
        gridLine->SetLineColor(kGray + 2); // Light gray
        gridLine->SetLineStyle(2);         // Dashed
        gridLine->Draw("same");
    }
    // Draw horizontal grid lines
    for (int i = 0; i < nGridLines; i++)
    {
        TLine *gridLine = new TLine(1., yGrid[i], 2.5, yGrid[i]);
        gridLine->SetLineColor(kGray + 2); // Light gray
        gridLine->SetLineStyle(2);         // Dashed
        gridLine->Draw("same");
    }
    // Save canvas
    canvas->SaveAs(("/Users/krishnaneupane/Documents/AA_PhD_useful_documents/My_PhD_work/" + name + "w_vs_q2.png").c_str());

    delete expFile;
    delete simFile;
    delete canvas;
}

// Main function to loop over histograms
struct HistogramInfo
{
    std::string histName;
    std::string name;
};

void pid_plots()
{
    // const std::string expFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/resIncl_EXP_pass2_126_runs_after_all_pid_cuts_for_det_ineff_with_cdfd_with_EB_ID_with_QADB_cuts.root";
    const std::string simFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/pass2_741_files_for_det_ineff_cuts_with_cdfd_cuts_after_twoPi_event_after_mmsq.root";
    // ///// with circular inner and outer cuts
    // const std::string expFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/resIncl_EXP_pass2_126_runs_afteer_circular_cuts_for_det_mod_ineff_with_cdfd_with_EB_ID_with_QADB_cuts.root";
    // const std::string simFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/output_Pass2_sim_twoPi_rga_fall2018_tor-1_sol-1_flagrad_2_bg_45nA_job_85_all_60.root";

    // ///// with circular and tringular inner and outer cuts
    // const std::string expFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/resIncl_EXP_pass2_126_runs_afteer_circular_cuts_for_det_mod_new_ineff_with_cdfd_with_EB_ID_with_QADB_cuts.root";

    const std::string expFilePath = "/Users/krishnaneupane/Downloads/2025/Feb_2025/resIncl_EXP_pass2_all_126_runs_for_w_q2_hist_with_EB_ID_with_QADB_cuts.root";
    // const std::string simFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/output_Pass2_sim_twoPi_rga_fall2018_tor-1_sol-1_flagrad_2_bg_45nA_job_85_all_with_new_mod.root";

    // // // ////////////////////// Elec Cuts ///////////////////
    // // // ////////////////////// Elec Cuts ///////////////////
    const std::string folderName = "W vs Q2";
    // // // // List of histograms to process

    std::vector<HistogramInfo> histograms = {
        {"W_vs_q2", "W versus Q^{2} Distribution"},

    };

    for (const auto &hist : histograms)
    {
        plotNormalizedHistograms(expFilePath, simFilePath, folderName, hist.histName, hist.name);

        //    plotNormalizedHistograms(expFilePath, simFilePath, folderName, hist.histName, hist.name, hist.x1Limit, hist.x2Limit, hist.y1Limit, hist.y2Limit, hist.elLimit, hist.ehLimit, hist.slLimit, hist.shLimit, hist.efflLimit, hist.effhLimit, hist.r1, hist.r2);
    }
}
