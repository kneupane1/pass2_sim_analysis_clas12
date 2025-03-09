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
                              const std::string &folderName, const std::string &histName, const std::string &name,
                              double x1, double x2, double y1, double y2, double el, double eh, double sl, double sh, double eff1, double eff2, double r1, double r2)
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
    // TCanvas *canvas = new TCanvas("canvas", "Efficiency Analysis", 1500, 400);
    // canvas->Divide(3, 1, 0.01, 0.01); // 3 rows, with spacing
    TCanvas *canvas = new TCanvas("canvas", "Efficiency Analysis", 1500, 1500);
    canvas->Divide(1, 1); // 3 rows, with spacing

    gStyle->SetOptStat(0);
    canvas->SetMargin(0.15, 0.2, 0.1, 0.1); // left, right, top, buttom

    // canvas->cd(1);
    // canvas->SetMargin(0.2, 0.05, 0.2, 0.1); // Wider margins
    // expHist->SetTitle((name + " Exp. data").c_str());
    // expHist->SetTitleSize(0.06, "t"); // Increase plot title size

    // expHist->GetXaxis()->SetRangeUser(x1, x2);
    // expHist->GetYaxis()->SetRangeUser(y1, y2);
    // expHist->GetXaxis()->SetTitleSize(0.05);
    // expHist->GetYaxis()->SetTitleSize(0.05);
    // expHist->GetXaxis()->SetLabelSize(0.04);
    // expHist->GetYaxis()->SetLabelSize(0.04);
    // expHist->GetXaxis()->SetTitleOffset(1.);
    // expHist->GetYaxis()->SetTitleOffset(1.);
    // expHist->SetMinimum(el);
    // expHist->SetMaximum(eh);
    // expHist->Draw("COLZ");
    // std::cout << "  r1 = " << r1 << "  r2 =  " << r2 << std::endl;
    // // Draw circular boundaries
    // TEllipse *circle1 = new TEllipse(0, 0, r1); // Center at (0, 0), radius r1
    // circle1->SetLineColor(kBlack);
    // circle1->SetLineStyle(2); // Dashed line
    // circle1->SetFillStyle(1); // Transparent
    // circle1->Draw("same");

    // TEllipse *circle2 = new TEllipse(0, 0, r2); // Center at (0, 0), radius r2
    // circle2->SetLineColor(kBlack);
    // circle2->SetLineStyle(2); // Dashed line
    // circle2->SetFillStyle(1); // Transparent
    // circle2->Draw("same");

    // // Draw circular boundaries
    // TEllipse *circle1 = new TEllipse(0, 0, r1 + 1); // Center at (0, 0), radius r1
    // circle1->SetLineColor(kRed);
    // circle1->SetLineStyle(2); // Dashed line
    // circle1->SetFillStyle(1); // Transparent
    // circle1->Draw("same");

    // TEllipse *circle2 = new TEllipse(0, 0, r2 - 1); // Center at (0, 0), radius r2
    // circle2->SetLineColor(kRed);
    // circle2->SetLineStyle(2); // Dashed line
    // circle2->SetFillStyle(1); // Transparent
    // circle2->Draw("same");
    // // Draw circular boundaries
    // TEllipse *circle1 = new TEllipse(0, 0, r1 - 1); // Center at (0, 0), radius r1
    // circle1->SetLineColor(kGreen);
    // circle1->SetLineStyle(2); // Dashed line
    // circle1->SetFillStyle(1); // Transparent
    // circle1->Draw("same");

    // TEllipse *circle2 = new TEllipse(0, 0, r2 + 1); // Center at (0, 0), radius r2
    // circle2->SetLineColor(kGreen);
    // circle2->SetLineStyle(2); // Dashed line
    // circle2->SetFillStyle(1); // Transparent
    // circle2->Draw("same");

    // // // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, -50)
    // // TLine *line = new TLine(280, -125, 335, -50);
    // // line->SetLineColor(kRed);
    // // line->SetLineStyle(1); // Solid line
    // // line->SetLineWidth(2); // Make it more visible
    // // line->Draw("same");

    // // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, 0)
    // TLine *line1 = new TLine(280, 125, 335, 50);
    // line1->SetLineColor(kRed);
    // line1->SetLineStyle(1); // Solid line
    // line1->SetLineWidth(2); // Make it more visible
    // line1->Draw("same");

    // // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, -50)
    // TLine *line2 = new TLine(330, -400, 330, 400);
    // line2->SetLineColor(kRed);
    // line2->SetLineStyle(1); // Solid line
    // line2->SetLineWidth(2); // Make it more visible
    // line2->Draw("same");

    // // // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, -50)
    // // TLine *line3 = new TLine(128, -400, 128, 400);
    // // line3->SetLineColor(kRed);
    // // line3->SetLineStyle(1); // Solid line
    // // line3->SetLineWidth(2); // Make it more visible
    // // line3->Draw("same");
    // // Define the parabola function
    // // TF1 *parabola = new TF1("parabola", "-0.21875*(x - 300)^2 - 37.5", 270, 330);
    // // TF1 *parabola = new TF1("parabola", " -0.01 * (x - 274) * (x - 274) - 75", 250, 400);

    // // parabola->SetLineColor(kBlue);
    // // parabola->SetLineWidth(2); // Thicker line for better visibility
    // // parabola->SetLineStyle(1); // Solid line

    // // Draw the parabola
    // // parabola->Draw("same");

    // // Plot simulation histogram
    // canvas->cd(2);
    // canvas->SetMargin(0.2, 0.2, 0.2, 0.1); // Extra space on the right for color axis
    // simHist->SetTitle((name + " MC data").c_str());
    // simHist->SetTitleSize(0.06, "t"); // Increase plot title size

    // simHist->GetXaxis()->SetRangeUser(x1, x2);
    // simHist->GetYaxis()->SetRangeUser(y1, y2);
    // simHist->GetXaxis()->SetTitleSize(0.05);
    // simHist->GetYaxis()->SetTitleSize(0.05);
    // simHist->GetXaxis()->SetLabelSize(0.04);
    // simHist->GetYaxis()->SetLabelSize(0.04);
    // simHist->GetXaxis()->SetTitleOffset(1.);
    // simHist->GetYaxis()->SetTitleOffset(1.);
    // simHist->SetMinimum(sl);
    // simHist->SetMaximum(sh);

    // // // Customize the color bar
    // // simHist->GetZaxis()->SetLabelSize(0.04);  // Increase color bar label size
    // // simHist->GetZaxis()->SetTitleSize(0.05);  // Increase color bar title size
    // simHist->GetZaxis()->SetTitleOffset(0.01); // Adjust color bar title position
    // simHist->Draw("COLZ");
    // // circle1->Draw("same");
    // // circle2->Draw("same");
    // // line->Draw("same");
    // // line1->Draw("same");
    // // line2->Draw("same");

    // // parabola->Draw("same");

    // // Plot efficiency histogram
    // canvas->cd(3);
    canvas->SetMargin(0.2, 0.05, 0.2, 0.1); // Similar margins
    setCustomPalette();
    effHist->SetTitle((name + " Efficiency").c_str());
    effHist->SetTitleSize(0.06, "t"); // Increase plot title size
    effHist->GetXaxis()->SetRangeUser(x1, x2);
    effHist->GetYaxis()->SetRangeUser(y1, y2);
    effHist->GetXaxis()->SetTitleSize(0.05);
    effHist->GetYaxis()->SetTitleSize(0.05);
    effHist->GetXaxis()->SetLabelSize(0.04);
    effHist->GetYaxis()->SetLabelSize(0.04);
    effHist->GetXaxis()->SetTitleOffset(1.);
    effHist->GetYaxis()->SetTitleOffset(1.);
    effHist->SetMinimum(eff1);
    effHist->SetMaximum(eff2);
    effHist->Draw("COLZ");

    // Define the X range where the lines should be drawn
    double xMin = 280;
    double xMax = 335;

    // Original Y-intercepts
    double intercept1 = -505.8;
    double intercept2 = 505.8;

    // Compute 10% outward and inward variations
    double delta1 = 5; // 10% of intercept1
    double delta2 = 5; // 10% of intercept2

    // Compute Y values for original lines
    double yMin1 = 1.36 * xMin + intercept1;
    double yMax1 = 1.36 * xMax + intercept1;
    double yMin2 = -1.36 * xMin + intercept2;
    double yMax2 = -1.36 * xMax + intercept2;

    // Compute Y values for outward shifted lines (+10%)
    double yMin1_tight = 1.36 * xMin + (-505.8 + delta1);
    double yMax1_tight = 1.36 * xMax + (-505.8 + delta1);

    double yMin2_loose = -1.36 * xMin + (505.8 + delta2);
    double yMax2_loose = -1.36 * xMax + (505.8 + delta2);

    // Compute Y values for inward shifted lines (-10%)
    double yMin1_loose = 1.36 * xMin + (-505.8 - delta1);
    double yMax1_loose = 1.36 * xMax + (-505.8 - delta1);

    double yMin2_tight = -1.36 * xMin + (505.8 - delta2);
    double yMax2_tight = -1.36 * xMax + (505.8 - delta2);

    // Draw original lines
    TLine *line1 = new TLine(xMin, yMin1, xMax, yMax1);
    line1->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line1->Draw("same");

    TLine *line2 = new TLine(xMin, yMin2, xMax, yMax2);
    line2->SetLineColor(kBlue);
    line2->SetLineWidth(2);
    line2->Draw("same");

    /////////////////////  lower cut ///////////////////
    // Draw outward shifted lines (+10%)
    TLine *line1_out = new TLine(xMin, yMin1_tight, xMax, yMax1_tight);
    line1_out->SetLineColor(kRed);
    line1_out->SetLineStyle(2); // Dashed line
    line1_out->SetLineWidth(2);
    line1_out->Draw("same");

    // Draw inward shifted lines (-10%)
    TLine *line1_in = new TLine(xMin, yMin1_loose, xMax, yMax1_loose);
    line1_in->SetLineColor(kBlue);
    line1_in->SetLineStyle(2); // Dashed line
    line1_in->SetLineWidth(2);
    line1_in->Draw("same");

    // ///////////////// upper cut
    TLine *line2_in = new TLine(xMin, yMin2_tight, xMax, yMax2_tight);
    line2_in->SetLineColor(kRed);
    line2_in->SetLineStyle(2); // Dashed line
    line2_in->SetLineWidth(2);
    line2_in->Draw("same");
    TLine *line2_out = new TLine(xMin, yMin2_loose, xMax, yMax2_loose);
    line2_out->SetLineColor(kBlue);
    line2_out->SetLineStyle(2); // Dashed line
    line2_out->SetLineWidth(2);
    line2_out->Draw("same");

    // // Define the X range where the lines should be drawn
    // double xMin = 260;
    // double xMax = 335;

    // // Compute corresponding Y values using the equations
    // double yMin1 = 1.36 * xMin - 505.8;
    // double yMax1 = 1.36 * xMax - 505.8;
    // double yMin2 = -1.36 * xMin + 505.8;
    // double yMax2 = -1.36 * xMax + 505.8;

    // // Draw the first line (Y3 = 1.36 * X3 - 505.8)
    // TLine *line1 = new TLine(xMin, yMin1, xMax, yMax1);
    // line1->SetLineColor(kRed);
    // line1->SetLineStyle(1); // Solid line
    // line1->SetLineWidth(2);
    // line1->Draw("same");

    // // Draw the second line (Y3 = -1.36 * X3 + 505.8)
    // TLine *line2 = new TLine(xMin, yMin2, xMax, yMax2);
    // line2->SetLineColor(kRed);
    // line2->SetLineStyle(1); // Solid line
    // line2->SetLineWidth(2);
    // line2->Draw("same");

    // // Define the X range where the lines should be drawn
    // double xMin1 = 250;
    // double xMax1 = 345;

    // // Compute corresponding Y values using the equations
    // double yMin11 = 1.36 * xMin - 505.8 * 1.1;
    // double yMax11 = 1.36 * xMax - 505.8 * 1.1;
    // double yMin21 = -1.36 * xMin + 505.8 * 1.1;
    // double yMax21 = -1.36 * xMax + 505.8 * 1.1;

    // // Draw the first line (Y3 = 1.36 * X3 - 505.8)
    // TLine *line3 = new TLine(xMin1, yMin11, xMax1, yMax11);
    // line3->SetLineColor(kBlue);
    // line3->SetLineStyle(1); // Solid line
    // line3->SetLineWidth(2);
    // line3->Draw("same");

    // // Draw the second line (Y3 = -1.36 * X3 + 505.8)
    // TLine *line4 = new TLine(xMin1, yMin2, xMax1, yMax2);
    // line4->SetLineColor(kBlue);
    // line4->SetLineStyle(1); // Solid line
    // line4->SetLineWidth(2);
    // line4->Draw("same");

    // // // circle1->Draw("same");
    // // circle2->Draw("same");
    // // line->Draw("same");
    // // line1->Draw("same");
    // // line2->Draw("same");
    // // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, -50)
    // TLine *line = new TLine(280, -125, 335, -50);
    // line->SetLineColor(kRed);
    // line->SetLineStyle(1); // Solid line
    // line->SetLineWidth(2); // Make it more visible
    // line->Draw("same");
    // // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, 0)
    // TLine *line1 = new TLine(280, 125, 335, 50);
    // line1->SetLineColor(kRed);
    // line1->SetLineStyle(1); // Solid line
    // line1->SetLineWidth(2); // Make it more visible
    // line1->Draw("same");

    // // Draw the line from (x1, y1) = (280, -125) to (x2, y2) = (335, -50)
    // TLine *line2 = new TLine(330, -400, 330, 400);
    // line2->SetLineColor(kRed);
    // line2->SetLineStyle(1); // Solid line
    // line2->SetLineWidth(2); // Make it more visible
    // line2->Draw("same");

    float a_a = 1;
    std ::cout << "  r1 = " << r1 << "  r2 =  " << r2 << std::endl;
    // Draw circular boundaries
    TEllipse *circle1 = new TEllipse(0, 0, r1); // Center at (0, 0), radius r1
    circle1->SetLineColor(kBlack);
    circle1->SetLineStyle(2); // Dashed line
    circle1->SetFillStyle(1); // Transparent
    circle1->Draw("same");

    TEllipse *circle2 = new TEllipse(0, 0, r2); // Center at (0, 0), radius r2
    circle2->SetLineColor(kBlack);
    circle2->SetLineStyle(2); // Dashed line
    circle2->SetFillStyle(1); // Transparent
    circle2->Draw("same");

    // Draw circular boundaries
    TEllipse *circle3 = new TEllipse(0, 0, r1 + a_a); // Center at (0, 0), radius r1
    circle3->SetLineColor(kRed);
    circle3->SetLineStyle(2); // Dashed line
    circle3->SetFillStyle(1); // Transparent
    circle3->Draw("same");

    TEllipse *circle4 = new TEllipse(0, 0, r2 - a_a); // Center at (0, 0), radius r2
    circle4->SetLineColor(kRed);
    circle4->SetLineStyle(2); // Dashed line
    circle4->SetFillStyle(1); // Transparent
    circle4->Draw("same");
    // Draw circular boundaries
    TEllipse *circle5 = new TEllipse(0, 0, r1 - a_a); // Center at (0, 0), radius r1
    circle5->SetLineColor(kGreen);
    circle5->SetLineStyle(2); // Dashed line
    circle5->SetFillStyle(1); // Transparent
    circle5->Draw("same");

    TEllipse *circle6 = new TEllipse(0, 0, r2 + a_a); // Center at (0, 0), radius r2
    circle6->SetLineColor(kGreen);
    circle6->SetLineStyle(2); // Dashed line
    circle6->SetFillStyle(1); // Transparent
    circle6->Draw("same");

    // parabola->Draw("same");

    // Save canvas
    canvas->SaveAs(("/Users/krishnaneupane/Documents/AA_PhD_useful_documents/My_PhD_work/detector_ineff_plots/tight_med_loose/" + name + "_test_Ineff_cuts.png").c_str());

    delete expFile;
    delete simFile;
    delete canvas;
}

// Main function to loop over histograms
struct HistogramInfo
{
    std::string histName;
    std::string name;
    double x1Limit;
    double x2Limit;
    double y1Limit;
    double y2Limit;

    double elLimit;
    double ehLimit;
    double slLimit;
    double shLimit;

    double efflLimit;
    double effhLimit;

    double r1;
    double r2;
};

void det_ineff_plots()
{
    // const std::string expFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/resIncl_EXP_pass2_126_runs_after_all_pid_cuts_for_det_ineff_with_cdfd_with_EB_ID_with_QADB_cuts.root";
    // const std::string simFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/pass2_741_files_for_det_ineff_cuts_with_cdfd_cuts_after_twoPi_event_after_mmsq.root";
    // ///// with circular inner and outer cuts
    // const std::string expFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/resIncl_EXP_pass2_126_runs_afteer_circular_cuts_for_det_mod_ineff_with_cdfd_with_EB_ID_with_QADB_cuts.root";
    // const std::string simFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/output_Pass2_sim_twoPi_rga_fall2018_tor-1_sol-1_flagrad_2_bg_45nA_job_85_all_60.root";

    // ///// with circular and tringular inner and outer cuts

    const std::string expFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/resIncl_EXP_pass2_126_runs_afteer_circular_cuts_for_det_mod_new_ineff_with_cdfd_with_EB_ID_with_QADB_cuts.root";
    const std::string simFilePath = "/Users/krishnaneupane/Downloads/2025/Det_ineff_cuts/output_Pass2_sim_twoPi_rga_fall2018_tor-1_sol-1_flagrad_2_bg_45nA_job_85_all_with_new_mod.root";

    // // // ////////////////////// Elec Cuts ///////////////////
    // // // // ////////////////////// Elec Cuts ///////////////////
    // const std::string folderName = "Electron_Cuts";
    // // // // // List of histograms to process
    // std::vector<HistogramInfo> histograms = {
    //     {"pcal_sec_after_all_cuts", "Electron PCAL y versus x", -350, 350, -350, 350, 1, 200, 0.005, 1.95, 0.05, 3, 80, 260},
    //     {"dcr1_sec_after_all_cuts", "Electron DCR1 y versus x", -150, 150, -150, 150, 1, 200, 0.005, 1.95, 0.05, 3, 38, 101},
    //     {"dcr2_sec_after_all_cuts", "Electron DCR2 y versus x", -220, 220, -220, 220, 1, 200, 0.005, 1.95, 0.05, 3, 54, 155},
    //     {"dcr3_sec_after_all_cuts", "Electron DCR3 y versus x", -300, 300, -300, 300, 1, 200, 0.005, 1.95, 0.05, 3, 65, 210},
    // };
    // //     double min_radious = 75;
    // // double max_radious = 270;
    // // {{101, 101, 101, 101, 101, 101},
    // //  {155, 155, 155, 155, 155, 155},
    // //  {210, 210, 210, 210, 210, 210}},

    ////////////////////// FD theta vs Mom ///////////////////
    ////////////////////// FD theta vs Mom ///////////////////
    ////////////////////// FD theta vs Mom ///////////////////

    // // const std::string folderName = "elec_theta_vs_mom_fd";
    // // const std::string folderName = "proton_theta_vs_mom_fd_sec";
    // const std::string folderName = "pip_th_vs_mom_fd_sec";

    // std::vector<HistogramInfo> histograms;
    // for (int sec = 1; sec <= 6; ++sec)
    // {
    //     // histograms.push_back({"Theta_fd_elec_lab_vs_mom_elec_sec_" + std::to_string(sec),
    //     //                       "electron_theta_vs_mom_fd_sec" + std::to_string(sec), 3, 10, 5, 30});
    //     // histograms.push_back({"Theta_fd_prot_lab_vs_mom_prot_sec_" + std::to_string(sec),
    //     //                       "proton_theta_vs_mom_fd_sec" + std::to_string(sec), 0, 10, 0, 50});
    //     histograms.push_back({"Theta_fd_pip_lab_vs_mom_pip_sec_" + std::to_string(sec),
    //                           "pip_theta_vs_mom_fd_sec" + std::to_string(sec), 0, 5, 0, 50});
    // }

    // //////////////////// FD prot ///////////////////
    // //////////////////// FD prot ///////////////////

    // const std::string folderName = "proton_cuts_fd";

    // std::vector<HistogramInfo> histograms = {
    //     // {"Theta_fd_prot_lab_vs_mom_prot_after_all_cuts", "Theta vs mom prot fd all", -350, 350, -350, 350},
    //     {"dcr1_sec_prot_after_all_cuts", "Proton DCR1 y versus x", -170, 170, -170, 170, 1, 50, 0.005, 0.1, 0.05, 3, 38, 145},
    //     {"dcr2_sec_prot_after_all_cuts", "Proton DCR2 y versus x", -270, 270, -270, 270, 1, 50, 0.005, 0.1, 0.05, 3, 62, 220},
    //     {"dcr3_sec_prot_after_all_cuts", "Proton DCR3 y versus x", -420, 420, -420, 420, 1, 55, 0.005, 0.1, 0.05, 3, 135, 350},
    // };
    // //////////////////// FD pip ///////////////////
    // //////////////////// FD pip ///////////////////
    //     upper radius for protons
    //     {{145, 145, 145, 145, 145, 145},
    //      {220, 220, 220, 220, 220, 220},
    //      {350, 350, 350, 350, 350, 350}},
    //     upper radius for pip
    //     {
    //         {140, 145, 145, 145, 145, 140},
    //             {210, 220, 220, 220, 220, 210},
    //             {325, 340, 340, 340, 340, 325}
    //     }
    // };
    const std::string folderName = "pip_cuts_fd";

    std::vector<HistogramInfo> histograms = {
        // {"Theta_fd_prot_lab_vs_mom_prot_after_all_cuts", "Theta vs mom prot fd all", -350, 350, -350, 350},
        {"dcr1_sec_pip_after_all_cuts", "Pip DCR1 y versus x", -170, 170, -170, 170, 1, 30, 0.005, 0.05, 0.0, 2.5, 36, 145},
        {"dcr2_sec_pip_after_all_cuts", "Pip DCR2 y versus x", -270, 270, -270, 270, 1, 30, 0.005, 0.05, 0.0, 2.5, 66, 220},
        {"dcr3_sec_pip_after_all_cuts", "Pip DCR3 y versus x", -420, 420, -420, 420, 1, 30, 0.005, 0.07, 0.0, 2.5, 135, 340},
    };

    // //////////////////// CD had ///////////////////
    // //////////////////// CD had ///////////////////
    // const std::string folderName = "proton_cuts_cd";

    // std::vector<HistogramInfo> histograms = {
    //     {"Theta_cd_prot_lab_vs_mom_prot_after_all_cuts", "Theta vs mom prot cd", 0, 4, 0, 60, 1, 80, 0.005, 2.55, 0.5, 5.5},
    // };

    // const std::string folderName = "pip_cuts_cd";

    // std::vector<HistogramInfo> histograms = {
    //     {"Theta_cd_pip_lab_vs_mom_pip_after_all_cuts", "Theta vs mom pip cd", 0, 4, 0, 150, 1, 80, 0.005, 2.55, 0.5, 5.5},
    // };

    // Process all histograms
    for (const auto &hist : histograms)
    {
        plotNormalizedHistograms(expFilePath, simFilePath, folderName, hist.histName, hist.name, hist.x1Limit, hist.x2Limit, hist.y1Limit, hist.y2Limit, hist.elLimit, hist.ehLimit, hist.slLimit, hist.shLimit, hist.efflLimit, hist.effhLimit, hist.r1, hist.r2);
    }
}
