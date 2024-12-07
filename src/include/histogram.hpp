
#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TThread.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"
#include <mutex>

using namespace std;

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;
using THnSparse_ptr = std::shared_ptr<THnSparse>;
using TGraph_ptr = std::shared_ptr<TGraph>;

class Histogram
{
protected:
    std::shared_ptr<TFile> RootOutputFile;
    std::shared_ptr<TCanvas> def;

    int bins = 500;
    double p_min = 0.0;
    double p_max = 6.0;
    double Dt_max = 10.0;
    double Dt_min = -Dt_max;
    double q2_min = 0.0;
    double q2_max = 10.0;

    double w_max = 3.5;
    double w_min = 1.0;

    double zero = 0.0;

    static const short particle_num = 3; // 0-e 1-Pi 2-P 3-K
    std::string particle_name[particle_num] = {"e", "pi", "P"};
    static const short charge_num = 2; // 0-pos 1-neg
    std::string charge_name[charge_num] = {"positive", "negative"};
    static const short with_id_num = 3; // 0-without 1-with 2-anti
    std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};

    static const short num_sectors = 6;
    std::string sec_name[num_sectors] = {"1", "2", "3", "4", "5", "6"};

    // std::string cut_name[2] = {"bc", "ac"};

    static const short CUTS = 4;
    enum cuts
    {
        before_any_cuts,
        with_one_cut,
        outside_one_cut,
        after_all_cuts
    };
    std::mutex mutex;

    static const short w_range_num = 3;
    std::string w_range_name[w_range_num] = {
        "all_W_range", " W < 2.5 ",
        " W > 2.5 "}; //{" W<2.0 ", " 2.0<W<2.5 ", " 2.5<W<3.0 ", " 3.0<W<3.5 "};
    static const short q2_range_num = 3;
    std::string q2_range_name[q2_range_num] = /*{" Q2<1.0 ",     " 1.0<Q2<2.0 ",
                                               " 2.0<Q2<3.0 ", " 3.0<Q2<4.0 ",
                                               " 4.0<Q2<5.0 ", " 5.0<Q2<6.0 ",
                                               " 6.0<Q2<7.0 ", " 7.0<Q2<8.0 ",
                                               " 8.0<Q2<9.0 ", " Q2>9.0 "};*/
        {"all_Q2_range ", " Q2 < 4.5 ", " Q2 > 4.5 "};

    static const short inv_Ppip_range_num = 2;
    std::string inv_Ppip_range_name[inv_Ppip_range_num] = {" M[Ppip]<1.5 GeV ",
                                                           " M[Ppip]>1.5 GeV "};

    static const short inv_pip_pim_range_num = 2;
    std::string inv_pip_pim_range_name[inv_pip_pim_range_num] = {
        " M[pi+pi-]<0.9 GeV ", " M[pi+pi-]>0.9 GeV "};

    static const short theta_pim_range_num = 2;
    std::string theta_pim_range_name[theta_pim_range_num] = {
        " theta[pi-]<90 deg ", " theta[pi-]>90 deg "};

    static const short phi_pim_range_num = 2;
    std::string phi_pim_range_name[phi_pim_range_num] = {" phi[pi-]<90 deg ",
                                                         " phi[pi-]>90 deg "};
    static const short alpha_pim_range_num = 2;
    std::string alpha_pim_range_name[alpha_pim_range_num] = {
        " alpha[pi-]<180 deg ", " alpha[pi-]>180 deg "};

    static const short NUM_CONDITIONS = 1;
    std::string NUM_CONDITIONS_NAME[NUM_CONDITIONS] = {
        /*"twoPi_event" ,*/ "missingPim events"};
    // Kinematics

    static const short q2_bin = 11;
    float q2_low_values[10] = {1.0, 2.0, 2.40, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0};
    float q2_up_values[10] = {2.0, 2.40, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0, 9.0};
    int q2_bin_size = 10;
    int w_lower_bin = 8;
    int w_higher_bin = 24;

    // ////////////// backgraound multiplication factors obtained from AND logic in the background fitting using exclusive topology data /////////
    // float background_fact[16][7] = {
    //     {2.3, 2.5, 2.47, 2.5, 3.0, 2.09, 2.61},
    //     {2.49, 2.5, 2.5, 2.26, 2.93, 2.37, 2.54},
    //     {2.5, 2.46, 2.5, 2.5, 2.46, 3.17, 2.25},
    //     {2.49, 2.5, 2.47, 2.41, 2.46, 2.72, 3.18},
    //     {2.48, 2.49, 2.47, 2.46, 3.44, 3.41, 3.33},
    //     {2.38, 2.48, 2.49, 2.5, 3.34, 3.36, 4.05},
    //     {2.49, 2.5, 2.5, 2.5, 3.21, 3.09, 3.32},
    //     {3.19, 2.5, 2.5, 2.48, 2.95, 3.18, 2.83},
    //     {3.2, 3.13, 3.19, 3.19, 3.38, 3.04, 3.32},
    //     {3.19, 3.18, 3.2, 3.2, 3.22, 3.2, 3.37},
    //     {3.2, 3.2, 3.19, 3.2, 3.25, 2.95, 3.25},
    //     {3.15, 3.2, 3.2, 3.18, 3.3, 3.12, 3.57},
    //     {3.19, 3.2, 3.2, 3.02, 3.24, 3.24, 3.51},
    //     {3.16, 3.2, 3.18, 3.2, 3.39, 3.47, 3.55},
    //     {3.19, 3.2, 3.19, 3.2, 3.49, 3.42, 3.66}};
    // // // //sim FD+CD

    int bin_val = -1;

    // Define q2 bin ranges
    std::vector<float> q2_bins = {2.0, 2.4, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0, 9.0};

    // Find the bin based on q2

    int q2_bining(float q2)
    {

        for (int i = 1; i < q2_bins.size(); ++i)
        {
            if (q2 < q2_bins[i])
            {
                bin_val = i;
                break;
            }
        }
        return bin_val;
    }

    // /// MMSQ cuts
    int q2_bin_val = -1;
    /////////////////////////// exp, sim data mmsq cuts [Q2][up/down][a,b,c]   /////////////////  updated nov 2024 with 0.65 sm * w dependent sm
    double mmsq_cuts[2][9][2][3] = {
        {{{-0.1330, 0.5196, -0.4026},
          {0.0719, -0.3014, 0.2581}},
         {{-0.1176, 0.4645, -0.3538},
          {0.0410, -0.1844, 0.1470}},
         {{-0.1461, 0.5680, -0.4461},
          {0.0863, -0.3422, 0.2793}},
         {{-0.1028, 0.4131, -0.3115},
          {0.0642, -0.2602, 0.2022}},
         {{-0.0591, 0.2539, -0.1706},
          {0.0615, -0.2450, 0.1804}},
         {{-0.0358, 0.1777, -0.1085},
          {0.0511, -0.2008, 0.1339}},
         {{-0.0256, 0.1207, -0.0358},
          {0.0627, -0.2115, 0.1135}},
         {{-0.0500, 0.2273, -0.1382},
          {-0.0289, 0.1211, -0.1872}},
         {{-0.0822, 0.3223, -0.2004},
          {0.1153, -0.3838, 0.2522}}},
        {{{-0.1648, 0.6717, -0.5615},
          {0.0845, -0.3438, 0.2921}},
         {{-0.1332, 0.5590, -0.4635},
          {0.0880, -0.3502, 0.2887}},
         {{-0.1135, 0.4879, -0.3995},
          {0.0647, -0.2682, 0.2166}},
         {{-0.1316, 0.5516, -0.4541},
          {0.0614, -0.2510, 0.1951}},
         {{-0.1120, 0.4788, -0.3874},
          {0.0496, -0.1989, 0.1396}},
         {{-0.0901, 0.3937, -0.3071},
          {0.0307, -0.1256, 0.0684}},
         {{-0.0512, 0.2545, -0.1847},
          {0.0150, -0.0601, 0.0004}},
         {{-0.0414, 0.2194, -0.1526},
          {0.0160, -0.0525, -0.0182}},
         {{-0.0117, 0.1203, -0.0703},
          {-0.0040, 0.0208, -0.0896}}}};

    //////////////////////////////////////////////////
    int inv_mass_binning(float inv_mass, float inv_pPip_llim, float bin_size_inv)
    {
        for (int i = 0; i < 11; ++i)
        {
            float lower_limit = inv_pPip_llim + i * bin_size_inv;
            float upper_limit = inv_pPip_llim + (i + 1) * bin_size_inv;

            if (lower_limit <= inv_mass && inv_mass < upper_limit)
            {
                return i;
            }
        }
        return -1; // Return -1 if the mass is outside of the bins
    }
    /////////////////////////////////////////////////

    int inv_binning(float w, float inv_mass, bool hp)
    {
        // Find the index of the bin that w falls into
        int w_bin_index = static_cast<int>((w - 1.4) / 0.05);

        // Calculate the midpoint of the bin
        float w_mid = 1.4 + (w_bin_index + 0.5) * 0.05;
        float inv_ulim;
        float inv_llim;
        // Double_t xmin_5D[ndims_5D] = {((0.938272 + 0.13957) - 2 * Bin_size_pPip), (0.13957 + 0.13957) - 2 * Bin_size_pipPim, 0., 0., 0.};
        // Double_t xmax_5D[ndims_5D] = {((1.0 + 0.05 * w + 0.025 - MASS_PIM) + 2 * Bin_size_pPip), ((1.0 + 0.05 * w + 0.025 - MASS_P) + 2 * Bin_size_pipPim), 180, 360, 360};

        if (hp)
        {
            inv_ulim = w_mid - MASS_PIM;
            inv_llim = MASS_P + MASS_PIP;
            // std::cout << " w  " << w << "  w_bin_index  " << w_bin_index << "   w mid  " << w_mid << " inv_llim  " << inv_llim << "   inv_ulim  " << inv_ulim << std::endl;
        }
        else
        {
            inv_ulim = w_mid - MASS_P;
            inv_llim = MASS_PIP + MASS_PIM;
        }
        int inv_bin_val = -1;

        for (int j = 0; j < 7; ++j)
        {
            float lim_value = inv_llim + (j + 1) * (inv_ulim - inv_llim) / 7.0;
            if (inv_mass <= lim_value)
            {
                inv_bin_val = j;
                break;
            }
        }
        return inv_bin_val;
    }

    ////////////////////////////////////////////
    int alpha_bin_val = -1;

    // Define q2 bin ranges
    std::vector<float> alpha_bins = {36.0, 72.0, 108.0, 144.0, 180.0, 216.0, 252.0, 288.0, 324.0, 360.0};

    // Find the bin based on q2

    int alpha_bining(float alpha)
    {

        for (int i = 0; i < alpha_bins.size(); ++i)
        {
            if (alpha < alpha_bins[i])
            {
                alpha_bin_val = i;
                break;
            }
        }
        return alpha_bin_val;
    }

    // double dCosTh(float theta)
    // {
    //     double theta_res = 18.0;
    //     int theta_bin_index = static_cast<int>((theta) / theta_res);
    //     double lower_edge = theta_res * theta_bin_index;
    //     double upper_edge = lower_edge + theta_res;
    //     // Calculate the cosine values of the upper and lower bin edges
    //     double cos_lower = TMath::Cos(lower_edge * TMath::Pi() / 180.0);
    //     double cos_upper = TMath::Cos(upper_edge * TMath::Pi() / 180.0);

    //     double difference = std::abs(cos_upper - cos_lower);

    //     return difference;
    // }
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////

    static const short w_bin = 64;
    THnSparse *threeDHist[q2_bin][w_bin];
    THnSparse *sevenDHist_pim[q2_bin][w_bin];
    THnSparse *sevenD_Hist_thrown_pim[q2_bin][w_bin];
    THnSparse *sevenDHist_pip[q2_bin][w_bin];
    THnSparse *sevenD_Hist_thrown_pip[q2_bin][w_bin];
    THnSparse *sevenDHist_prot[q2_bin][w_bin];
    THnSparse *sevenD_Hist_thrown_prot[q2_bin][w_bin];

    THnSparse *h_5dim_prot_evt[q2_bin][w_bin];
    THnSparse *h_5dim_pip_evt[q2_bin][w_bin];
    THnSparse *h_5dim_pim_evt[q2_bin][w_bin];

    THnSparse *h_5dim_thrown_prot_evt[q2_bin][w_bin];
    THnSparse *h_5dim_thrown_pip_evt[q2_bin][w_bin];
    THnSparse *h_5dim_thrown_pim_evt[q2_bin][w_bin];

    TH1D_ptr w_gen_hist[q2_bin][w_bin];
    TH1D_ptr q2_gen_hist[q2_bin][w_bin];

    TH1D_ptr w_gen_hist_inv_pPip[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_inv_pPip[q2_bin][w_bin][11];

    TH1D_ptr w_gen_hist_inv_pPim[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_inv_pPim[q2_bin][w_bin][11];

    TH1D_ptr w_gen_hist_inv_pipPim[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_inv_pipPim[q2_bin][w_bin][11];

    TH1D_ptr inv_pPip_hist[q2_bin][w_bin][11];
    TH1D_ptr inv_pPim_hist[q2_bin][w_bin][11];
    TH1D_ptr inv_pipPim_hist[q2_bin][w_bin][11];
    // TH1D *histogram = new TH1D("histogram", "Title", 100, 0, 100);

    TH1D_ptr w_gen_hist_th_prot[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_th_prot[q2_bin][w_bin][11];
    TH1D_ptr w_gen_hist_th_pip[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_th_pip[q2_bin][w_bin][11];
    TH1D_ptr w_gen_hist_th_pim[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_th_pim[q2_bin][w_bin][11];

    TH1D_ptr prot_theta_hist[q2_bin][w_bin][11];
    TH1D_ptr pip_theta_hist[q2_bin][w_bin][11];
    TH1D_ptr pim_theta_hist[q2_bin][w_bin][11];

    TH1D_ptr w_gen_hist_al_prot[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_al_prot[q2_bin][w_bin][11];
    TH1D_ptr w_gen_hist_al_pip[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_al_pip[q2_bin][w_bin][11];
    TH1D_ptr w_gen_hist_al_pim[q2_bin][w_bin][11];
    TH1D_ptr q2_gen_hist_al_pim[q2_bin][w_bin][11];

    TH1D_ptr prot_alpha_hist[q2_bin][w_bin][11];
    TH1D_ptr pip_alpha_hist[q2_bin][w_bin][11];
    TH1D_ptr pim_alpha_hist[q2_bin][w_bin][11];

    static const short W_BIN_CHECK_NUM = 11;

    // std::string W_BIN_CHECK_NAME[W_BIN_CHECK_NUM] = {" All_W_range "," <1.30W<1.35 ",     " 1.35<W<1.40 ",
    //                                            " 1.40<W<1.45 ", " 1.45<W<1.50 ",
    //                                            " 1.50<W<1.55 ", " 1.55<W<1.60 ",
    //                                            " 1.60<W<1.65 ", " 1.65<W<1.70 ",
    //                                            " 1.70<W<1.75 ", " 1.75<W<1.80 "};

    // for inclusive
    std::string W_BIN_CHECK_NAME[W_BIN_CHECK_NUM] = {" ALL Q2 ", " Q2= 2.37 GeV2 ", " Q2= 2.77 GeV2 ", " Q2= 3.24 GeV2 ",
                                                     " Q2= 3.79 GeV2 ", " Q2= 4.44 GeV2 ",
                                                     " Q2= 5.19 GeV2 ", " Q2= 6.07 GeV2 ",
                                                     " Q2= 7.09 GeV2 ", " Q2= 8.29 GeV2 ",
                                                     " Q2= 9.70 GeV2 "};
    TH1D_ptr mm2_mPim_hist_check_th[11];
    TH1D_ptr W_hist_check_th[11];
    TH2D_ptr W_VS_Q2_hist_check_th[11];

    TH1D_ptr mm2_mPim_hist_check[11];
    TH1D_ptr W_hist_check[11];
    TH2D_ptr W_VS_Q2_hist_check[11];

    static const short NUM_CUT = 2;

    static const short theta_bin_NUM = 18;
    TH1D_ptr theta_pim_measured_3_sigma[theta_bin_NUM];
    std::string theta_bin_NAME[theta_bin_NUM] = {
        "0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100", "100-110", "110-120", "120-130", "130-140", "140-150", "150-160", "160-170", "170-180"};

    //  TH2D_ptr sf_hist = std::make_shared<TH2D>("SF", "SF", 500, 0, 10.5, 500,
    //  0, 0.5);
    //
    TH1D_ptr missing_Energy_hist = std::make_shared<TH1D>("e#pi^{+}#pi^{-}pX", "missing Energy", 500, -0.6, 0.6);

    TH1D_ptr diff_E_P_x_mu_hist_ = std::make_shared<TH1D>(
        "diff_E_P_x_mu ", "diff_E_P_x_mu ", 500, -1.0, 1.0);
    TH1D_ptr P_x_mu = std::make_shared<TH1D>("Mom P ", "Mom P ", 500, -2.0, 10.0);

    TH2D_ptr theta_vs_mom_elec[num_sectors];
    TH2D_ptr theta_vs_mom_prot[num_sectors];
    TH2D_ptr theta_vs_mom_pip[num_sectors];
    TH2D_ptr theta_vs_mom_pim[num_sectors];
    TH2D_ptr pim_phi_vs_theta_rec_FD_sec[num_sectors];
    TH2D_ptr pim_phi_vs_theta_rec_FD_after_exclusive_sec[num_sectors];
    TH2D_ptr pim_phi_vs_theta_measured_FD_sec[num_sectors];
    TH2D_ptr pim_phi_vs_theta_measured_FD_after_exclusive_sec[num_sectors];

    // missingPiP
    static const short HADRON_NUM = 3;
    std::string HADRON_NAME[HADRON_NUM] = {
        "mProt", "mPip", "mPim"};
    static const short EFF_CONDITIONS_NUM_ALL = 3;
    std::string EFF_CONDITIONS_NAME_ALL[EFF_CONDITIONS_NUM_ALL] = {
        "missing", "exclusive", "after_MMSQ_exclusive_cuts"};
    TH2D_ptr pip_theta_rec_vs_mom[EFF_CONDITIONS_NUM_ALL];
    TH2D_ptr prot_theta_rec_vs_mom[EFF_CONDITIONS_NUM_ALL];
    TH2D_ptr pip_theta_measured_vs_mom[EFF_CONDITIONS_NUM_ALL];
    TH2D_ptr prot_theta_measured_vs_mom[EFF_CONDITIONS_NUM_ALL];

    TH1D_ptr weight_hist;
    TH1D_ptr pid_at_zero;
    TH1D_ptr mc_pid_at_zero;
    TH1D_ptr inv_mass_pPip;
    TH1D_ptr inv_mass_pPim;
    TH1D_ptr inv_mass_pipPim;
    TH1D_ptr theta_Prot_cm;
    TH1D_ptr theta_Pip_cm;
    TH1D_ptr theta_Pim_cm;
    TH1D_ptr phi_Prot_cm;
    TH1D_ptr phi_Pip_cm;
    TH1D_ptr phi_Pim_cm;
    TH1D_ptr alpha_Prot_cm;
    TH1D_ptr alpha_Pip_cm;
    TH1D_ptr alpha_Pim_cm;

    TH1D_ptr inv_mass_pPip_swapped;
    TH1D_ptr inv_mass_pPim_swapped;
    TH1D_ptr inv_mass_pipPim_swapped;
    TH1D_ptr theta_Prot_cm_swapped;
    TH1D_ptr theta_Pip_cm_swapped;
    TH1D_ptr theta_Pim_cm_swapped;
    TH1D_ptr phi_Prot_cm_swapped;
    TH1D_ptr phi_Pip_cm_swapped;
    TH1D_ptr phi_Pim_cm_swapped;
    TH1D_ptr alpha_Prot_cm_swapped;
    TH1D_ptr alpha_Pip_cm_swapped;
    TH1D_ptr alpha_Pim_cm_swapped;

    TH1D_ptr dp_prot_hist;
    TH1D_ptr dp_pip_hist;
    TH1D_ptr dp_ambi_prot_all_hist;
    TH1D_ptr dp_ambi_pip_all_hist;
    TH1D_ptr dp_ambi_prot_hist;
    TH1D_ptr dp_ambi_pip_hist;
    TH1D_ptr dp_prot_for_pip_hist;
    TH1D_ptr dp_pip_for_prot_hist;

    TH1D_ptr entries_in_each_event;
    TH1D_ptr entries_prot;
    TH1D_ptr entries_pip;
    TH1D_ptr MM2_mPim_all_comb;
    TH1D_ptr MM2_mPim_1_comb;
    TH1D_ptr MM2_mPim_2_comb;
    TH1D_ptr MM2_mPim_3_comb;
    TH1D_ptr MM2_mPim_4_or_more_comb;

    TH1D_ptr p_gen_prot_hist;
    TH1D_ptr p_gen_pip_hist;
    TH1D_ptr p_gen_ambi_prot_hist;
    TH1D_ptr p_gen_ambi_pip_hist;
    TH1D_ptr p_gen_prot_for_pip_hist;
    TH1D_ptr p_gen_pip_for_prot_hist;

    TH1D_ptr p_rec_prot_hist;
    TH1D_ptr p_rec_pip_hist;
    TH1D_ptr p_rec_ambi_prot_hist;
    TH1D_ptr p_rec_ambi_pip_hist;
    TH1D_ptr p_rec_prot_for_pip_hist;
    TH1D_ptr p_rec_pip_for_prot_hist;

    TH1D_ptr scalar_triple_product_hist;
    TH1D_ptr scalar_triple_product_hist_with_mmsq_cuts;
    TH2D_ptr W_vs_sclar_product;
    TH2D_ptr W_vs_scalar_product_after_mmsq_cuts;

    TH1D_ptr W_hist;
    TH1D_ptr Q2_hist;
    TH2D_ptr W_vs_q2;
    TH1D_ptr W_P2pi_hist;

    TH1D_ptr W_thrown;
    TH2D_ptr W_vs_Q2_thrown;
    TH1D_ptr Q2_thrown;

    TH1D_ptr vz_position[CUTS];
    TH2D_ptr pcal_sec[CUTS];
    TH2D_ptr pcal_hx_hy_sec[CUTS];
    TH2D_ptr dcr1_sec[CUTS];
    TH2D_ptr dcr2_sec[CUTS];
    TH2D_ptr dcr3_sec[CUTS];

    // EC Sampling Fraction
    TH2D_ptr EC_sampling_fraction[CUTS];
    TH2D_ptr ECin_sf_vs_PCAL_sf[CUTS];
    TH1D_ptr momentum[CUTS];

    //// Hadron pid
    TH1D_ptr prot_Delta_vz_cut_fd[CUTS];
    TH1D_ptr prot_Chi2pid_cut_fd[CUTS];
    TH1D_ptr pip_Delta_vz_cut_fd[CUTS];
    TH1D_ptr pip_Chi2pid_cut_fd[CUTS];
    TH1D_ptr pim_Delta_vz_cut[CUTS];
    TH1D_ptr pim_Chi2pid_cut[CUTS];

    TH1D_ptr prot_Delta_vz_cut_cd[CUTS];
    TH1D_ptr prot_Chi2pid_cut_cd[CUTS];
    TH1D_ptr pip_Delta_vz_cut_cd[CUTS];
    TH1D_ptr pip_Chi2pid_cut_cd[CUTS];

    TH2D_ptr dcr1_sec_prot[CUTS];
    TH2D_ptr dcr2_sec_prot[CUTS];
    TH2D_ptr dcr3_sec_prot[CUTS];
    TH2D_ptr dcr1_sec_pip[CUTS];
    TH2D_ptr dcr2_sec_pip[CUTS];
    TH2D_ptr dcr3_sec_pip[CUTS];
    TH2D_ptr dcr1_sec_pim[CUTS];
    TH2D_ptr dcr2_sec_pim[CUTS];
    TH2D_ptr dcr3_sec_pim[CUTS];
    TH2D_ptr phi_vs_mom_prot_fd[CUTS];
    TH2D_ptr phi_vs_mom_pip_fd[CUTS];
    TH2D_ptr phi_vs_mom_pim_fd[CUTS];
    TH1D_ptr theta_prot_fd[CUTS];
    TH1D_ptr theta_pip_fd[CUTS];
    TH1D_ptr theta_pim_fd[CUTS];
    TH2D_ptr Theta_prot_lab_vs_mom_prot_fd[CUTS];
    TH2D_ptr Theta_pip_lab_vs_mom_pip_fd[CUTS];
    TH2D_ptr Theta_pim_lab_vs_mom_pim_fd[CUTS];

    TH2D_ptr phi_vs_momT_prot_cd[CUTS];
    TH2D_ptr phi_vs_momT_pip_cd[CUTS];
    TH2D_ptr phi_vs_momT_pim_cd[CUTS];

    TH1D_ptr theta_prot_cd[CUTS];
    TH1D_ptr theta_pip_cd[CUTS];
    TH1D_ptr theta_pim_cd[CUTS];
    TH2D_ptr Theta_prot_lab_vs_mom_prot_cd[CUTS];
    TH2D_ptr Theta_pip_lab_vs_mom_pip_cd[CUTS];
    TH2D_ptr Theta_pim_lab_vs_mom_pim_cd[CUTS];

    TH2D_ptr Theta_prot_cm_vs_mom_prot;
    TH2D_ptr Theta_pip_cm_vs_mom_pip;
    TH2D_ptr Theta_pim_cm_vs_mom_pim;

    TH2D_ptr Theta_prot_thrown_cm_vs_mom_prot;
    TH2D_ptr Theta_pip_thrown_cm_vs_mom_pip;
    TH2D_ptr Theta_pim_thrown_cm_vs_mom_pim;

    TH2D_ptr Theta_prot_thrown_lab_vs_mom_prot;
    TH2D_ptr Theta_pip_thrown_lab_vs_mom_pip;
    TH2D_ptr Theta_pim_thrown_lab_vs_mom_pim;

    TH1D_ptr Phi_gamma;
    TH1D_ptr Phi_prot;
    TH1D_ptr Phi_pip;
    TH1D_ptr Phi_pim;

    TH1D_ptr alpha_pim;
    TH1D_ptr alpha_pip;
    TH1D_ptr alpha_prot;

    TH1D_ptr theta_prot_mc;
    TH1D_ptr theta_pip_mc;
    TH1D_ptr theta_pim_mc;

    TH1D_ptr Phi_gamma_mc;
    TH1D_ptr Phi_prot_mc;
    TH1D_ptr Phi_pip_mc;
    TH1D_ptr Phi_pim_mc;

    TH1D_ptr alpha_pim_mc;
    TH1D_ptr alpha_pip_mc;
    TH1D_ptr alpha_prot_mc;

    TH1D_ptr theta_prot_thrown;
    TH1D_ptr theta_pip_thrown;
    TH1D_ptr theta_pim_thrown;

    TH1D_ptr Phi_gamma_thrown;
    TH1D_ptr Phi_prot_thrown;
    TH1D_ptr Phi_pip_thrown;
    TH1D_ptr Phi_pim_thrown;

    TH1D_ptr alpha_pim_thrown;
    TH1D_ptr alpha_pip_thrown;
    TH1D_ptr alpha_prot_thrown;

    TH2D_ptr W_vs_q2_sec[num_sectors];
    TH1D_ptr W_sec[num_sectors];

    TH1D_ptr W_det[3];
    TH2D_ptr WQ2_det[3];

    TH1D_ptr W_hist_singlePi;
    TH1D_ptr Q2_hist_singlePi;
    TH2D_ptr W_vs_q2_singlePi;
    TH2D_ptr W_vs_q2_singlePi_sec[num_sectors];
    TH1D_ptr W_singlePi_sec[num_sectors];
    TH2D_ptr W_vs_MM_singlePi[num_sectors];

    TH2D_ptr W_vs_q2_Npip_sec[num_sectors];
    TH1D_ptr W_Npip_sec[num_sectors];
    TH1D_ptr MM_Npip_sec[num_sectors];

    TH1D_ptr MM_neutron;
    TH1D_ptr MM_neutron_sec[num_sectors];

    TH2D_ptr W_vs_MM;
    TH2D_ptr W_vs_MM2;

    //////////////////////////   electron pid cuts ///////////////
    TH1D_ptr htcc_nphe_sec[CUTS][num_sectors];
    TH1D_ptr elec_Chi2pid_sec[CUTS][num_sectors];
    TH1D_ptr vz_sec[CUTS][num_sectors];
    TH2D_ptr ECAL_VS_PCAL[CUTS][num_sectors][9];
    TH2D_ptr SF_VS_MOM[CUTS][num_sectors];
    std::string ECIN_ECOUT_MOM_NAME[9] = {" p<2 GeV ", " 2<p<=3 GeV ", " 3<p<=4 GeV ",
                                          " 4<p<=5 GeV ", " 5<=p<6 GeV ",
                                          " 6<p<=7 GeV ", " 7<p<=8 GeV ",
                                          " 8<p<=9 GeV ", " p>9 GeV "};

    TH2D_ptr W_vs_q2_twoPi_sec[num_sectors];
    TH1D_ptr MM_twoPi;
    TH1D_ptr MM_mPim_twoPi_sec[num_sectors];
    TH1D_ptr MM2_twoPi_excl;
    TH1D_ptr MM_twoPi_excl;

    TH1D_ptr MM2_twoPi_mPim;
    TH1D_ptr MM_twoPi_mPim;
    TH1D_ptr MM2_mPim_twoPi_sec[num_sectors];
    TH1D_ptr MM2_twoPi_missingPip;
    TH1D_ptr MM2_twoPi_missingPip_sec[num_sectors];
    TH1D_ptr MM2_twoPi_missingProt;
    TH1D_ptr MM2_twoPi_missingProt_sec[num_sectors];
    TH1D_ptr W_hist_twoPi;
    TH1D_ptr Q2_hist_twoPi;
    TH2D_ptr W_vs_q2_twoPi;
    TH2D_ptr W_vs_q2_twoPi_thrown;

    TH1D_ptr MMSQ_mPim_hist[q2_bin][15];
    TH1D_ptr MMSQ_mPim_hist_with_cut[q2_bin][15];

    TH1D_ptr MMSQ_mPim_hist_1_comb[q2_bin][15];
    TH1D_ptr MMSQ_mPim_hist_2_comb[q2_bin][15];
    TH1D_ptr MMSQ_mPim_hist_3_comb[q2_bin][15];
    TH1D_ptr MMSQ_mPim_hist_4_or_more_comb[q2_bin][15];

    TH1D_ptr Inv_mass_pPip[q2_bin][15];
    TH1D_ptr Inv_mass_pPim[q2_bin][15];
    TH1D_ptr Inv_mass_pipPim[q2_bin][15];
    TH1D_ptr Alpha_Prot_cm[q2_bin][15];
    TH1D_ptr Alpha_Pip_cm[q2_bin][15];
    TH1D_ptr Alpha_Pim_cm[q2_bin][15];

    static const short FDmomArraySize = 24;
    static const short CDmomArraySize = 15;

    TH1D_ptr dt_prot_fd_hist[FDmomArraySize];
    TH1D_ptr dt_prot_cd_hist[CDmomArraySize];
    TH1D_ptr dt_pip_fd_hist[FDmomArraySize];
    TH1D_ptr dt_pip_cd_hist[CDmomArraySize];
    TH1D_ptr dt_pim_fd_hist[FDmomArraySize];
    TH1D_ptr dt_pim_cd_hist[CDmomArraySize];

    TH2D_ptr W_vs_q2_twoPi_sec_thrown[num_sectors];

    // Mom vs Beta
    TH2D_ptr momvsbeta_hist[particle_num][charge_num][with_id_num];
    TH2D_ptr momvsbeta_hist_prot[2][3];
    TH2D_ptr momvsbeta_hist_pip[2][3];
    TH2D_ptr momvsbeta_hist_pim[2][3];

    // Mom vs Beta

    // Delta T
    TH2D_ptr delta_t_hist[3][2][3];

    // Delta T

public:
    Histogram(const std::string &output_file);
    ~Histogram();

    bool MM_cut(float w, float q2, float mm2)
    {
        int is_mc = 0;
        if (_mc)
        {
            is_mc = 1;
        }
        for (int i = 1; i < q2_bins.size(); ++i)
        {
            if (q2 < q2_bins[i])
            {
                q2_bin_val = i;
                break;
            }
        }

        if ((mm2 < (mmsq_cuts[is_mc][q2_bin_val - 1][0][0] * pow(w, 2) + mmsq_cuts[is_mc][q2_bin_val - 1][0][1] * pow(w, 1) + mmsq_cuts[is_mc][q2_bin_val - 1][0][2])) &&
            (mm2 > (mmsq_cuts[is_mc][q2_bin_val - 1][1][0] * pow(w, 2) + mmsq_cuts[is_mc][q2_bin_val - 1][1][1] * pow(w, 1) + mmsq_cuts[is_mc][q2_bin_val - 1][1][2])))

        {
            // std::cout << "   w  = " << w << "  q2 = " << q2 << "  mm2 = " << mm2 << "  up lim mm2 is =  "
            //           << (mmsq_cuts[q2_bin_val - 1][0][0] * pow(w, 2) + mmsq_cuts[q2_bin_val - 1][0][2] * pow(w, 1) + mmsq_cuts[q2_bin_val - 1][0][2])
            //           << "  q2_bin_val = " << q2_bin_val << "   params  " << mmsq_cuts[q2_bin_val - 1][0][0] << " , " << mmsq_cuts[q2_bin_val - 1][0][1] << std::endl;
            return true;
        }
        else
            return false;
    }

    void populate_eff_check_mPim(const std::shared_ptr<Reaction> &_e, double min_w, double max_w, double min_theta, double max_theta, double min_phi, double max_phi, short index_w, short index_theta, short index_phi);
    // void populate_eff_check_exclusive(const std::shared_ptr<Reaction> &_e, double min, double max, short index_w);
    void Fill_eff_ckeck_mPim(const std::shared_ptr<Reaction> &_e);
    // void Fill_eff_ckeck_exclusive(const std::shared_ptr<Reaction> &_e);
    void writeHist_eff_check_mPim();
    // void writeHist_eff_check_mPim_after_MMSQ_cuts();
    // void writeHist_eff_check_exclusive();
    // void writeHist_eff_check_exclusive_MMSQ_cuts();

    void Fill_histthreeD(const std::shared_ptr<Reaction> &_e);
    void writeHists3D();
    void Fill_hists4D_background(const std::shared_ptr<Reaction> &_e);
    void writehists4D_background();

    void Fill_histSevenD_prot(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_thrown_prot(const std::shared_ptr<MCReaction> &_e);

    void Fill_histSevenD_pim(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_thrown_pim(const std::shared_ptr<MCReaction> &_e);
    void Fill_histSevenD_pip(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_thrown_pip(const std::shared_ptr<MCReaction> &_e);

    void Fill_histSevenD_prot_evt(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_thrown_prot_evt(const std::shared_ptr<MCReaction> &_e);
    void Fill_histSevenD_pim_evt(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_thrown_pim_evt(const std::shared_ptr<MCReaction> &_e);
    void Fill_histSevenD_pip_evt(const std::shared_ptr<Reaction> &_e);
    void Fill_histSevenD_thrown_pip_evt(const std::shared_ptr<MCReaction> &_e);

    void writeHists7D_prot();
    void writeHists7D_thrown_prot();
    void writeHists7D_pim();
    void writeHists7D_thrown_pim();
    void writeHists7D_pip();
    void writeHists7D_thrown_pip();
    void writeHists7D_prot_evt();
    void writeHists7D_thrown_prot_evt();
    void writeHists7D_pip_evt();
    void writeHists7D_thrown_pip_evt();
    void writeHists7D_pim_evt();
    void writeHists7D_thrown_pim_evt();

    //// bin centering corrections
    void Fill_hist1D_thrown_w_q2(const std::shared_ptr<MCReaction> &_e);
    void Fill_hist1D_thrown_inv_mass(const std::shared_ptr<MCReaction> &_e);
    void Fill_hist1D_thrown_theta(const std::shared_ptr<MCReaction> &_e);
    void Fill_hist1D_thrown_alpha(const std::shared_ptr<MCReaction> &_e);

    void writeHists1D_thrown_w_gen();
    void writeHists1D_thrown_q2_gen();

    void writeHists1D_thrown_w_gen_inv_pPip();
    void writeHists1D_thrown_q2_gen_inv_pPip();

    void writeHists1D_thrown_w_gen_inv_pPim();
    void writeHists1D_thrown_q2_gen_inv_pPim();

    void writeHists1D_thrown_w_gen_inv_pipPim();
    void writeHists1D_thrown_q2_gen_inv_pipPim();

    void writeHists1D_thrown_w_gen_th_prot();
    void writeHists1D_thrown_q2_gen_th_prot();

    void writeHists1D_thrown_w_gen_th_pip();
    void writeHists1D_thrown_q2_gen_th_pip();

    void writeHists1D_thrown_w_gen_th_pim();
    void writeHists1D_thrown_q2_gen_th_pim();

    void writeHists1D_thrown_w_gen_al_prot();
    void writeHists1D_thrown_q2_gen_al_prot();

    void writeHists1D_thrown_w_gen_al_pip();
    void writeHists1D_thrown_q2_gen_al_pip();

    void writeHists1D_thrown_w_gen_al_pim();
    void writeHists1D_thrown_q2_gen_al_pim();

    void writeHists1D_thrown_protPip();
    void writeHists1D_thrown_protPim();
    void writeHists1D_thrown_pipPim();
    void writeHists1D_thrown_th_prot();
    void writeHists1D_thrown_th_pip();
    void writeHists1D_thrown_th_pim();
    void writeHists1D_thrown_alpha_prot();
    void writeHists1D_thrown_alpha_pip();
    void writeHists1D_thrown_alpha_pim();
    double CosTheta(int theta_bin_);

    // W and Q^2
    void makeHists_sector();
    void Fill_WvsQ2(const std::shared_ptr<Reaction> &_e);
    // void Fill_WvsQ2(const std::shared_ptr<MCReaction> &_e);
    void Fill_WvsQ2_singlePi(const std::shared_ptr<Reaction> &_e);
    void Fill_WvsQ2_Npip(const std::shared_ptr<Reaction> &_e);
    void Fill_WvsQ2_twoPi(const std::shared_ptr<Reaction> &_e);
    // void Fill_WvsQ2_twoPi(const std::shared_ptr<MCReaction>& _e);
    void Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<MCReaction> &_e);
    void Write_WvsQ2();

    void Fill_W_vs_Q2_all_sec();
    void Fill_W_vs_Q2_thrown();
    void Fill_inv_mass_hist();

    void makeHistMMSQ_mPim();
    void Fill_MMSQ_mPim(const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim();

    void Fill_MMSQ_mPim_1_comb(const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim_1_comb();

    void Fill_MMSQ_mPim_2_comb(const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim_2_comb();

    void Fill_MMSQ_mPim_3_comb(const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim_3_comb();

    void Fill_MMSQ_mPim_4_or_more_comb(float dv2, const std::shared_ptr<Reaction> &_e);
    void writeMMSQ_mPim_4_or_more_comb();

    void write_Inv_Mass_hist();

    // P and E
    // ecectron cuts
    void makeHists_electron_cuts();
    void FillHists_electron_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e);
    void FillHists_electron_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e);

    void FillHists_prot_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_prot_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &prot);

    void FillHists_pip_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_pip_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &pip);
    void FillHists_pim_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_pim_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);

    void Write_Electron_cuts();
    void Write_Hadrons_cuts();
    void write_hist_cd_fid();
    void Fill_hist_cd_fid(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Reaction> &_e, int i);

    void Fill_pi0(const std::shared_ptr<Reaction> &_e);

    void makeHists_MomVsBeta();
    void Fill_MomVsBeta(const std::shared_ptr<Branches12> &data, int part, const std::shared_ptr<Reaction> &_e);
    void Write_MomVsBeta();

    // Delta T
    void makeHists_deltat();
    void Fill_deltat_pip(const std::shared_ptr<Branches12> &data,
                         const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_pim(const std::shared_ptr<Branches12> &data,
                         const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_prot(const std::shared_ptr<Branches12> &data,
                          const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_before_cut(const std::shared_ptr<Branches12> &data,
                                const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_prot_after_cut(const std::shared_ptr<Branches12> &data,
                                    const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_pip_after_cut(const std::shared_ptr<Branches12> &data,
                                   const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_pim_after_cut(const std::shared_ptr<Branches12> &data,
                                   const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Write_deltat();

    //////////////
    void Fill_deltaP_prot(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_pip(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_prot_for_pip(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_pip_for_prot(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_ambi_prot(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_ambi_pip(const std::shared_ptr<Reaction> &_e, double dp);

    void Fill_deltaP_ambi_all_prot(const std::shared_ptr<Reaction> &_e, double dp);
    void Fill_deltaP_ambi_all_pip(const std::shared_ptr<Reaction> &_e, double dp);
    void Write_deltaP();

    void Fill_Entries(int num_entries);
    void Fill_Entries_prot(int num_entries);
    void Fill_Entries_pip(int num_entries);
    void Fill_all_Combi(const std::shared_ptr<Reaction> &_e);
    void Fill_1_Combi(const std::shared_ptr<Reaction> &_e);
    void Fill_2_Combi(const std::shared_ptr<Reaction> &_e);
    void Fill_3_Combi(const std::shared_ptr<Reaction> &_e);
    void Fill_4_or_more_Combi(float dv2, const std::shared_ptr<Reaction> &_e);

    ///////////////

    void makeHistTheta_pim_measured();
    void populate_theta_pim_measured(const std::shared_ptr<Reaction> &_e, double min, double max, short index_theta_pim);
    void Fill_theta_pim_measured(const std::shared_ptr<Reaction> &_e);
    void write_hist_theta_pim_measured();

    void Write();

    // Function to get the momentum range index based on the value of p

    int getMomRange(double p)
    {
        const double boundaries[] = {2, 3, 4, 5, 6, 7, 8, 9};
        const int numBoundaries = sizeof(boundaries) / sizeof(boundaries[0]);

        for (int i = 0; i < numBoundaries; ++i)
        {
            if (p < boundaries[i])
            {
                return i;
            }
        }
        return numBoundaries; // For p >= 9
    }
};

#endif
