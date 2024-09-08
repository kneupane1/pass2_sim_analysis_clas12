#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "mom_corr.hpp"
#include "physics.hpp"
#include "eff_corr.hpp"
#include <TRandom3.h>
#include <cmath>
#include "boost_cms.hpp"
#include <vector>
class Reaction
{
protected:
        std::shared_ptr<Branches12> _data;

        double _beam_energy = 10.6041;

        std::unique_ptr<TLorentzVector> _beam;
        std::unique_ptr<TLorentzVector> _elec;
        std::unique_ptr<TLorentzVector> _gamma;
        std::unique_ptr<TLorentzVector> _target;

        // mom corrections
        std::unique_ptr<TLorentzVector> _mom_corr_elec;

        std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_pim;
        std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_pip;
        std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_prot;

        // std::unique_ptr<TLorentzVector> _x_mu;
        // std::unique_ptr<TLorentzVector> _photons;
        std::vector<std::shared_ptr<TLorentzVector>> _photons;

        // std::unique_ptr<TLorentzVector> _prot;
        // std::unique_ptr<TLorentzVector> _pip;
        // std::unique_ptr<TLorentzVector> _pim;
        std::unique_ptr<TLorentzVector> _boosted_gamma;
        std::unique_ptr<TLorentzVector> _boosted_prot;
        std::unique_ptr<TLorentzVector> _boosted_pip;
        std::unique_ptr<TLorentzVector> _boosted_pim;
        std::unique_ptr<TLorentzVector> _other;
        std::unique_ptr<TLorentzVector> _neutron;

        std::unique_ptr<TLorentzVector> _elecUnSmear;
        std::unique_ptr<TLorentzVector> _protUnSmear;
        std::unique_ptr<TLorentzVector> _pipUnSmear;
        std::unique_ptr<TLorentzVector> _pimUnSmear;

        std::vector<std::unique_ptr<TLorentzVector>> _prot;
        std::vector<std::unique_ptr<TLorentzVector>> _pip;
        std::vector<std::unique_ptr<TLorentzVector>> _pim;

        std::vector<std::unique_ptr<TLorentzVector>> _mom_corr_prot;
        std::vector<std::unique_ptr<TLorentzVector>> _mom_corr_pip;
        std::vector<std::unique_ptr<TLorentzVector>> _mom_corr_pim;

        std::vector<int> _prot_indices;
        std::vector<int> _pip_indices;
        std::vector<int> _pim_indices;

        // TVector3 _prot_Vect3;
        // TVector3 _pip_Vect3;
        // TVector3 _pim_Vect3;

        // std::unique_ptr<TLorentzVector> _missingPim;
        std::unique_ptr<TLorentzVector> _boosted_pim_measured;

        float _weight = NAN;

        bool _mc = false;
        bool _is_eff_corrected = false;

        bool _is_boosted = false;

        bool _hasE = false;
        bool _hasP = false;

        bool _hasPip = false;
        bool _hasPim = false;
        bool _hasOther = false;
        bool _hasNeutron = false;

        short _numPart = 0;
        short _numProt = 0;
        short _numPip = 0;
        short _numPim = 0;
        short _numPos = 0;
        short _numNeg = 0;
        short _numNeutral = 0;
        short _numOther = 0;
        short _numPhoton = 0;
        short _sector = -1;

        bool _is_FD_Prot = false;
        bool _is_CD_Prot = false;

        bool _is_FD_Pip = false;
        bool _is_CD_Pip = false;

        bool _is_FD_Pim = false;
        bool _is_CD_Pim = false;
        bool _is_FD = false;
        bool _is_CD = false;

        bool _is_lower_band = false;
        float _MM_mPim = NAN;
        float _MM2_mPim = NAN;

        float _MM_mpip = NAN;
        float _MM2_mpip = NAN;
        float _MM_mprot = NAN;
        float _MM2_mprot = NAN;
        float _pi0_mass = NAN;
        float _MM2_exclusive = NAN;
        float _MM_exclusive = NAN;
        float _W = NAN;
        float _Q2 = NAN;

        float _W_e_Prot = NAN;
        float _Q2_e_Prot = NAN;

        float _inv_Ppip = NAN;
        float _inv_Ppim = NAN;
        float _inv_pip_pim = NAN;
        float _W_P2pi = NAN;

        float _phi_gamma = NAN;
        float _phi_prot = NAN;
        float _phi_pip = NAN;
        float _phi_pim = NAN;

        float _alpha_ppip_pipim = NAN;
        float _alpha_pippim_pipf = NAN;
        float _alpha_ppim_pipip = NAN;

        float _x_mu_Px = NAN;
        float _x_mu_Py = NAN;
        float _x_mu_Pz = NAN;

        float _x_mu_m = NAN;
        float _x_mu_m2 = NAN;
        float _x_mu_P = NAN;
        float _x_mu_E = NAN;
        float _x_mu_theta = NAN;
        float _beam_theta = NAN;
        float _elec_theta = NAN;
        float _E_elec = NAN;
        float _pim_theta_measured = NAN;
        short _pim_sec = -9999;

        float _theta_e = NAN;

        float _P_elec = NAN;
        float _elec_status;
        float _prot_status = NAN;
        float _pip_status = NAN;
        float _pim_status = NAN;

        int _sectorElec = -1;
        int _sectorPim = -1;
        int _sectorPip = -1;
        int _sectorProt = -1;
        // elec mom corr
        double _cx = NAN;
        double _cy = NAN;
        double _cz = NAN;

        double _px_prime_elec = NAN;
        double _py_prime_elec = NAN;
        double _pz_prime_elec = NAN;

        // energy loss corr

        double _px_prime_prot_E = NAN;
        double _py_prime_prot_E = NAN;
        double _pz_prime_prot_E = NAN;

        double _prot_mom = NAN;
        double _prot_mom_uncorr = NAN;
        float _E_corr_val_prot = NAN;
        double _prot_theta_uncorr = NAN;
        float _prot_phi_uncorr = NAN;
        double _prot_mom_tmt = NAN;
        double _prot_mom_prime = NAN;
        double _prot_theta = NAN;
        double _prot_phi = NAN;

        double _px_prime_pip_E = NAN;
        double _py_prime_pip_E = NAN;
        double _pz_prime_pip_E = NAN;

        double _pip_mom = NAN;
        double _pip_mom_tmt = NAN;
        double _pip_mom_uncorr = NAN;
        float _E_corr_val_pip = NAN;
        double _pip_theta_uncorr = NAN;
        double _pip_phi_uncorr = NAN;
        double _pip_mom_prime = NAN;

        double _pip_theta = NAN;
        double _pip_phi = NAN;

        double _px_prime_pim_E = NAN;
        double _py_prime_pim_E = NAN;
        double _pz_prime_pim_E = NAN;

        double _pim_mom = NAN;
        double _pim_mom_prime = NAN;
        double _pim_mom_tmt = NAN;
        double _pim_mom_uncorr = NAN;
        float _E_corr_val_pim = NAN;
        double _pim_theta_uncorr = NAN;
        double _pim_phi_uncorr = NAN;
        double _pim_theta = NAN;
        double _pim_phi = NAN;

        double _px_prime_prot_mom = NAN;
        double _py_prime_prot_mom = NAN;
        double _pz_prime_prot_mom = NAN;

        double _px_prime_pip_mom = NAN;
        double _py_prime_pip_mom = NAN;
        double _pz_prime_pip_mom = NAN;

        double _px_prime_pim_mom = NAN;
        double _py_prime_pim_mom = NAN;
        double _pz_prime_pim_mom = NAN;

        float _excl_Energy = NAN;
        float lcut_me = -0.11;
        float rcut_me = 0.25;
        float lcut_pim = -0.18;
        float rcut_pim = 0.20;
        float lcut_pip = -0.20;
        float rcut_pip = 0.28;
        float lcut_Prot = 0.80; // original is 0.80
        float rcut_Prot = 1.06; // original is 1.06

        double fe = NAN;
        double fpro = NAN;
        double fpip = NAN;
        double fpim = NAN;
        float _thetaDC_r1_Prot = NAN;
        float _thetaDC_r1_Pip = NAN;
        float _thetaDC_r1_Pim = NAN;

        float _pr_p = NAN;
        float _pr_th = NAN;
        float _pr_ph_eff = NAN;

        float _pip_p = NAN;
        float _pip_th = NAN;
        float _pip_ph_eff = NAN;

        float _pim_p = NAN;
        float _pim_th = NAN;
        float _pim_ph_eff = NAN;

        float _eff_corr_fact_mPim = NAN;
        float _eff_corr_fact_Excl = NAN;
        //
        static const int CD_SEC = 3;
        static const int FD_SEC = 6;

        void SetElec();

        double _elec_mom_corrected = NAN;
        double _elec_mom = NAN;

public:
        Reaction(){};
        Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
        ~Reaction();

        inline bool mc()
        {
                return _mc;
        }

        void SetProton(int i);
        void SetPip(int i);
        void SetPim(int i);
        void SetOther(int i);
        void SetNeutron(int i);

        const std::vector<std::unique_ptr<TLorentzVector>> &GetProtons() const { return _prot; }
        const std::vector<std::unique_ptr<TLorentzVector>> &GetPips() const { return _pip; }
        const std::vector<std::unique_ptr<TLorentzVector>> &GetPims() const { return _pim; }
        const std::vector<int> &GetProtonIndices() const { return _prot_indices; }
        const std::vector<int> &GetPipIndices() const { return _pip_indices; }
        const std::vector<int> &GetPimIndices() const { return _pim_indices; }

        short pim_sec();
        bool ctof_prot();
        bool ftof_prot();
        bool ctof_pip();
        bool ftof_pip();
        bool ctof_pim();
        bool ftof_pim();

        // void Prot_HMom_corr(int status_prot, int ststus_pip, int ststus_pim, int sector_Prot, float alFD[4], float alCD[3]);
        // void Pip_HMom_corr(int status_prot, int ststus_pip, int ststus_pim, int sector_Pip, float alFD[4], float alCD[3]);
        // void Pim_HMom_corr(int status_prot, int ststus_pip, int ststus_pim, int sector_Pim, float alFD[4], float alCD[3]);

        // momentum correction
        // void SetMomCorrElec();

        // double dpp(float px, float py, float pz, int sec_mom_corr, int ivec);
        // double Corr_elec_mom();
        double elec_mom();

        // float EffCorrFactor();

        void boost(const TLorentzVector &prot, const TLorentzVector &pip);
        // void CalcMissMass();
        void CalcMissMassPim(const TLorentzVector &prot, const TLorentzVector &pip);
        void CalcMissMassExcl(const TLorentzVector &prot, const TLorentzVector &pip, const TLorentzVector &pim);
        float MM_mPim();
        float MM2_mPim();
        float MM2_mpip();
        float MM2_mprot();
        void CalcMassPi0();
        float pi0_mass();
        float MM2_exclusive();
        float MM_exclusive();
        // float weight();
        inline float weight()
        {
                return _data->mc_weight(); //
                // return 1.0;
        }

        // Check lists when you swich from mc to exp or vice-versa
        // 1. inline weight function above and also in class mc, careful when you do wt =1 for sim to calculate erreors
        // 2. check the deltat cuts and sampling fraction cuts

        // 3. gamma, _w, _q2 and dpp function in electron four vector set up at reaction.cpp because of momentum corrections
        // for elec included only for exp data

        // event->SetMomCorrElec(); please adjust this one also
        // remember to include or exclude momentum corrections for all four particles and our hadron corr too

        //4. check if you have momentum corrected four vectors different
        //5. Change mmsq_low and mmsq_high for MM_cut function in histogram.hpp file

        //6. check it at boosted four vector, if you need mom corrected
        // also:  TVector3 ux = ((_beam->Vect()).Cross(_mom_corr_elec->Vect())).Unit(); // unit vector along e cross e'

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /// smearing fx's function
        void SmearingFunc(int part_id, int status_part, double p, double theta, double phi, double &pNew, double &thetaNew,
                          double &phiNew)
        {
                // Constants
                const double pS1 = 0.0184291 - 0.0110083 * theta + 0.00227667 * pow(theta, 2) - 0.000140152 * pow(theta, 3) +
                                   3.07424e-6 * pow(theta, 4);
                const double pR = 0.02 * sqrt(pow(pS1 * p, 2) + pow(0.02 * theta, 2));
                const double thetaR = 2.5 * sqrt(pow((0.004 * theta + 0.1) * (pow(p, 2) + 0.13957 * 0.13957) / pow(p, 2), 2));
                const double phiS1 = 0.85 - 0.015 * theta;
                const double phiS2 = 0.17 - 0.003 * theta;
                const double phiR = 3.5 * sqrt(pow(phiS1 * sqrt(pow(p, 2) + 0.13957 * 0.13957) / pow(p, 2), 2) + pow(phiS2, 2));

                // Generate new values
                if (part_id == ELECTRON)
                {
                        phiNew = phi + 0.4 * phiR * gRandom->Gaus(0, 1);
                        thetaNew = theta + 0.4 * thetaR * gRandom->Gaus(0, 1);
                        pNew = p + 0.4 * pR * gRandom->Gaus(0, 1) * p;
                }
                else if (part_id == PROTON)
                {
                        double fact_cd = 0;
                        double fact_fd = 0;
                        double fact_cd1 = 0;
                        double fact_fd1 = 0;
                        if (status_part > 4000)
                        {
                                fact_cd = (0.000821) * pow(p, 3) + (-0.016500) * pow(p, 2) + (0.103611) * p + (1.393237);
                                fact_cd1 = 1.0; //(0.001536) * pow(p, 3) + (-0.024778) * pow(p, 2) + (0.119853) * p + (0.832939);

                                phiNew = phi + 1 / (fact_cd)*phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_cd)*thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_cd)*pR * gRandom->Gaus(0, 1) * p;
                                // std::cout << "mom " << p << "prot fact_cd : " << 1 / fact_cd << std::endl;
                        }
                        else if (status_part <= 4000)
                        {
                                fact_fd = (0.000264) * pow(p, 3) + (-0.006454) * pow(p, 2) + (0.032683) * p + (1.658142);
                                fact_fd1 = 1.0; //(0.000051) * pow(p, 3) + (-0.001569) * pow(p, 2) + (0.015891) * p + (0.966351);

                                phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;

                                // std::cout << "mom " << p << "prot fact_fd : " << 1 / fact_fd << std::endl;
                        }
                }

                else if (part_id == PIP)
                {
                        double fact_cd = 0;
                        double fact_fd = 0;
                        double fact_cd1 = 0;
                        double fact_fd1 = 0;
                        if (status_part > 4000)
                        {
                                fact_cd = (0.000981) * pow(p, 3) + (-0.016882) * pow(p, 2) + (0.046752) * p + (1.720426);
                                fact_cd1 = 1.0; //(-0.000104) * pow(p, 3) + (0.000998) * pow(p, 2) + (-0.008019) * p + (1.105314);

                                // std::cout << "mom " << p << "pip fact_cd : " << 1 / fact_cd << std::endl;

                                phiNew = phi + 1 / (fact_cd * fact_cd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_cd * fact_cd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_cd * fact_cd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                        else if (status_part <= 4000)
                        {
                                fact_fd = (0.000085) * pow(p, 3) + (-0.003096) * pow(p, 2) + (0.023553) * p + (1.509910);
                                fact_fd1 = 1.0; //(-0.000006) * pow(p, 3) + (-0.001310) * pow(p, 2) + (0.023171) * p + (0.890554);

                                // std::cout << "mom " << p << "pip fact_fd : " << 1 / fact_fd << std::endl;

                                phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                }

                else if (part_id == PIM)
                {
                        double fact_cd = 0;
                        double fact_fd = 0;
                        double fact_cd1 = 0;
                        double fact_fd1 = 0;
                        if (status_part > 4000)
                        {
                                fact_cd = (-0.001788) * pow(p, 3) + (0.025796) * pow(p, 2) + (-0.136577) * p + (2.007917);
                                fact_cd1 = 1.0; //(-0.001327) * pow(p, 3) + (0.019826) * pow(p, 2) + (-0.097667) * p + (1.308904);
                                // std::cout << "mom " << p << "pim fact_cd : " << 1 / fact_cd << std::endl;

                                phiNew = phi + 1 / (fact_cd * fact_cd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_cd * fact_cd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_cd * fact_cd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                        else if (status_part <= 4000)
                        {
                                fact_fd = (0.000760) * pow(p, 3) + (-0.021295) * pow(p, 2) + (0.171180) * p + (1.238299);
                                fact_fd1 = 1.0; //(0.000249) * pow(p, 3) + (-0.007461) * pow(p, 2) + (0.067686) * p + (0.817653);

                                // std::cout << "mom " << p << "pim fact_cd : " << 1 / fact_fd << std::endl;

                                phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                }
        }

        float prot_theta();
        float pip_theta();
        float pim_theta();
        float inv_Ppip();
        float inv_Ppim();
        float inv_pip_pim();
        float w_P2pi_rec();

        void W_2pi_P();
        void invMassPpip();
        void invMassPpim();
        void invMasspippim();
        float gamma_Phi();
        float prot_Phi();
        float pip_Phi();
        float pim_Phi();

        float alpha_ppip_pipim();
        float alpha_pippim_pipf();
        float alpha_ppim_pipip();
        ///////////// related to lab system before boost ////////////
        float prot_momentum(const TLorentzVector &prot);
        float prot_theta_lab(const TLorentzVector &prot);
        float prot_Phi_lab(const TLorentzVector &prot);
        float prot_momT(const TLorentzVector &prot);

        float pip_momentum(const TLorentzVector &pip);
        float pip_theta_lab(const TLorentzVector &pip);
        float pip_Phi_lab(const TLorentzVector &pip);
        float pip_momT(const TLorentzVector &pip);

        float pim_momentum(const TLorentzVector &prot, const TLorentzVector &pip);
        float pim_momentum_measured(const TLorentzVector &pim);
        float pim_theta_lab(const TLorentzVector &prot, const TLorentzVector &pip);
        float pim_theta_lab_measured(const TLorentzVector &pim);
        float pim_E(const TLorentzVector &prot, const TLorentzVector &pip);
        float pim_E_measured(const TLorentzVector &prot);
        float pim_Phi_lab(const TLorentzVector &prot, const TLorentzVector &pip);
        float pim_Phi_lab_measured(const TLorentzVector &pim);
        float pim_momT(const TLorentzVector &pim);

        float elec_momentum();
        float theta_beam();
        float theta_elec();
        float Phi_elec();
        float E_elec();

        float Energy_excl();

        float AlphaCalc();

        inline float W()
        {
                return _W;
        }
        inline float Q2()
        {
                return _Q2;
        }
        //
        // float theta_x_mu();
        // float E_x_mu();
        // float P_x_mu();
        // float theta_beam();

        inline short sec()
        {
                return _data->dc_sec(0);
        }
        inline int cc_tot_nphe()
        {
                return _data->cc_nphe_tot(0);
        }
        inline int vz()
        {
                return _data->vz(0);
        }
        inline int det()
        {
                return abs(_data->status(0) / 1000);
        }

        inline bool Fixed_MM_cut()
        {
                //exp
                return (Reaction::MM2_mPim() < 0.08 && Reaction::MM2_mPim() > -0.06);
        }

        inline bool TwoPion_missingPim()
        {
                bool _channelTwoPi = true;
                _channelTwoPi &= ((_numProt >= 1 && _numPip >= 1) &&
                                  (_hasE && _hasP && _hasPip
                                   /*&&!_hasPim && !_hasNeutron
                              &&!_hasOther*/
                                   ));
                return _channelTwoPi;
        }
        inline bool TwoPion_exclusive()
        {
                bool _channelTwoPi_excl = true;

                _channelTwoPi_excl &= ((_numProt >= 1 && _numPip >= 1 && _numPim >= 1) &&
                                       (_hasE && _hasP && _hasPip && _hasPim /*&& !_hasNeutron && !_hasOther*/));
                return _channelTwoPi_excl;
        }
        inline bool TwoPion_missingPip()
        {
                bool _channelTwoPi_mpip = true;

                _channelTwoPi_mpip &= ((_numProt >= 1 && _numPim >= 1) &&
                                       (_hasE && _hasP &&
                                        _hasPim /*&&!_hasPip && !_hasNeutron && !_hasOther*/));
                return _channelTwoPi_mpip;
        }
        inline bool TwoPion_missingProt()
        {
                bool _channelTwoPi_mprot = true;
                _channelTwoPi_mprot &= ((_numPip >= 1 && _numPim >= 1) &&
                                        (_hasE && _hasPip && _hasPim /*&&!_hasP  && !_hasOther*/));
                return _channelTwoPi_mprot;
        }
        inline bool inclusive()
        {
                return (_hasE);
        }
        //
        // inline bool NeutronPip() {
        //   bool _channel = true;
        //   _channel &= ((_numPip == 1) && (_hasE && !_hasP && _hasPip && !_hasPim &&
        //   _hasNeutron)) ||
        //               (Reaction::SinglePip() && Reaction::MM() >= 0.85 &&
        //               Reaction::MM() <= 1.0);
        //   return _channel;
        // }

        const TLorentzVector &e_mu()
        {
                return *_beam;
        }
        const TLorentzVector &e_mu_prime()
        {
                return *_elec;
        }
        const TLorentzVector &gamma()
        {
                return *_gamma;
        }
};

class MCReaction : public Reaction
{
private:
        float _weight_mc = NAN;
        float _W_mc = NAN;
        float _Q2_mc = NAN;

        // float _MCinv_Ppip = NAN;
        // float _MCinv_Ppim = NAN;
        // float _MCinv_pip_pim = NAN;

        std::unique_ptr<TLorentzVector> _beam_mc;
        std::unique_ptr<TLorentzVector> _elec_mc;
        std::unique_ptr<TLorentzVector> _gamma_mc;
        std::unique_ptr<TLorentzVector> _other_mc;

        // std::unique_ptr<TLorentzVector> _prot_mc;
        // std::unique_ptr<TLorentzVector> _pip_mc;
        // std::unique_ptr<TLorentzVector> _pim_mc;

        std::vector<std::unique_ptr<TLorentzVector>> _prot_mc;
        std::vector<std::unique_ptr<TLorentzVector>> _pip_mc;
        std::vector<std::unique_ptr<TLorentzVector>> _pim_mc;

        std::vector<int> _prot_mc_indices;
        std::vector<int> _pip_mc_indices;
        std::vector<int> _pim_mc_indices;

        std::unique_ptr<TLorentzVector> _boosted_gamma_mc;
        std::unique_ptr<TLorentzVector> _boosted_prot_mc;
        std::unique_ptr<TLorentzVector> _boosted_pip_mc;
        std::unique_ptr<TLorentzVector> _boosted_pim_mc;

        bool _is_boosted_mc = false;

        float _MM_mc = NAN;
        float _MM2_mc = NAN;
        float _MM2_exclusive_mc = NAN;

        float _alpha_ppip_pipim_mc = NAN;
        float _alpha_pippim_pipf_mc = NAN;
        float _alpha_ppim_pipip_mc = NAN;

public:
        MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
        void SetMCElec();
        inline float weight()
        {
                return _data->mc_weight();
                // return 1.0;
        }
        inline float W_mc()
        {
                return _W_mc;
        }
        inline float Q2_mc()
        {
                return _Q2_mc;
        }

        inline bool MM_cut_mc()
        {
                //sim
                return (MCReaction::MM2_mPim_mc() < 0.08 && MCReaction::MM2_mPim_mc() > -0.06);
        }
        void CalcMissMass_mc();
        float MM_mPim_mc();
        float MM2_mPim_mc();
        float MM2_exclusive_mc();

        void SetMCProton(int i);
        void SetMCPip(int i);
        void SetMCPim(int i);
        void SetMCOther(int i);

        const std::vector<std::unique_ptr<TLorentzVector>> &GetMcProtons() const { return _prot_mc; }
        const std::vector<std::unique_ptr<TLorentzVector>> &GetMcPips() const { return _pip_mc; }
        const std::vector<std::unique_ptr<TLorentzVector>> &GetMcPims() const { return _pim_mc; }

        const std::vector<int> &GetProtonMcIndices() const { return _prot_mc_indices; }
        const std::vector<int> &GetPipMcIndices() const { return _pip_mc_indices; }
        const std::vector<int> &GetPimMcIndices() const { return _pim_mc_indices; }

        void boost_mc(const TLorentzVector &prot_mc, const TLorentzVector &pip_mc, const TLorentzVector &pim_mc);

        float MCinv_Ppip();
        float MCinv_Ppim();
        float MCinv_pip_pim();
        //
        // float MCprot_theta();
        // float MCpip_theta();
        // float MCpim_theta();
        //
        float MCprot_theta_lab();
        float MCpip_theta_lab();
        float MCpim_theta_lab();

        float prot_momentum_thrown();
        float pip_momentum_thrown();
        float pim_momentum_thrown();
        // float MCgamma_Phi();
        // float MCprot_Phi();
        // float MCpip_Phi();
        // float MCpim_Phi();
        //
        // float MCalpha_ppip_pipim();
        // float MCalpha_pippim_pipf();
        // float MCalpha_ppim_pipip();
        float MCprot_theta_thrown();
        float MCpip_theta_thrown();
        float MCpim_theta_thrown();

        float MCgamma_Phi_thrown();
        float MCprot_Phi_thrown();
        float MCpip_Phi_thrown();
        float MCpim_Phi_thrown();

        void MCAlphaCalc();

        float MCalpha_ppip_pipim_thrown();
        float MCalpha_pippim_pipf_thrown();
        float MCalpha_ppim_pipip_thrown();
};

#endif