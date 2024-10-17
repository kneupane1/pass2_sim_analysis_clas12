
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"
// #include "QADB.h"
bool is_match(std::shared_ptr<Branches12> data, int rec_part, int mc_part)
{
        // Example matching based on momentum similarity
        TVector3 rec_momentum(data->px(rec_part), data->py(rec_part), data->pz(rec_part));
        TVector3 mc_momentum(data->mc_px(mc_part), data->mc_py(mc_part), data->mc_pz(mc_part));

        // Check if the momentum vectors are approximately equal
        if (rec_momentum.Mag() - mc_momentum.Mag() < 0.1)
        { // Adjust tolerance as needed
                return true;
        }

        return false;
}
template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram> &_hists, int thread_id)
// size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram> &_hists,
//            const std::shared_ptr<QA::QADB> &_qa = nullptr, int thread_id)
{
        if (_mc)
        {
                // Monte Carlo data processing (no QA)
                std::cout << "Processing without QA (MC)" << std::endl;
        }
        else
        {
                // Real data processing with QA checks
                std::cout << "Processing with QA" << std::endl;
        }

        // Get the number of events in this thread
        size_t num_of_events = (int)_chain->GetEntries();

        float beam_energy = 10.6041;

        if (std::is_same<CutType, Pass2_Cuts>::value)
        {
                if (thread_id == 0)
                        std::cout << BLUE << "Using Pass2 RGA Cuts" << DEF << std::endl;

                beam_energy = 10.6041;
        }

        if (getenv("CLAS12_E") != NULL)
                beam_energy = atof(getenv("CLAS12_E"));

        // Print some information for each thread
        std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
                  << num_of_events << " Events " << DEF << "===============\n";

        // // // Make a data object which all the branches can be accessed from
        // auto data = std::make_shared<Branches12>(_chain, true);
        auto data = std::make_shared<Branches12>(_chain, _mc);

        //         // if (!_qa->Golden(data->getRun(), data->getEvent()))
        //         //         continue;

        // Total number of events "Processed"
        size_t total = 0;
        float no_of_events = 0, miss_prot = 0, miss_pim = 0, miss_pip = 0, other = 0, excl_events = 0,
              twopi = 0;
        int two_pion_mPim_events = 0, two_pion_Excl_events = 0;
        int prot = 0, pip = 0, pim = 0, elec = 0, failed_prot = 0;
        int more_prot = 0, more_pip = 0, more_pim = 0, prot_pip = 0;
        int no_prot_pip = 0, no_prot_pip_mc = 0;
        int Pip_pid_mc = -9999;
        int Prot_pid_mc = -9999;
        int Pim_pid_mc = -9999;
        int Pip_pid_rec = -9999;
        int Prot_pid_rec = -9999;
        int Pim_pid_rec = -9999;
        double dp_prot1 = NAN;
        double dp_pip1 = NAN;
        double dp_prot2 = NAN;
        double dp_pip2 = NAN;
        double dp_prot3 = NAN;
        double dp_pip3 = NAN;
        double dp_pim1 = NAN;
        double dp_pim2 = NAN;
        double dp_pim3 = NAN;
        // For each event
        for (size_t current_event = 0; current_event < num_of_events; current_event++)
        // for (size_t current_event = 0; current_event < 10; current_event++)
        {

                // Get current event
                _chain->GetEntry(current_event);
                // If we are the 0th thread print the progress of the thread every 1000 events
                if (thread_id == 0 && current_event % 1000 == 0)
                        std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

                /////////// This part is only for mc events /////////////
                /////////// This part is only for mc events /////////////
                /////////// This part is only for mc events /////////////
                /////////// This part is only for mc events /////////////
                /////////// This part is only for mc events /////////////
                // if (_mc)
                // {
                if (data->mc_npart() < 1 || data->mc_weight() <= 0)
                        continue;

                // Make a reaction class from the data given
                auto mc_event = std::make_shared<MCReaction>(data, beam_energy);
                // if (mc_event->weight() <= 0.0)
                //         continue;
                // {
                // std::cout << "  mc  pid at : 3  " << data->mc_pid(3) << std::endl;

                for (int part = 0; part < data->mc_npart(); part++)
                // for (int part = 0; part < 2; part++)

                {
                        // std::cout << part << "  mc_pip  pid  " << data->mc_pid(part) << std::endl;

                        // Check particle ID's and fill the reaction class
                        if (data->mc_pid(part) == PIP)
                        {
                                mc_event->SetMCPip(part);
                        }
                        if (data->mc_pid(part) == PROTON)
                        {
                                mc_event->SetMCProton(part);
                        }
                        if (data->mc_pid(part) == PIM)
                        {
                                mc_event->SetMCPim(part);
                                // } else {
                                //   mc_event->SetMCOther(part);
                        }
                }
                // if (data->mc_pid(3) != PIP)
                // std::cout << "mc_pip  pid  " << data->mc_pid(1) << std::endl;
                //         // if (mc_event->W_mc() < 2.15 && mc_event->W_mc() > 1.3 && mc_event->Q2_mc() < 9.0 && mc_event->Q2_mc() > 1.5 && mc_event->weight() > 0.0)
                {

                        // // Retrieve the number of protons and pions in the event
                        size_t num_protons_mc = mc_event->GetMcProtons().size();
                        size_t num_pips_mc = mc_event->GetMcPips().size();
                        size_t num_pims_mc = mc_event->GetMcPims().size();

                        for (size_t i = 0; i < num_protons_mc; ++i)
                        {
                                for (size_t j = 0; j < num_pips_mc; ++j)
                                {

                                        if (mc_event->GetProtonMcIndices()[i] == mc_event->GetPipMcIndices()[j])
                                                no_prot_pip_mc++;

                                        // Exclude the case where the same particle is assigned as both proton and pip
                                        if (mc_event->GetProtonMcIndices()[i] != mc_event->GetPipMcIndices()[j])
                                        {

                                                for (size_t k = 0; k < num_pims_mc; ++k)
                                                {

                                                        mc_event->boost_mc(*mc_event->GetMcProtons()[i], *mc_event->GetMcPips()[j], *mc_event->GetMcPims()[k]);

                                                        _hists->Fill_WvsQ2_twoPi_thrown(mc_event);
                                                        _hists->Fill_histSevenD_thrown_pim(mc_event);
                                                        _hists->Fill_histSevenD_thrown_prot(mc_event);
                                                        _hists->Fill_histSevenD_thrown_pip(mc_event);
                                                        _hists->Fill_histSevenD_thrown_pim_evt(mc_event);
                                                        _hists->Fill_histSevenD_thrown_prot_evt(mc_event);
                                                        _hists->Fill_histSevenD_thrown_pip_evt(mc_event);
                                                        //         // // // _hists->Fill_W_bin_check_th(mc_event);

                                                        // // // ///// bin centering corr
                                                        // _hists->Fill_hist1D_thrown_w_q2(mc_event);
                                                        // _hists->Fill_hist1D_thrown_inv_mass(mc_event);
                                                        // _hists->Fill_hist1D_thrown_theta(mc_event);
                                                        // _hists->Fill_hist1D_thrown_alpha(mc_event);
                                                }
                                        }
                                }
                        }
                }
                // }
                ///////// This part is only for Rec events both mc and exp/////////////
                ///////// This part is only for Rec events both mc and exp/////////////
                ///////// This part is only for Rec events both mc and exp/////////////
                ///////// This part is only for Rec events both mc and exp/////////////

                auto event = std::make_shared<Reaction>(data, beam_energy);
                no_of_events++;
                if (data->charge(0) == -1)
                        _hists->FillHists_electron_cuts(data, event);

                // auto cuts = std::make_unique<rga_Cuts>(data);
                auto cuts = std::make_unique<Pass2_Cuts>(data);

                if (!cuts->ElectronCuts())
                        continue;
                // event->SetMomCorrElec();

                _hists->FillHists_electron_with_cuts(data, event);

                // If we pass electron cuts the event is processed
                elec++;
                // total++;
                int statusPim = -99;
                int statusPip = -99;
                int statusProt = -99;
                int sectorPim = -99;
                int sectorPip = -99;
                int sectorProt = -99;
                int prot_idx = -1;
                int pip_idx = -1;

                // // Make a reaction class from the data
                auto dt = std::make_shared<Delta_T>(data);

                // For each particle in the event
                int more_than_1_prot = 0, more_than_1_pip = 0, more_than_1_pim = 0, both_prot_pip = 0;

                for (int part = 1; part < data->gpart(); part++)
                {
                        // std::cout << "mc_pip  pid  " << data->mc_pid(5) << std::endl;
                        // std::cout << "  rec particle no is :   " << part << std::endl;
                        if (data->charge(part) != 0)
                        {
                                // std::cout << "  rec particle no is :   " << part << std::endl;

                                // total_part++;
                                dt->dt_calc(part);
                                // //dt->dt_calc_ctof(part);
                                // _hists->Fill_MomVsBeta(data, part, event);
                                // // //if (event->TwoPion_missingPim()) {
                                // // // if (event->TwoPion_missingPim()) {
                                // // if(data->status(part) >=2000 && data->status(part) < 4000) {
                                // // _hists->Fill_deltat_pi(data, dt, part, event);
                                // // //}
                                // _hists->Fill_deltat_before_cut(data, dt, part, event);
                                // // _hists->Fill_deltat_prot(data, dt, part, event);
                                // // //}
                                // // // }
                                // // Check particle ID's and fill the reaction class
                                if (data->charge(part) > 0)
                                {
                                        _hists->FillHists_prot_pid_cuts(data, event, part);
                                        _hists->FillHists_pip_pid_cuts(data, event, part);

                                        if ((cuts->IsProton(part)) && (cuts->IsPip(part)))
                                        {

                                                // std::cout << "event is :  " << current_event << "   part is  " << part << std::endl;

                                                // dp_prot1 = (pow((mc_event->prot_momX_mc_gen() - data->px(part)), 2) +
                                                //             pow((mc_event->prot_momY_mc_gen() - data->py(part)), 2) +
                                                //             pow((mc_event->prot_momZ_mc_gen() - data->pz(part)), 2));

                                                // dp_pip1 = (pow((mc_event->pip_momX_mc_gen() - data->px(part)), 2) +
                                                //            pow((mc_event->pip_momY_mc_gen() - data->py(part)), 2) +
                                                //            pow((mc_event->pip_momZ_mc_gen() - data->pz(part)), 2));

                                                // _hists->Fill_deltaP_ambi_all_prot(event, dp_prot1);
                                                // _hists->Fill_deltaP_ambi_all_pip(event, dp_pip1);
                                                // if ((dp_prot1 < 0.1) || (dp_pip1 < 0.1))
                                                {
                                                        // if (dp_prot1 < dp_pip1)
                                                        {
                                                                // _hists->Fill_deltaP_ambi_prot(event, dp_prot1);

                                                                event->SetProton(part);
                                                                // prot_idx++;
                                                                statusProt = abs(data->status(part));
                                                                sectorProt = data->dc_sec(part);
                                                                // _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                                                                // _hists->FillHists_prot_pid_with_cuts(data, event, part, *event->GetProtons()[prot_idx]);
                                                                prot++;
                                                        }
                                                        // else
                                                        {
                                                                // _hists->Fill_deltaP_ambi_pip(event, dp_pip1);
                                                                event->SetPip(part);
                                                                // pip_idx++;
                                                                statusPip = abs(data->status(part));
                                                                sectorPip = data->dc_sec(part);
                                                                // _hists->Fill_deltat_pip_after_cut(data, dt, part, event);
                                                                // _hists->FillHists_pip_pid_with_cuts(data, event, part, *event->GetPips()[pip_idx]);
                                                        }
                                                }
                                        }

                                        else if (cuts->IsProton(part))
                                        {

                                                prot++;
                                                // more_than_1_prot++;
                                                // // std::cout << "event is :  " << current_event << "   part is  " << part << "  mre than 1 protons  =  " << more_than_1_prot << std::endl;

                                                // dp_prot2 = (pow((mc_event->prot_momX_mc_gen() - data->px(part)), 2) +
                                                //             pow((mc_event->prot_momY_mc_gen() - data->py(part)), 2) +
                                                //             pow((mc_event->prot_momZ_mc_gen() - data->pz(part)), 2));
                                                // dp_pip2 = (pow((mc_event->pip_momX_mc_gen() - data->px(part)), 2) +
                                                //            pow((mc_event->pip_momY_mc_gen() - data->py(part)), 2) +
                                                //            pow((mc_event->pip_momZ_mc_gen() - data->pz(part)), 2));
                                                // // std::cout << "event is :  " << current_event << "   part is  " << part << "  diff = " << dp_pip1 << "  && "
                                                // //           << dp_pip2 << std::endl;

                                                // _hists->Fill_deltaP_prot(event, dp_prot2);
                                                // _hists->Fill_deltaP_pip_for_prot(event, dp_pip2);

                                                event->SetProton(part);
                                                // prot_idx++;

                                                statusProt = abs(data->status(part));
                                                sectorProt = data->dc_sec(part);

                                                // _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                                                // _hists->FillHists_prot_pid_with_cuts(data, event, part, *event->GetProtons()[prot_idx]);
                                        }

                                        else if (cuts->IsPip(part))
                                        {
                                                pip++;
                                                // more_than_1_pip++;

                                                // dp_prot3 = (pow((mc_event->prot_momX_mc_gen() - data->px(part)), 2) + pow((mc_event->prot_momY_mc_gen() - data->py(part)), 2) +
                                                //             pow((mc_event->prot_momZ_mc_gen() - data->pz(part)), 2));
                                                // dp_pip3 = (pow((mc_event->pip_momX_mc_gen() - data->px(part)), 2) +
                                                //            pow((mc_event->pip_momY_mc_gen() - data->py(part)), 2) +
                                                //            pow((mc_event->pip_momZ_mc_gen() - data->pz(part)), 2));

                                                // // if (more_than_1_pip > 1)
                                                // // {
                                                // // std::cout << "event is :  " << current_event << "   part is  " << part << "  diff = " << dp_pip1 << "  && "
                                                // //           << dp_pip2 << std::endl;
                                                // _hists->Fill_deltaP_prot_for_pip(event, dp_prot3);
                                                // _hists->Fill_deltaP_pip(event, dp_pip3);
                                                // }

                                                event->SetPip(part);

                                                // pip_idx++;

                                                statusPip = abs(data->status(part));
                                                sectorPip = data->dc_sec(part);

                                                // _hists->Fill_deltat_pip_after_cut(data, dt, part, event);
                                                // _hists->FillHists_pip_pid_with_cuts(data, event, part, *event->GetPips()[pip_idx]);
                                        }
                                }
                                else
                                {
                                        _hists->FillHists_pim_pid_cuts(data, event, part);

                                        if (cuts->IsPim(part))

                                        {
                                                event->SetPim(part);

                                                pim++;
                                                more_than_1_pim++;

                                                statusPim = abs(data->status(part));
                                                sectorPim = data->dc_sec(part);
                                                // _hists->FillHists_pim_pid_cuts(data, event, part);

                                                {
                                                        _hists->Fill_deltat_pim_after_cut(data, dt, part, event);
                                                        _hists->FillHists_pim_pid_with_cuts(data, event, part);
                                                }
                                        }
                                }

                                // } else if (cuts->IsmissingPim(part)) {
                                //   if (event->MM_cut()) event->SetmissingPim(part);

                                // else
                                // {
                                //         other++;
                                //         event->SetOther(part);
                                // }
                                // // }

                                if ((cuts->IsPip(part)) && (cuts->IsProton(part)))
                                {
                                        both_prot_pip++;
                                }
                        }
                }
                // /////////////////////////
                // if (more_than_1_prot > 1)
                // {
                //         more_prot++;
                // }
                // if (more_than_1_pip > 1)
                // {
                //         more_pip++;
                // }
                // if (more_than_1_pim > 1)
                // {
                //         more_pim++;
                // }
                // if (both_prot_pip >= 1)
                // {
                //         prot_pip++;
                // }
                // std::cout << "no of events = " << no_of_events << '\n'
                //           << " elec " << elec << " prot " << prot << " pip " << pip << " pim " << pim << '\n';
                // //  if (event->TwoPion_missingPim()) {
                // for (int part = 1; part < data->gpart(); part++) {
                //   if (event->MM_cut()) event->SetmissingPim(part);
                //   //  }
                // }

                // Check the reaction class what kind of event it is and fill the appropriate histograms

                if (event->W() > 1.35 && event->W() <= 2.15 && event->Q2() <= 9.0 && event->Q2() >= 1.95 && event->weight() > 0.0)
                {
                        //  if ((Prot_pid_mc == PROTON) && (Pip_pid_mc == PIP))
                        {

                                // // if (event->TwoPion_exclusive())
                                // // {
                                // //         for (size_t i = 0; i < event->GetProtons().size(); ++i)
                                // //         {
                                // //                 for (size_t j = 0; j < event->GetPips().size(); ++j)
                                // //                 {
                                // for (size_t k = 0; k < event->GetPims().size(); ++k)
                                // {
                                //         //                                 event->CalcMissMassExcl(*event->GetProtons()[i], *event->GetPips()[j], *event->GetPims()[k]);

                                //         //                                 two_pion_Excl_events++;
                                //         //                                 _hists->Fill_WvsQ2(event);

                                //         // // You should have a similar method for π⁻ if applicable
                                //         dt->dt_calc(event->GetPimIndices()[k]);
                                //         _hists->Fill_deltat_pim_after_cut(data, dt, event->GetPimIndices()[k], event);
                                //         _hists->FillHists_pim_pid_with_cuts(data, event, event->GetPimIndices()[k]);
                                // }
                                // //                 }
                                // //         }
                                // // }

                                if (event->TwoPion_missingPim())
                                { // // Retrieve the number of protons and pions in the event
                                        int entries_in_this_event = 0;

                                        int num_protons = event->GetProtons().size();
                                        int num_pips = event->GetPips().size();

                                        _hists->Fill_Entries_prot(num_protons);
                                        _hists->Fill_Entries_pip(num_pips);

                                        // std::cout << "Number of particles: " << data->gpart() << std::endl;

                                        // if (num_protons > 1)
                                        //         std::cout << "no of protons  : " << num_protons << std::endl;
                                        // if (num_pips > 1)
                                        //         std::cout << "no of pip      : " << num_pips << std::endl;

                                        // // Calculate the total number of combinations
                                        // size_t num_combinations = num_protons * num_pips;

                                        // // Get the vector of protons
                                        // const auto &protons = event->GetProtons();
                                        // const auto &pip = event->GetPips();

                                        for (size_t i = 0; i < num_protons; ++i)
                                        {
                                                for (size_t j = 0; j < num_pips; ++j)
                                                {
                                                        if (event->GetProtonIndices()[i] == event->GetPipIndices()[j])
                                                                no_prot_pip++;

                                                        // // // Print//////////////////////////////////////
                                                        // std::cout << "Event " << current_event << ": "
                                                        //           << num_protons << " proton(s), "
                                                        //           << num_pips << " pip(s), "
                                                        //           << num_combinations << " combination(s)." << std::endl;

                                                        // if (both_prot_pip >= 1)
                                                        // Access the i-th proton and j-th pip

                                                        // Exclude the case where the same particle is assigned as both proton and pip
                                                        if (event->GetProtonIndices()[i] != event->GetPipIndices()[j])
                                                        {
                                                                event->CalcMissMassPim(*event->GetProtons()[i], *event->GetPips()[j]);
                                                                event->boost(*event->GetProtons()[i], *event->GetPips()[j]);
                                                                entries_in_this_event++;

                                                                if (entries_in_this_event == 1)
                                                                {
                                                                        _hists->Fill_MMSQ_mPim(event);

                                                                        // if (event->Fixed_MM_cut())
                                                                        // if (MM_cut(event->W(), event->Q2(), event->MM2_mPim()))
                                                                        // if ((data->p(event->GetProtonIndices()[i]) > 3.0) || (data->p(event->GetPipIndices()[j]) > 3.0))

                                                                        {
                                                                                // _hists->FillHists_electron_with_cuts(data, event);

                                                                                two_pion_mPim_events++;
                                                                                // {
                                                                                // _hists->Fill_MMSQ_mPim(event);
                                                                                // if (entries_in_this_event == 1)
                                                                                _hists->Fill_WvsQ2(event);

                                                                                // // if (data->p(event->GetProtonIndices()[i]) > 3.0)
                                                                                // {
                                                                                //         dt->dt_calc(event->GetProtonIndices()[i]);
                                                                                //         _hists->Fill_deltat_prot_after_cut(data, dt, event->GetProtonIndices()[i], event);
                                                                                //         _hists->FillHists_prot_pid_with_cuts(data, event, event->GetProtonIndices()[i], *event->GetProtons()[i]);
                                                                                // }
                                                                                // // if (data->p(event->GetPipIndices()[j]) > 3.0)
                                                                                // {
                                                                                //         dt->dt_calc(event->GetPipIndices()[j]);
                                                                                //         _hists->Fill_deltat_pip_after_cut(data, dt, event->GetPipIndices()[j], event);
                                                                                //         _hists->FillHists_pip_pid_with_cuts(data, event, event->GetPipIndices()[j], *event->GetPips()[j]);
                                                                                // }

                                                                                _hists->Fill_histSevenD_prot(event);
                                                                                _hists->Fill_histSevenD_pip(event);
                                                                                _hists->Fill_histSevenD_pim(event);

                                                                                // _hists->Fill_histSevenD_prot_evt(event);
                                                                                // _hists->Fill_histSevenD_pip_evt(event);
                                                                                // _hists->Fill_histSevenD_pim_evt(event);

                                                                                // // //         }
                                                                                // // // }
                                                                                // // //         else if (event->TwoPion_missingPip())
                                                                                // // //                 miss_pip++;
                                                                                // // //         else if (event->TwoPion_missingProt())
                                                                                // // //                 miss_prot++;

                                                                                // // //         //{
                                                                                // // //         twopi++;
                                                                                for (int part = 1; part < data->gpart(); part++)
                                                                                {

                                                                                        if (data->charge(part) != 0)
                                                                                        {

                                                                                                dt->dt_calc(part);
                                                                                                // //dt->dt_calc_ctof(part);
                                                                                                _hists->Fill_MomVsBeta(data, part, event);

                                                                                                _hists->Fill_deltat_before_cut(data, dt, part, event);

                                                                                                if (data->charge(part) > 0)
                                                                                                {
                                                                                                        _hists->FillHists_prot_pid_cuts(data, event, part);
                                                                                                        _hists->FillHists_pip_pid_cuts(data, event, part);

                                                                                                        if ((cuts->IsProton(part)) && (cuts->IsPip(part)))
                                                                                                        {

                                                                                                                dp_prot1 = (pow((mc_event->prot_momX_mc_gen() - data->px(part)), 2) +
                                                                                                                            pow((mc_event->prot_momY_mc_gen() - data->py(part)), 2) +
                                                                                                                            pow((mc_event->prot_momZ_mc_gen() - data->pz(part)), 2));

                                                                                                                dp_pip1 = (pow((mc_event->pip_momX_mc_gen() - data->px(part)), 2) +
                                                                                                                           pow((mc_event->pip_momY_mc_gen() - data->py(part)), 2) +
                                                                                                                           pow((mc_event->pip_momZ_mc_gen() - data->pz(part)), 2));

                                                                                                                // std::cout << "dp_prot1 :  " << dp_prot1 << "   dp_pip1  " << dp_pip1 << std::endl;

                                                                                                                _hists->Fill_deltaP_ambi_all_prot(event, dp_prot1);
                                                                                                                _hists->Fill_deltaP_ambi_all_pip(event, dp_pip1);
                                                                                                                // if ((dp_prot1 < 0.1) || (dp_pip1 < 0.1))
                                                                                                                {
                                                                                                                        // if (dp_prot1 < dp_pip1)
                                                                                                                        {
                                                                                                                                _hists->Fill_deltaP_ambi_prot(event, dp_prot1);

                                                                                                                                event->SetProton(part);
                                                                                                                                prot_idx++;
                                                                                                                                statusProt = abs(data->status(part));
                                                                                                                                sectorProt = data->dc_sec(part);
                                                                                                                                // std::cout << "   dt  :  " << dt << std::endl;

                                                                                                                                _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                                                                                                                                // _hists->FillHists_prot_pid_with_cuts(data, event, part, *event->GetProtons()[prot_idx]);
                                                                                                                        }
                                                                                                                        // else
                                                                                                                        {
                                                                                                                                _hists->Fill_deltaP_ambi_pip(event, dp_pip1);
                                                                                                                                event->SetPip(part);
                                                                                                                                pip_idx++;
                                                                                                                                statusPip = abs(data->status(part));
                                                                                                                                sectorPip = data->dc_sec(part);
                                                                                                                                _hists->Fill_deltat_pip_after_cut(data, dt, part, event);
                                                                                                                                // _hists->FillHists_pip_pid_with_cuts(data, event, part, *event->GetPips()[pip_idx]);
                                                                                                                        }
                                                                                                                }
                                                                                                        }

                                                                                                        else if (cuts->IsProton(part))
                                                                                                        {
                                                                                                                prot_idx++;
                                                                                                                dp_prot2 = (pow((mc_event->prot_momX_mc_gen() - data->px(part)), 2) +
                                                                                                                            pow((mc_event->prot_momY_mc_gen() - data->py(part)), 2) +
                                                                                                                            pow((mc_event->prot_momZ_mc_gen() - data->pz(part)), 2));
                                                                                                                dp_pip2 = (pow((mc_event->pip_momX_mc_gen() - data->px(part)), 2) +
                                                                                                                           pow((mc_event->pip_momY_mc_gen() - data->py(part)), 2) +
                                                                                                                           pow((mc_event->pip_momZ_mc_gen() - data->pz(part)), 2));

                                                                                                                _hists->Fill_deltaP_prot(event, dp_prot2);
                                                                                                                _hists->Fill_deltaP_pip_for_prot(event, dp_pip2);

                                                                                                                _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                                                                                                                // _hists->FillHists_prot_pid_with_cuts(data, event, part, *event->GetProtons()[prot_idx]);
                                                                                                        }

                                                                                                        else if (cuts->IsPip(part))
                                                                                                        {
                                                                                                                pip_idx++;

                                                                                                                dp_prot3 = (pow((mc_event->prot_momX_mc_gen() - data->px(part)), 2) + pow((mc_event->prot_momY_mc_gen() - data->py(part)), 2) +
                                                                                                                            pow((mc_event->prot_momZ_mc_gen() - data->pz(part)), 2));
                                                                                                                dp_pip3 = (pow((mc_event->pip_momX_mc_gen() - data->px(part)), 2) +
                                                                                                                           pow((mc_event->pip_momY_mc_gen() - data->py(part)), 2) +
                                                                                                                           pow((mc_event->pip_momZ_mc_gen() - data->pz(part)), 2));
                                                                                                                _hists->Fill_deltaP_prot_for_pip(event, dp_prot3);
                                                                                                                _hists->Fill_deltaP_pip(event, dp_pip3);

                                                                                                                _hists->Fill_deltat_pip_after_cut(data, dt, part, event);
                                                                                                                // _hists->FillHists_pip_pid_with_cuts(data, event, part, *event->GetPips()[pip_idx]);
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                        _hists->Fill_Entries(entries_in_this_event);
                                }
                        }
                }
        }
        std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
        std::cout << "   there are  " << prot_pip << "  particles identified as both prot and pip" << std::endl;
        std::cout << "   % of events having more than one proton  " << float(more_prot) / float(prot) << std::endl;
        std::cout << "   % of events having more than one pip  " << float(more_pip) / float(pip) << std::endl;
        std::cout << "   % of events having more than one pim  " << float(more_pim) / float(pim) << std::endl;
        std::cout << "   # of particles as prot as well as pip all events  " << float(prot_pip) << std::endl;
        std::cout << "   # of particles as prot as well as pip  for selected twoPi mPim events " << no_prot_pip << std::endl;
        std::cout << "   # of particles as prot as well as pip  for thrown  events " << no_prot_pip_mc << std::endl;

        std::cout.precision(3);
        std::cout
            << "no of twoPion events -->  mPim  = " << two_pion_mPim_events
            << "   Excl = " << two_pion_Excl_events << "  ratio : " << float(two_pion_Excl_events) / float(two_pion_mPim_events) << '\n'
            << " elec " << elec << " prot " << prot << " pip " << pip << " pim " << pim << '\n';
        //   << "miss_pim  " << miss_pim << " mPim ratio = " << float(miss_pim) / float(total) << '\n'
        //   << "miss_prot  " << miss_prot << " mProt ratio = " << float(miss_prot) / float(total) << '\n'
        //   << "miss_pip  " << miss_pip << " mPip ratio = " << float(miss_pip) / float(total) << '\n'
        //   //   << "other  " << other << '\n'
        //   << "excl twopion = " << excl_events << "excl ratio = " << float(excl_events) / float(total) << '\n'
        //   << "twoPions = " << twopi << '\n'
        //   << " total " << total << " total twoPi events ratio = " << twopi / total << '\n';

        // Return the total number of events

        return num_of_events;
}
#endif
