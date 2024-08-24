#ifndef MOM_CORR_H_GUARD
#define MOM_CORR_H_GUARD
#include <cstdlib>
#include <ctime> // add this line to include <ctime> header
#include <iostream>
#include "TROOT.h"
#include "constants.hpp"
// using namespace std;
class mom_corr
{
private:
public:
  mom_corr(){};
  ~mom_corr();

  bool is_FD(int prot_status);
  bool is_CD(int prot_status);

  // hadron mom corrections
  double dppC(float Px, float Py, float Pz, int sec, int ivec);
  double elossPipFD(double pion_p, double pip_theta);
  double elossPipCD(double pion_p, double pip_theta);

  // 4-vector method dp corrections:
  float CD_prot_Hmom_corr(float mom_, float phi_, float alpha_prot);
  float FD_prot_Hmom_corr(float mom_, float dc_sec, float alpha_prot);

  float CD_pip_Hmom_corr(float mom_, float phi_, float alpha_pip);
  float FD_pip_Hmom_corr(float mom_, float dc_sec, float alpha_pip);

  float CD_pim_Hmom_corr(float mom_, float phi_, float alpha_pim);
  float FD_pim_Hmom_corr(float mom_, float dc_sec, float alpha_pim);

  // float alpha_prot_mom_corr_FD[4] = {0.5, 0.6, 0.5, 0.5};

  // void random_no_gen() {
  //   std::srand(std::time(nullptr));  // seed the random number generator
  //                                    // rest of your code here
  //   for (int i = 0; i < 4; i++) {
  //     float variation =
  //         -0.1 + (static_cast<float>(std::rand()) / RAND_MAX) * 0.2;  // generate random percentage between -10%
  //         and +10%
  //     alpha_prot_mom_corr_FD[i] *= (1.0 + variation);  // apply the random percentage to the current element

  //     // std::cout << "New value in element " << i << " of array: " << alpha_prot_mom_corr_FD[i] << std::endl;
  //   }
  // }
};

#endif
