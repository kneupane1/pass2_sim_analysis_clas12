
#include "physics.hpp"

namespace physics {
// double IsMomCorrected(double pe, double px, double py){
//
//         //Correction function: Pnew = Pe - dP;
//
//         double epars[] = {-0.07448192048426636, 0.007612567755187672, 0.002418059256619887, 5.228934861536853e-05, -0.00014492558875767224, 1.993876401764548e-05, -0.02064213056508147, 0.0029545131471652403, 0.00027791314405979585, 0.00010380511562495014, -0.00010469092452846252, 1.8897742939189947e-05, -0.08588081509553264, 0.015878617692830884, -0.00033117742120785244, -0.00031359008739008707, -2.0350707393024823e-05, 9.728030236830449e-06, -0.05629662730516873, 0.00684233309175524, -0.002089924807259924, 1.7596057592940903e-05, 9.729280093911233e-05, -1.4531945568558965e-05, -0.03749122961267899, 0.005419067236090319, -0.0030724195868331856, 0.0005080717477174384, 0.00010349899351904892, -2.0918174366357658e-05, -0.07210983308002632, 0.00966169711164072, 0.001448822048331462, -4.9070424768743856e-05, -4.371204991147332e-05, 6.286676815939829e-06};
//
//         double p0 = epars[isec*6] + epars[isec*6+1]*pe;
//
//         double p1 = epars[isec*6+2] + epars[isec*6+3]*pe;
//
//         double p2 = epars[isec*6+4] + epars[isec*6+5]*pe;
//
//         double fi = fi = TMath::RadToDeg()*atan2(py,px);
//
//         double fie = fi + (fi<0 && isec>1)*360 - isec*60;
//
//         double dP = p0 + p1*fie + p2*fie*fie;
//
//         pe = pe -dP;
//
//         return dP;
// }


// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
        TLorentzVector q_mu = (e_mu - e_mu_prime);
        return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
        TLorentzVector q_mu = (e_mu - e_mu_prime);
        TVector3 p_mu_3(0, 0, 0);
        TLorentzVector p_mu;
        p_mu.SetVectM(p_mu_3, MASS_P);
        return ((p_mu + q_mu).Mag());
}

// overload with 4 vectors
double xb_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
        double Q2 = Q2_calc(e_mu, e_mu_prime);
        TLorentzVector q = e_mu - e_mu_prime;
        TLorentzVector target(0, 0, 0, MASS_P);
        return (Q2 / (2 * (q.Dot(target))));
}
double vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta) {
        return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

double deltat(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r) {
        double relatavistic_beta = 1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
        return electron_vertex_time - vertex_time(sc_t, sc_r, relatavistic_beta);
}
}  // namespace physics
