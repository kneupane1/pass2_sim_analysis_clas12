

// energy loss corr
int _prot_status = NAN;
int _pim_status = NAN;

double _px_prime_prot_E = NAN;
double _py_prime_prot_E = NAN;
double _pz_prime_prot_E = NAN;

double _prot_mom_uncorr = NAN;
double _prot_theta_uncorr = NAN;
float _prot_phi_uncorr = NAN;
double _prot_mom_tmt = NAN;

double _px_prime_pim_E = NAN;
double _py_prime_pim_E = NAN;
double _pz_prime_pim_E = NAN;

double _pim_mom_tmt = NAN;
double _pim_mom_uncorr = NAN;
double _pim_theta_uncorr = NAN;
double _pim_phi_uncorr = NAN;

void Reaction::SetProton(int i)
{
    auto _prot = std::make_unique<TLorentzVector>();
    auto _Energy_loss_uncorr_prot = std::make_unique<TLorentzVector>();

    _prot_status = abs(_data->status(i));

    _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

    _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();
    _prot_theta_uncorr = _Energy_loss_uncorr_prot->Theta() * 180 / PI;
    if (_Energy_loss_uncorr_prot->Phi() > 0)
        _prot_phi_uncorr = _Energy_loss_uncorr_prot->Phi() * 180 / PI;
    else if (_Energy_loss_uncorr_prot->Phi() < 0)
        _prot_phi_uncorr = (_Energy_loss_uncorr_prot->Phi() + 2 * PI) * 180 / PI;

    if (_prot_status >= 4000)
    {
        _prot_mom_tmt = _prot_mom_uncorr;
    }
    if (_prot_status < 4000)
    {
        if (_prot_theta_uncorr < 27)
        {
            if (_prot_mom_tmt < 2.4)
                _prot_mom_tmt = _prot_mom_uncorr + (0.001046) * pow(_prot_mom_uncorr, 4) +
                                (-0.010446) * pow(_prot_mom_uncorr, 3) + (0.036945) * pow(_prot_mom_uncorr, 2) +
                                (-0.055368) * _prot_mom_uncorr + 0.034539;
            else
                _prot_mom_tmt = _prot_mom_uncorr + 0.004741;
        }
        else
        {
            if (_prot_mom_tmt < 2.4)
                _prot_mom_tmt = _prot_mom_uncorr + (0.005519) * pow(_prot_mom_uncorr, 4) +
                                (-0.046289) * pow(_prot_mom_uncorr, 3) + (0.137504) * pow(_prot_mom_uncorr, 2) +
                                (-0.177027) * _prot_mom_uncorr + 0.094555;
            else
                _prot_mom_tmt = _prot_mom_uncorr + 0.004899;
        }
    }

    _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
    _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
    _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

    _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P); // energy loss corrected
}

void Reaction::SetPim(int i)
{

    auto _pim = std::make_unique<TLorentzVector>();
    auto _Energy_loss_uncorr_pim = std::make_unique<TLorentzVector>();

    _pim_status = abs(_data->status(i));

    _Energy_loss_uncorr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);

    _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
    _pim_theta_uncorr = _Energy_loss_uncorr_pim->Theta() * 180 / PI;
    if (_Energy_loss_uncorr_pim->Phi() > 0)
        _pim_phi_uncorr = _Energy_loss_uncorr_pim->Phi() * 180 / PI;
    else if (_Energy_loss_uncorr_pim->Phi() < 0)
        _pim_phi_uncorr = (_Energy_loss_uncorr_pim->Phi() + 2 * PI) * 180 / PI;

    if (_pim_status >= 4000)
    {
        _pim_mom_tmt = _pim_mom_uncorr;
    }
    if (_pim_status < 4000)
    {
        if (_pim_theta_uncorr < 27)
        {
            _pim_mom_tmt = _pim_mom_uncorr + 0.00046571 * _pim_mom_uncorr + 0.00322164;
        }
        else
        {
            if (_pim_mom_uncorr < 1.7)
                _pim_mom_tmt = _pim_mom_uncorr + (-0.0024313) * pow(_pim_mom_uncorr, 3) +
                               (0.0094416) * pow(_pim_mom_uncorr, 2) + (-0.01257967) * pow(_pim_mom_uncorr, 1) + 0.0122432;
            else
                _pim_mom_tmt = _pim_mom_uncorr + 0.006199071;
        }
    }
    _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
    _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
    _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

    _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM); // energy loss corrected
}
