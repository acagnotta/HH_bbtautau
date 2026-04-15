#include "correction.h"

using correction::CorrectionSet;
using ROOT::VecOps::RVec;

using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;
using rvec_b = const RVec<bool> &;
using rvec_rvec_i = const RVec<RVec<int>> &;
using rvec_rvec_f = const RVec<RVec<float>> &;

TLorentzVector GetDiHiggsP4LHE(rvec_f LHEPart_pt, rvec_f LHEPart_eta, rvec_f LHEPart_phi, rvec_f LHEPart_mass, rvec_i LHEPart_pdgId){
    TLorentzVector dihiggs_p4(0.,0.,0.,0.);
    RVec<TLorentzVector> higgs_p4s;
    TLorentzVector tmp_vec(0.,0.,0.,0.);
    for (size_t i = 0; i < LHEPart_pt.size(); ++i) {
        tmp_vec.SetPtEtaPhiM(0.,0.,0.,0.);
        if (abs(LHEPart_pdgId[i]) == 25) {
            tmp_vec.SetPtEtaPhiM(LHEPart_pt[i], LHEPart_eta[i], LHEPart_phi[i], LHEPart_mass[i]);
            higgs_p4s.push_back(tmp_vec);
        }
    }
    if (higgs_p4s.size() == 2) {
        dihiggs_p4 = higgs_p4s[0] + higgs_p4s[1];
    } else {
        cout << "Warning: Found " << higgs_p4s.size() << " Higgs bosons instead of 2." << endl;
        dihiggs_p4 = TLorentzVector(0.,0.,0.,0.);
    }
    return dihiggs_p4;
}

float GetMhhLHE(rvec_f LHEPart_pt, rvec_f LHEPart_eta, rvec_f LHEPart_phi, rvec_f LHEPart_mass, rvec_i LHEPart_pdgId){
    float mhhGen = 0.;
    auto dihiggs_p4 = GetDiHiggsP4LHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId);
    mhhGen = dihiggs_p4.M();
    return mhhGen;
}

float GetPthhLHE(rvec_f LHEPart_pt, rvec_f LHEPart_eta, rvec_f LHEPart_phi, rvec_f LHEPart_mass, rvec_i LHEPart_pdgId){
    float pthhGen = 0.;
    auto dihiggs_p4 = GetDiHiggsP4LHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId);
    pthhGen = dihiggs_p4.Pt();
    return pthhGen;
}


float GetCosThetaStarLHE(rvec_f LHEPart_pt, rvec_f LHEPart_eta, rvec_f LHEPart_phi, rvec_f LHEPart_mass, rvec_i LHEPart_pdgId){
    float cosThetaStarGen = 0.;
    RVec<TLorentzVector> higgs_p4s;
    TLorentzVector tmp_vec(0.,0.,0.,0.);
    for (size_t i = 0; i < LHEPart_pt.size(); ++i) {
        tmp_vec.SetPtEtaPhiM(0.,0.,0.,0.);
        if (abs(LHEPart_pdgId[i]) == 25) {
            tmp_vec.SetPtEtaPhiM(LHEPart_pt[i], LHEPart_eta[i], LHEPart_phi[i], LHEPart_mass[i]);
            higgs_p4s.push_back(tmp_vec);
        }
    }
    if (higgs_p4s.size() == 2) {
        auto dihiggs_p4 = higgs_p4s[0] + higgs_p4s[1];
        auto boost_vector = dihiggs_p4.BoostVector();
        auto higgs1_boosted = higgs_p4s[0];
        higgs1_boosted.Boost(-boost_vector);
        cosThetaStarGen = std::cos(higgs1_boosted.Theta());
    } else {
        cosThetaStarGen = -2.; // Invalid value
    }
    return cosThetaStarGen;
}

float GetHEFTweight(float mhh_lhe, float pthh_lhe, float costhetastar_lhe, float kl_i, float kt_i, float c2_i, float cg_i, float c2g_i, float kl_o, float kt_o, float c2_o, float cg_o, float c2g_o, const std::string& file){
    static auto cset = CorrectionSet::from_file(file);
    static auto poly = cset->at("HEFT_poly");
    float weight = 1.0;
    if (pthh_lhe >= 0 && std::abs(costhetastar_lhe) >= 0 && mhh_lhe >= 250) {
        weight = poly->evaluate({pthh_lhe, std::abs(costhetastar_lhe), mhh_lhe, kl_i, kt_i, c2_i, cg_i, c2g_i}) / poly->evaluate({pthh_lhe, std::abs(costhetastar_lhe), mhh_lhe, kl_o, kt_o, c2g_o, cg_o, c2g_o});
    }
    return weight;
}

float GetHEFTweightError(float mhh_lhe, float pthh_lhe, float costhetastar_lhe, float kl_i, float kt_i, float c2_i, float cg_i, float c2g_i, float kl_o, float kt_o, float c2_o, float cg_o, float c2g_o, const std::string& file_err_poly, const std::string& file_val_poly, const std::string& file_cov_poly){
    static auto csetpoly = CorrectionSet::from_file(file_val_poly);
    static auto csetcov = CorrectionSet::from_file(file_cov_poly);
    static auto cseterr = CorrectionSet::from_file(file_err_poly);
    static auto poly_error = cseterr->at("HEFT_poly_error");
    static auto poly = csetpoly->at("HEFT_poly");
    static auto cov = csetcov->at("HEFT_poly_params_covariance");

    float weight = 1.0;
    int num_params = 23;

    if (pthh_lhe >= 0 && std::abs(costhetastar_lhe) >= 0 && mhh_lhe >= 250) {
        float val_target = poly->evaluate({pthh_lhe, std::abs(costhetastar_lhe), mhh_lhe, kl_i, kt_i, c2_i, cg_i, c2g_i});
        float val_input = poly->evaluate({pthh_lhe, std::abs(costhetastar_lhe), mhh_lhe, kl_o, kt_o, c2_o, cg_o, c2g_o});
        float err_input = poly_error->evaluate({pthh_lhe, std::abs(costhetastar_lhe), mhh_lhe, kl_i, kt_i, c2_i, cg_i, c2g_i});
        float err_target = poly_error->evaluate({pthh_lhe, std::abs(costhetastar_lhe), mhh_lhe, kl_o, kt_o, c2_o, cg_o, c2g_o});

        float cab = cov->evaluate({pthh_lhe, std::abs(costhetastar_lhe), mhh_lhe, kl_i, kt_i, c2_i, cg_i, c2g_i, kl_o, kt_o, c2_o, cg_o, c2g_o});
        // Error propagation for ratio: (a/b) -> sqrt((da/b)^2 + (a*db/b^2)^2 - 2*a*cov/(b^3))) 
        weight = std::sqrt(std::pow(err_target / val_input, 2) + std::pow(val_target * err_input / (val_input * val_input), 2) - 2 * val_target * cab / (val_input * val_input * val_input));
    }
    return weight;
}