#include "correction.h"

using correction::CorrectionSet;
using ROOT::VecOps::RVec;

using RNode = ROOT::RDF::RNode;
using rvec_f = const RVec<float> &;
using rvec_i = const RVec<int> &;
using rvec_b = const RVec<bool> &;
using rvec_rvec_i = const RVec<RVec<int>> &;
// using TLorentzVector = const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiMVector<float>>;

TLorentzVector GetDiHiggsP4Gen(rvec_f GenPart_pt, rvec_f GenPart_eta, rvec_f GenPart_phi, rvec_f GenPart_mass, rvec_i GenPart_statusFlags, rvec_i GenPart_pdgId){
    TLorentzVector dihiggs_p4(0.,0.,0.,0.);
    RVec<TLorentzVector> higgs_p4s;
    TLorentzVector tmp_vec(0.,0.,0.,0.);
    for (size_t i = 0; i < GenPart_pt.size(); ++i) {
        tmp_vec.SetPtEtaPhiM(0.,0.,0.,0.);
        if (abs(GenPart_pdgId[i]) == 25 && (GenPart_statusFlags[i] & 1<<13)) {
            tmp_vec.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
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

float GetMhhGen(rvec_f GenPart_pt, rvec_f GenPart_eta, rvec_f GenPart_phi, rvec_f GenPart_mass, rvec_i GenPart_statusFlags, rvec_i GenPart_pdgId){
    float mhhGen = 0.;
    auto dihiggs_p4 = GetDiHiggsP4Gen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId);
    mhhGen = dihiggs_p4.M();
    return mhhGen;
}

float GetPthhGen(rvec_f GenPart_pt, rvec_f GenPart_eta, rvec_f GenPart_phi, rvec_f GenPart_mass, rvec_i GenPart_statusFlags, rvec_i GenPart_pdgId){
    float pthhGen = 0.;
    auto dihiggs_p4 = GetDiHiggsP4Gen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId);
    pthhGen = dihiggs_p4.Pt();
    return pthhGen;
}

float GetCosThetaStarGen(rvec_f GenPart_pt, rvec_f GenPart_eta, rvec_f GenPart_phi, rvec_f GenPart_mass, rvec_i GenPart_statusFlags, rvec_i GenPart_pdgId){
    float cosThetaStarGen = 0.;
    RVec<TLorentzVector> higgs_p4s;
    TLorentzVector tmp_vec(0.,0.,0.,0.);
    for (size_t i = 0; i < GenPart_pt.size(); ++i) {
        tmp_vec.SetPtEtaPhiM(0.,0.,0.,0.);
        if (abs(GenPart_pdgId[i]) == 25 && (GenPart_statusFlags[i] & 1<<13)) {
            tmp_vec.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
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

float GetWeight(float mhh_gen, float pthh_gen, float costhetastar_gen, const std::string& file){
    auto cset = CorrectionSet::from_file(file);
    auto weights = cset->at("HEFT_reweighting");
    float weight = 1.0;
    if (pthh_gen >= 0 && std::abs(costhetastar_gen) >= 0 && mhh_gen >= 250 && mhh_gen <= 1000) {
        weight = weights->evaluate({pthh_gen, std::abs(costhetastar_gen), mhh_gen});
    }
    return weight;
}