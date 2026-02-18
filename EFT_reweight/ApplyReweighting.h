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

float GetInputXSec(const std::string& sample_name) {{
    auto it = xSec.find(sample_name);
    if (it != xSec.end()) {{
        return it->second;
    }}
    return 1.0;
}}

float GetNtotalEvents(const std::string& sample_name) {{
    auto it = NtotalEvents.find(sample_name);
    if (it != NtotalEvents.end()) {{
        return it->second;
    }}
    return -1;
}}

float GetMhhGen(rvec_f GenPart_pt, rvec_f GenPart_eta, rvec_f GenPart_phi, rvec_f GenPart_mass, rvec_i GenPart_statusFlags, rvec_i GenPart_pdgId){
    float mhhGen = 0.;
    auto dihiggs_p4 = GetDiHiggsP4Gen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId);
    mhhGen = dihiggs_p4.M();
    return mhhGen;
}

float GetMhhLHE(rvec_f LHEPart_pt, rvec_f LHEPart_eta, rvec_f LHEPart_phi, rvec_f LHEPart_mass, rvec_i LHEPart_pdgId){
    float mhhGen = 0.;
    auto dihiggs_p4 = GetDiHiggsP4LHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId);
    mhhGen = dihiggs_p4.M();
    return mhhGen;
}

float GetPthhGen(rvec_f GenPart_pt, rvec_f GenPart_eta, rvec_f GenPart_phi, rvec_f GenPart_mass, rvec_i GenPart_statusFlags, rvec_i GenPart_pdgId){
    float pthhGen = 0.;
    auto dihiggs_p4 = GetDiHiggsP4Gen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId);
    pthhGen = dihiggs_p4.Pt();
    return pthhGen;
}

float GetPthhLHE(rvec_f LHEPart_pt, rvec_f LHEPart_eta, rvec_f LHEPart_phi, rvec_f LHEPart_mass, rvec_i LHEPart_pdgId){
    float pthhGen = 0.;
    auto dihiggs_p4 = GetDiHiggsP4LHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId);
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

float GetWeight(float mhh_gen, float pthh_gen, float costhetastar_gen, const std::string& sampleName, const std::string& file, const std::string& target){
    auto cset = CorrectionSet::from_file(file);
    auto weights = cset->at("HEFT_reweighting_" + target);
    float weight = 1.0;
    if (pthh_gen >= 0 && std::abs(costhetastar_gen) >= 0 && mhh_gen >= 250) {
        weight = weights->evaluate({sampleName, pthh_gen, std::abs(costhetastar_gen), mhh_gen});
    }
    return weight;
}

float GetWeightFromPoly(float mhh_gen, float pthh_gen, float costhetastar_gen, const std::string& sampleName, const std::string& file, const std::string& target){
    auto cset = CorrectionSet::from_file(file);
    auto poly_target = cset->at("HEFT_poly_" + target);
    auto poly_input = cset->at("HEFT_poly_" + sampleName);
    float weight = 1.0;
    if (pthh_gen >= 0 && std::abs(costhetastar_gen) >= 0 && mhh_gen >= 250) {
        weight = poly_target->evaluate({pthh_gen, std::abs(costhetastar_gen), mhh_gen}) / poly_input->evaluate({pthh_gen, std::abs(costhetastar_gen), mhh_gen});
    }
    return weight;
}

float GetWeightFromPolyPDF(float mhh_gen, float pthh_gen, float costhetastar_gen, const std::string& sampleName, const std::string& file, const std::string& target, const std::string& filePDF){
    auto cset = CorrectionSet::from_file(file);
    auto cset_pdf = CorrectionSet::from_file(filePDF);
    auto poly_target = cset->at("HEFT_poly_" + target);
    auto poly_input = cset->at("HEFT_poly_" + sampleName);
    auto pdf_input = cset_pdf->at("PDF_weighted_SM");
    float weight = 1.0;
    if (pthh_gen >= 0 && std::abs(costhetastar_gen) >= 0 && mhh_gen >= 250) {
        weight = pdf_input->evaluate({pthh_gen, std::abs(costhetastar_gen), mhh_gen}) * (poly_target->evaluate({pthh_gen, std::abs(costhetastar_gen), mhh_gen}) / poly_input->evaluate({pthh_gen, std::abs(costhetastar_gen), mhh_gen}));
    }
    return weight;
}