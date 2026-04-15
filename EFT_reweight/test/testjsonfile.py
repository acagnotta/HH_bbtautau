import os
import json
import re
import yaml
import numpy as np
import correctionlib.schemav2 as cs
import optparse

usage                   = 'python3 HEFT_generalizedweight.py '
parser                  = optparse.OptionParser(usage)
parser.add_option('-c', '--config', dest='config', type=str, default="../config/config.yaml", help='Path to the config file')
parser.add_option('--dryrun', dest='dryrun', action='store_true', default=False, help='If true, the script will not save the correctionlib file, but will just print the poly values for the target samples and exit')
(opt, args)             = parser.parse_args()

config_file = opt.config
dryrun = opt.dryrun

def read_coeff_binning_from_json(filepath):
    with open(filepath, 'r') as f:
        data = json.load(f)
    binning_dicts = {}
    for key in data.keys(): # pthh_0.0-20_theta_0.0-0.25
        elements = key.split('_')
        var1 = elements[0] # pthh
        bin_range1 = elements[1] # 0.0-20
        var2 = elements[2] # theta
        bin_range2 = elements[3] # 0.0-0.25
        if not var1 in binning_dicts:
            binning_dicts[var1] = []
        if not var2 in binning_dicts:
            binning_dicts[var2] = []
        binning_dicts[var1].append(list(map(float, bin_range1.split('-'))))
        binning_dicts[var2].append(list(map(float, bin_range2.split('-'))))
        for key2 in data[key].keys(): #250-270
            if not 'mhh' in binning_dicts:
                binning_dicts['mhh'] = []
            bin_range_mhh = key2
            binning_dicts['mhh'].append(list(map(float, bin_range_mhh.split('-'))))
    # Now sort and make unique
    for var in binning_dicts.keys():
        bins = sorted(set([item for sublist in binning_dicts[var] for item in sublist]))
        binning_dicts[var] = bins
    return binning_dicts, data

def order_coeffs(coeffs_dict, binning, key_ = "fitted_parameters"):
    nbin_var1 = len(binning['pthh']) - 1
    nbin_var2 = len(binning['theta']) - 1
    nbin_mhh = len(binning['mhh']) - 1
    ncoeffs = 23
    coeffs_array = np.zeros((nbin_var1, nbin_var2, nbin_mhh, ncoeffs)) if key_=="fitted_parameters" else np.zeros((nbin_var1, nbin_var2, nbin_mhh, ncoeffs, ncoeffs))
    for key in coeffs_dict.keys(): # pthh_0.0-20_theta_0.0-0.25
        elements = key.split('_')
        var1 = elements[0] # pthh
        bin_range1 = list(map(float, elements[1].split('-'))) # [0.0,20.0]
        var2 = elements[2] # theta
        bin_range2 = list(map(float, elements[3].split('-'))) # [0.0,0.25]
        idx_var1 = binning[var1].index(bin_range1[0])
        idx_var2 = binning[var2].index(bin_range2[0])
        for key2 in coeffs_dict[key].keys(): #250-270
            bin_range_mhh = list(map(float, key2.split('-')))
            idx_mhh = binning['mhh'].index(bin_range_mhh[0])
            if key_=="fitted_parameters":
                coeffs_line = list(map(float, coeffs_dict[key][key2][key_]))
                coeffs_array[idx_var1, idx_var2, idx_mhh, :] = coeffs_line
            elif key_=="covariance":
                cov_matrix = np.array(coeffs_dict[key][key2][key_])
                coeffs_array[idx_var1, idx_var2, idx_mhh, :, :] = cov_matrix
    return coeffs_array

def create_correctionlib_file(coeffs_HEFT, cov_HEFT, binning):
    formulas = [
        "pow(x[1],4)", #1
        "pow(x[2],2)", #2
        "(pow(x[1],2))*(pow(x[0],2))", #3
        "(pow(x[3],2))*(pow(x[0],2))", #4
        "pow(x[4],2)", #5
        "x[2]*(pow(x[1],2))", #6
        "(pow(x[1],3))*x[0]", #7
        "x[1]*x[0]*x[2]", #8
        "x[3]*x[0]*x[2]", #9
        "x[2]*x[4]", #10
        "x[3]*x[0]*(pow(x[1],2))", #11
        "x[4]*(pow(x[1],2))", #12
        "(pow(x[0],2))*x[3]*x[1]", #13
        "x[4]*x[1]*x[0]", #14
        "x[3]*x[4]*x[0]", #15
        "(pow(x[1],3))*x[3]", #16
        "x[1]*x[2]*x[3]", #17
        "x[1]*(pow(x[3],2))*x[0]", #18
        "x[1]*x[3]*x[4]", #19
        "(pow(x[1],2))*(pow(x[3],2))", #20
        "x[2]*(pow(x[3],2))", #21
        "(pow(x[3],3))*x[0]", #22
        "(pow(x[3],2))*x[4]" #23
    ]
    formulas_ = [f.replace("x[0]", "x[5]").replace("x[1]", "x[6]").replace("x[2]", "x[7]").replace("x[3]", "x[8]").replace("x[4]", "x[9]") for f in formulas]
    poly_expr = " + ".join([f"[{i}]*({form})" for i, form in enumerate(formulas)])


    poly_err_terms = [f"({fi}*[{i*len(formulas)+j}]*{fj})" for i,fi in enumerate(formulas) for j,fj in enumerate(formulas)]
    poly_err = " + ".join(poly_err_terms)

    poly_cov_terms = [f"({fi}*[{i*len(formulas)+j}]*{fj})" for i,fi in enumerate(formulas) for j,fj in enumerate(formulas_)]
    poly_cov = " + ".join(poly_cov_terms)

    corr = []
    cparams = cs.Correction(
        name="couplings",
        version="1.0",
        description="tools for poly derivation",
        inputs=[
            cs.Variable(name="kl", type="real", description="Higgs self-coupling modifier"),
            cs.Variable(name="kt", type="real", description="Higgs top-quark coupling modifier"),
            cs.Variable(name="c2", type="real", description="Higgs gluon coupling modifier"),
            cs.Variable(name="cg", type="real", description="Higgs photon coupling modifier"),
            cs.Variable(name="c2g", type="real", description="Higgs gluon-photon coupling modifier"),
            cs.Variable(name="idx", type="int", description="which entry of the vector")
        ],
        output=cs.Variable(name="x", type="real", description="entry of the vector to multiply by A param"),
        data=cs.Category(
            nodetype="category",
            input="idx",
            content=[
                {
                    "key": i,
                    "value": cs.Formula(
                                nodetype="formula",
                                parser="TFormula",
                                variables=["kl", "kt", "c2", "cg", "c2g"],
                                parameters=[],
                                expression=f"{f}",
                            )
                } for i, f in enumerate(formulas)
            ]
        )

    )
    polys = cs.Correction(
        name="HEFT_poly",
        version="1.0",
        description="polynomial value for for HEFT reweighting",
        inputs=[
            cs.Variable(name="pthh", type="real", description="transverse momentum of the Higgs boson pair LHE level"),
            cs.Variable(name="costhetastar", type="real", description="absolute value of cos Theta* LHE level"),
            cs.Variable(name="mhh", type="real", description="invariant mass of the Higgs boson pair LHE level"),
            cs.Variable(name="kl", type="real", description="Higgs self-coupling modifier"),
            cs.Variable(name="kt", type="real", description="Higgs top-quark coupling modifier"),
            cs.Variable(name="c2", type="real", description="Higgs gluon coupling modifier"),
            cs.Variable(name="cg", type="real", description="Higgs photon coupling modifier"),
            cs.Variable(name="c2g", type="real", description="Higgs gluon-photon coupling modifier"),
        ],
        output=cs.Variable(name="poly", type="real", description="poly value"),
        data=cs.Binning(
            nodetype="binning",
            input="pthh",
            edges=binning['pthh'],
            content=[
                cs.Binning(
                    nodetype="binning",
                    input="costhetastar",
                    edges=binning['theta'],
                    content=[
                        cs.Binning(
                            nodetype="binning",
                            input="mhh",
                            edges=binning['mhh'],
                            content=[
                                cs.Formula(
                                    nodetype="formula",
                                    parser="TFormula",
                                    variables=["kl", "kt", "c2", "cg", "c2g"],
                                    parameters=[n for n in coeffs_HEFT[i, j, k, :].tolist()],
                                    expression=f"{poly_expr}",
                                )
                                for k in range(len(binning["mhh"]) - 1)
                            ],
                            flow="clamp"
                        ) for j in range(len(binning["theta"]) - 1)
                    ],
                    flow="clamp"
                ) for i in range(len(binning["pthh"]) - 1)
            ],
            flow="clamp"
        )
    )

    polys_err = cs.Correction(
        name="HEFT_poly_error",
        version="1.0",
        description="sigma2 on poly value for HEFT reweighting",
        inputs=[
            cs.Variable(name="pthh", type="real", description="transverse momentum of the Higgs boson pair LHE level"),
            cs.Variable(name="costhetastar", type="real", description="absolute value of cos Theta* LHE level"),
            cs.Variable(name="mhh", type="real", description="invariant mass of the Higgs boson pair LHE level"),
            cs.Variable(name="kl", type="real", description="Higgs self-coupling modifier"),
            cs.Variable(name="kt", type="real", description="Higgs top-quark coupling modifier"),
            cs.Variable(name="c2", type="real", description="Higgs gluon coupling modifier"),
            cs.Variable(name="cg", type="real", description="Higgs photon coupling modifier"),
            cs.Variable(name="c2g", type="real", description="Higgs gluon-photon coupling modifier"),
        ],
        output=cs.Variable(name="poly", type="real", description="poly value"),
        data=cs.Binning(
            nodetype="binning",
            input="pthh",
            edges=binning['pthh'],
            content=[
                cs.Binning(
                    nodetype="binning",
                    input="costhetastar",
                    edges=binning['theta'],
                    content=[
                        cs.Binning(
                            nodetype="binning",
                            input="mhh",
                            edges=binning['mhh'],
                            content=[
                                cs.Formula(
                                    nodetype="formula",
                                    parser="TFormula",
                                    variables=["kl", "kt", "c2", "cg", "c2g"],
                                    parameters=[n for n in cov_HEFT[i, j, k, :, :].flatten().tolist()],
                                    expression=f"{poly_err}",
                                )
                                for k in range(len(binning["mhh"]) - 1)
                            ],
                            flow="clamp"
                        ) for j in range(len(binning["theta"]) - 1)
                    ],
                    flow="clamp"
                ) for i in range(len(binning["pthh"]) - 1)
            ],
            flow="clamp"
        )
    )

    cov = cs.Correction(
        name="poly_covariance",
        version="1.0",
        description="covariance between two polys for HEFT reweighting",
        inputs=[
            cs.Variable(name="pthh", type="real", description="transverse momentum of the Higgs boson pair LHE level"),
            cs.Variable(name="costhetastar", type="real", description="absolute value of cos Theta* LHE level"),
            cs.Variable(name="mhh", type="real", description="invariant mass of the Higgs boson pair LHE level"),
            cs.Variable(name="kl_i", type="real", description="Higgs self-coupling modifier input sample"),
            cs.Variable(name="kt_i", type="real", description="Higgs top-quark coupling modifier input sample"),
            cs.Variable(name="c2_i", type="real", description="Higgs gluon coupling modifier input sample"),
            cs.Variable(name="cg_i", type="real", description="Higgs photon coupling modifier input sample"),
            cs.Variable(name="c2g_i", type="real", description="Higgs gluon-photon coupling modifier input sample"),
            cs.Variable(name="kl_o", type="real", description="Higgs self-coupling modifier output sample"),
            cs.Variable(name="kt_o", type="real", description="Higgs top-quark coupling modifier output sample"),
            cs.Variable(name="c2_o", type="real", description="Higgs gluon coupling modifier output sample"),
            cs.Variable(name="cg_o", type="real", description="Higgs photon coupling modifier output sample"),
            cs.Variable(name="c2g_o", type="real", description="Higgs gluon-photon coupling modifier output sample"),
        ],
        output=cs.Variable(name="params", type="real", description="parameters of the polynomial fit for HEFT reweighting"),
        data=cs.Binning(
            nodetype="binning",
            input="pthh",
            edges=binning['pthh'],
            content=[
                cs.Binning(
                    nodetype="binning",
                    input="costhetastar",
                    edges=binning['theta'],
                    content=[
                        cs.Binning(
                            nodetype="binning",
                            input="mhh",
                            edges=binning['mhh'],
                            content=[
                                   cs.Formula(
                                    nodetype="formula",
                                    parser="TFormula",
                                    variables=["kl_i", "kt_i", "c2_i", "cg_i", "c2g_i", "kl_o", "kt_o", "c2_o", "cg_o", "c2g_o"],
                                    parameters=[n for n in cov_HEFT[i, j, k, :, :].flatten().tolist()],
                                    expression=f"{poly_cov}",
                                ) for k in range(len(binning['mhh']) - 1)],
                            flow="clamp"
                            ) for j in range(cov_HEFT.shape[1])
                        ],
                        flow="clamp"
                    ) for i in range(cov_HEFT.shape[0])
                ],
                flow="clamp"
            ),                      
        )
    corr.append(cparams)
    corr.append(polys)
    cset = cs.CorrectionSet(schema_version=2, corrections=corr)
    with open(f"poly.json", "w") as f:
        json.dump(cset.model_dump(), f, indent=2)
    corr=[]
    corr.append(cov)
    cset = cs.CorrectionSet(schema_version=2, corrections=corr)
    with open(f"cov.json", "w") as f:
        json.dump(cset.model_dump(), f, indent=2)
    corr = []
    corr.append(polys_err)
    cset = cs.CorrectionSet(schema_version=2, corrections=corr)
    with open(f"err.json", "w") as f:
        json.dump(cset.model_dump(), f, indent=2)

    return 0



if __name__ == "__main__":
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    main_dir = os.path.dirname(os.path.dirname(config_file))

    input_file_coeff = os.path.join(main_dir, config['2D_HEFT']['coeff'])
    input_cov_file = os.path.join(main_dir, config['2D_HEFT']['cov_matrix'])

    binning, coeffs_HEFT = read_coeff_binning_from_json(input_file_coeff)
    binning_, cov_HEFT = read_coeff_binning_from_json(input_cov_file)

    coeffs_HEFT = order_coeffs(coeffs_HEFT, binning, key_="fitted_parameters")
    cov_HEFT = order_coeffs(cov_HEFT, binning_, key_="covariance")
    create_correctionlib_file(coeffs_HEFT, cov_HEFT, binning)