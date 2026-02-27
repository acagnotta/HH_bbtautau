import os
import json
import re
import yaml
import numpy as np
import correctionlib.schemav2 as cs
import optparse

usage                   = 'python3 HEFT_generalizedweight.py '
parser                  = optparse.OptionParser(usage)
parser.add_option('-c', '--config', dest='config', type=str, default="./config/config.yaml", help='Path to the config file')
parser.add_option('--poly_ratio_output', dest='poly_ratio_output', action='store_true', default=False, help='If true, the script will save poly ratios in the correctionlib file instead to save polys')
parser.add_option('-o', '--output', dest='output', type=str, default="./HEFT_reweighting.json", help='Path to the output correctionlib file')
parser.add_option('--dryrun', dest='dryrun', action='store_true', default=False, help='If true, the script will not save the correctionlib file, but will just print the poly values for the target samples and exit')
(opt, args)             = parser.parse_args()

config_file = opt.config
output_name = opt.output
save_poly_ratio = opt.poly_ratio_output
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

def save_max_error_from_cov(cov_HEFT, binning):
    ncoeffs = 23
    max_errors = np.zeros((len(binning['pthh']) - 1, len(binning['theta']) - 1, len(binning['mhh']) - 1, ncoeffs))
    for i in range(cov_HEFT.shape[0]):
        for j in range(cov_HEFT.shape[1]):
            for k in range(cov_HEFT.shape[2]):
                cov_matrix = cov_HEFT[i, j, k, :, :]
                errors = np.diag(cov_matrix)
                max_errors[i, j, k, :] = errors
    return max_errors

def prepare_params_for_poly(params):
    x = np.array([params['kt']**4, #1
                params['c2']**2, #2
                (params['kt']**2)*(params['kl']**2), #3
                (params['cg']**2)*(params['kl']**2), #4
                params['c2g']**2, #5
                params['c2']*(params['kt']**2), #6
                (params['kt']**3)*params['kl'], #7
                params['kt']*params['kl']*params['c2'], #8
                params['cg']*params['kl']*params['c2'], #9
                params['c2']*params['c2g'], #10
                params['cg']*params['kl']*(params['kt']**2), #11
                params['c2g']*(params['kt']**2), #12
                (params['kl']**2)*params['cg']*params['kt'], #13
                params['c2g']*params['kt']*params['kl'], #14
                params['cg']*params['c2g']*params['kl'], #15
                (params['kt']**3)*params['cg'], #16
                params['kt']*params['c2']*params['cg'], #17
                params['kt']*(params['cg']**2)*params['kl'], #18
                params['kt']*params['cg']*params['c2g'], #19
                (params['kt']**2)*(params['cg']**2), #20
                params['c2']*(params['cg']**2), #21
                (params['cg']**3)*params['kl'], #22
                (params['cg']**2)*params['c2g'] #23
    ])
    return x

def CreateCorrectionLibfilePolyRatio(poly_ratio, errors, binning, var_names):
    corr = []
    for i, target_sample in enumerate(target_samples):
        correction_name = f"HEFT_polyratio_{target_sample}"
        corr.append(cs.Correction(
            name= correction_name,
            version=1,
            inputs=[
                cs.Variable(name="input sample", type="string", description="input sample name"),
                cs.Variable(name=var_names[0], type="real", description="HH system transverse momentum pTHH"),
                cs.Variable(name=var_names[1], type="real", description="cos(theta*)"),
                cs.Variable(name=var_names[2], type="real", description="HH invariant mass mHH"),
                cs.Variable(name="type", type="string", description="type of output, can be 'weight' or 'error'"),
            ],
            output=cs.Variable(
                name="weight", type="real", description="event-level weight"
            ),
            data=cs.Category(
                nodetype="category",
                input="input sample",
                content=[{
                    "key": key,
                    "value": cs.Binning(
                                nodetype="binning",
                                input=var_names[0],
                                edges=binning[var_names[0]],
                                content=[
                                    cs.Binning(
                                        nodetype="binning",
                                        input=var_names[1],
                                        edges=binning[var_names[1]],
                                        content=[
                                            cs.Binning(
                                                nodetype="binning",
                                                input=var_names[2],
                                                edges=binning[var_names[2]],
                                                content=[
                                                    cs.Category(
                                                        nodetype="category",
                                                        input="type",
                                                        content=[
                                                            {
                                                                "key": "weight",
                                                                "value": poly_ratio[target_sample][key][i, j, k]
                                                            },
                                                            {
                                                                "key": "error",
                                                                "value": errors[target_sample][key][i, j, k]
                                                            }
                                                        ]
                                                    )
                                                    for k in range(len(binning[var_names[2]]) - 1)],
                                                flow="clamp"
                                            ) for j in range(len(binning[var_names[1]]) - 1)
                                        ],
                                        flow="clamp"
                                    ) for i in range(len(binning[var_names[0]]) - 1)
                                ],
                                flow="clamp"
                            ), 
                } for key in poly_ratio[target_sample].keys()]
            )                     
        )
        )
    cset = cs.CorrectionSet(schema_version=2, corrections=corr)
    with open(f"{output_name}", "w") as f:
        json.dump(cset.model_dump(), f, indent=2)
    return 0


def CreateCorrectionLibfilePoly(polys, errors, binning, var_names):
    corr = []
    for i, sample in enumerate(polys.keys()):
        correction_name = f"HEFT_poly_{sample}"
        corr.append(cs.Correction(
            name= correction_name,
            version=1,
            inputs=[
                # cs.Variable(name="input sample", type="string", description="input sample name"),
                cs.Variable(name=var_names[0], type="real", description="HH system transverse momentum pTHH"),
                cs.Variable(name=var_names[1], type="real", description="cos(theta*)"),
                cs.Variable(name=var_names[2], type="real", description="HH invariant mass mHH"),
                cs.Variable(name="type", type="string", description="type of output, can be 'weight' or 'error'")
            ],
            output=cs.Variable(
                name="poly value", type="real", description="value of poly for the given sample and bin"
            ),
            data=cs.Binning(
                        nodetype="binning",
                        input=var_names[0],
                        edges=binning[var_names[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=var_names[1],
                                edges=binning[var_names[1]],
                                content=[
                                    cs.Binning(
                                        nodetype="binning",
                                        input=var_names[2],
                                        edges=binning[var_names[2]],
                                        content=[
                                            cs.Category(
                                                nodetype="category",
                                                input="type",
                                                content=[
                                                    {
                                                        "key": "weight",
                                                        "value": polys[sample][i, j, k]
                                                    },
                                                    {
                                                        "key": "error",
                                                        "value": errors[sample][i, j, k]
                                                    }
                                                ]
                                            )
                                            for k in range(len(binning[var_names[2]]) - 1)],
                                        flow="clamp"
                                    ) for j in range(len(binning[var_names[1]]) - 1)
                                ],
                                flow="clamp"
                            ) for i in range(len(binning[var_names[0]]) - 1)
                        ],
                        flow="clamp"
                    ), 
            )                     
        )

    cset = cs.CorrectionSet(schema_version=2, corrections=corr)
    with open(f"{output_name}", "w") as f:
        json.dump(cset.model_dump(), f, indent=2)
    return 0


with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

input_samples_file =  config['input_samples_file']
# input_samples = config['input_samples']
with open(input_samples_file, 'r') as f:
    samples_dict = yaml.safe_load(f)

input_file_coeff = os.path.join(".", config['2D_HEFT']['coeff'])
input_cov_file = os.path.join(".", config['2D_HEFT']['cov_matrix'])

binning, coeffs_HEFT = read_coeff_binning_from_json(input_file_coeff)
binning_, cov_HEFT = read_coeff_binning_from_json(input_cov_file)

coeffs_HEFT = order_coeffs(coeffs_HEFT, binning, key_="fitted_parameters")
cov_HEFT = order_coeffs(cov_HEFT, binning_, key_="covariance")

sigma2_coeffs = save_max_error_from_cov(cov_HEFT, binning)
print("sigma2 of coefficients from covariance matrix first bin:", sigma2_coeffs[0,0,0,:])

print(binning)
print("Coefficient shape:", coeffs_HEFT.shape)
print("coefficients, first bin (0,0,0):", coeffs_HEFT[0,0,0,:])

sm_parameters = samples_dict['signals']['GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00']
sm_parameters = prepare_params_for_poly(sm_parameters)
poly_SM = np.dot(coeffs_HEFT, sm_parameters)
error_SM = np.sqrt(np.dot(sigma2_coeffs, sm_parameters**2))

poly_signals = {}
error_signals = {}
for key in samples_dict['signals'].keys():
    signals_params = samples_dict['signals'][key]
    signals_params = prepare_params_for_poly(signals_params)
    poly_signals[key] = np.dot(coeffs_HEFT, signals_params)
    error_signals[key] = np.sqrt(np.dot(sigma2_coeffs, signals_params**2))

poly_target = {}
error_target = {}
target_samples = config['target_samples']
print("Target samples:", target_samples)
for target_sample in target_samples:
    target_parameters = samples_dict['signals'][target_sample]
    target_parameters = prepare_params_for_poly(target_parameters)
    poly_target[target_sample] = np.dot(coeffs_HEFT, target_parameters)
    error_target[target_sample] = np.sqrt(np.dot(sigma2_coeffs, target_parameters**2))

# if save_poly_ratio:
#     poly_ratio = {}
#     for target_sample in target_samples:
#         poly_ratio[target_sample] = {}
#         poly_ratio[target_sample]['GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00'] = poly_target[target_sample] / poly_SM
#         for key in poly_signals.keys():
#             if key == 'GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00' or key == target_sample or key == 'GluGlutoHHto2B2Tau_kl_0p00_kt_1p00_c2_0p00':
#                 continue
#             poly_ratio[target_sample][key] = poly_target[target_sample] / poly_signals[key]
# else:
polys = {}
errors = {}
for key in poly_signals.keys():
    if key == 'GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00' or key == 'GluGlutoHHto2B2Tau_kl_0p00_kt_1p00_c2_0p00':
        continue
    polys[key] = poly_signals[key]
    errors[key] = error_signals[key]
for key in poly_target.keys():
    polys[key] = poly_target[key]
    errors[key] = error_target[key]
    
if save_poly_ratio:
    poly_ratio = {}
    for target_sample in target_samples:
        poly_ratio[target_sample] = {}
        poly_ratio[target_sample]['GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00'] = poly_target[target_sample] / poly_SM
        errors_ratio[target_sample] = {}
        errors_ratio[target_sample]['GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00'] = np.sqrt((error_target[target_sample]/poly_SM)**2 + (poly_target[target_sample]*error_SM/poly_SM**2)**2)
        for key in poly_signals.keys():
            if key == 'GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00' or key == target_sample or key == 'GluGlutoHHto2B2Tau_kl_0p00_kt_1p00_c2_0p00':
                continue
            poly_ratio[target_sample][key] = poly_target[target_sample] / poly_signals[key]
            errors_ratio[target_sample][key] = np.sqrt((error_target[target_sample]/poly_signals[key])**2 + (poly_target[target_sample]*error_signals[key]/poly_signals[key]**2)**2)
    if not dryrun:
        CreateCorrectionLibfilePolyRatio(poly_ratio, binning, errors_ratio, var_names=['pthh','theta','mhh'])
else:
    if not dryrun:
        CreateCorrectionLibfilePoly(polys, errors, binning, var_names=['pthh','theta','mhh'])