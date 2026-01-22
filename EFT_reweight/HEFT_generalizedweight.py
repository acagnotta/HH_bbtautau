import os
import json
import re
import yaml
import numpy as np
import correctionlib.schemav2 as cs

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

def order_coeffs(coeffs_dict, binning):
    nbin_var1 = len(binning['pthh']) - 1
    nbin_var2 = len(binning['theta']) - 1
    nbin_mhh = len(binning['mhh']) - 1
    ncoeffs = 23
    coeffs_array = np.zeros((nbin_var1, nbin_var2, nbin_mhh, ncoeffs))
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
            coeffs_line = list(map(float, coeffs_dict[key][key2]["fitted_parameters"]))
            coeffs_array[idx_var1, idx_var2, idx_mhh, :] = coeffs_line
    return coeffs_array

def read_coefficients_1D(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()[:-1] # remove the last line which are the inclusive a coefficients
    coeffs = np.array([list(map(float, line.strip().split())) for line in lines])
    coeffs = coeffs[:, 1:] # Remove the first element of each line which is the central mhh of the bin
    coeffs = np.append(coeffs,[coeffs[-1]],axis=0) # Add the last line again to apply the last bin coefficienct to m_hh > 1400
    return coeffs

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

def CreateCorrectionLibfile(poly_ratio, binning, var_names):

    corr = cs.Correction(
        name="HEFT_reweighting",
        version=1,
        inputs=[
            cs.Variable(name=var_names[0], type="real", description="HH system transverse momentum pTHH"),
            cs.Variable(name=var_names[1], type="real", description="cos(theta*)"),
            cs.Variable(name=var_names[2], type="real", description="HH invariant mass mHH"),
        ],
        output=cs.Variable(
            name="weight", type="real", description="event-level weight"
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
                            content=[poly_ratio[i, j, k] for k in range(len(binning[var_names[2]]) - 1)],
                            flow="clamp"
                        ) for j in range(len(binning[var_names[1]]) - 1)
                    ],
                    flow="clamp"
                ) for i in range(len(binning[var_names[0]]) - 1)
            ],
            flow="clamp"
        ),        
    )
    cset = cs.CorrectionSet(schema_version=2, corrections=[corr])
    with open(f"HEFT_reweighting_{target_sample}.json", "w") as f:
        json.dump(cset.model_dump(), f, indent=2)
    return 0

with open("./config/mHH_binedges.yaml", 'r') as f:
    config = yaml.safe_load(f)

input_samples =  os.path.join(".",config['input_samples'])
with open(input_samples, 'r') as f:
    samples_dict = yaml.safe_load(f)

input_file_coeff = os.path.join(".", config['2D_HEFT']['coeff'])

binning, coeffs_HEFT = read_coeff_binning_from_json(input_file_coeff)

coeffs_HEFT = order_coeffs(coeffs_HEFT, binning)

print(binning)
print("Coefficient shape:", coeffs_HEFT.shape)
print("coefficients, first bin (0,0,0):", coeffs_HEFT[0,0,0,:])
# print("coefficients, bin (0,0,1):", coeffs_HEFT[0,0,1,:])
# print("coefficients, bin (0,1,0):", coeffs_HEFT[0,1,0,:])
# print("coefficients, bin (1,0,0):", coeffs_HEFT[1,0,0,:])


sm_parameters = samples_dict['SM_signal']['GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00']
sm_parameters = prepare_params_for_poly(sm_parameters)

# kl, kt, c2, cg, c2g = HEFT_points[EFT_point].get('kl', 1.), HEFT_points[EFT_point].get('kt', 1.), HEFT_points[EFT_point].get('c2', 0.), HEFT_points[EFT_point].get('cg', 0.), HEFT_points[EFT_point].get('c2g', 0.)

target_sample = config['target_samples']
target_parameters = samples_dict['signals'][target_sample]
target_parameters = prepare_params_for_poly(target_parameters)

poly_target = np.dot(coeffs_HEFT, target_parameters)
poly_SM = np.dot(coeffs_HEFT, sm_parameters)
poly_ratio = poly_target / poly_SM

# print("Polynomial ratio shape:", poly_ratio.shape)
# print("Polynomial ratio, first bin (0,0,0):", poly_ratio[0,0,0])
# print("Polynomial ratio, first bin (1,0,0):", poly_ratio[1,0,0])
# print("Polynomial ratio, first bin (0,1,0):", poly_ratio[0,1,0])
# print("Polynomial ratio, first bin (0,0,1):", poly_ratio[0,0,1])


CreateCorrectionLibfile(poly_ratio, binning, var_names=['pthh','theta','mhh'])

# with open(f"HEFT_reweighting_{target_sample}.json", "w") as f:
#     json.dump(corr_json.model_dump(),f, indent=4)