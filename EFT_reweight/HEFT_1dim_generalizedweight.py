#This script is aimed to evaluate general part of the weight, so the part that does not depend from non-SM samples
# formula used: bin_i = poli_frac(i) / C
# C= sum_over_j poli_frac(j)*FracSM(j)

import os
import pickle
import re
import yaml
import numpy as np

def read_coefficients_1D(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()[:-1] # remove the last line which are the inclusive a coefficients
    coeffs = np.array([list(map(float, line.strip().split())) for line in lines])
    coeffs = coeffs[:, 1:] # Remove the first element of each line which is the central mhh of the bin
    coeffs = np.append(coeffs,[coeffs[-1]],axis=0) # Add the last line again to apply the last bin coefficienct to m_hh > 1400
    return coeffs

with open("./config/mHH_binedges.yaml", 'r') as f:
    binning = yaml.safe_load(f)

mhh_binedges = binning['mHH']['HEFT']['1d']
nbins_mhh = len(mhh_binedges) - 1
coeffs_HEFT = read_coefficients_1D("./Coefficients_EFT/1D_mhh_arxiv_2304_01968/HEFT_dA_and_A_with_Binning_250_1050_41_Variable_Bins_1200_1400_muR_muF_1.txt")

# kl, kt, c2, cg, c2g = HEFT_points[EFT_point].get('kl', 1.), HEFT_points[EFT_point].get('kt', 1.), HEFT_points[EFT_point].get('c2', 0.), HEFT_points[EFT_point].get('cg', 0.), HEFT_points[EFT_point].get('c2g', 0.)
kl, kt, c2, cg, c2g = 1., 2.5, 0., 0., 0.

v_parameters = np.array([kt**4,           #1
                         c2**2,           #2
                         (kt**2)*(kl**2), #3
                         (cg**2)*(kl**2), #4
                         c2g**2,          #5
                         c2*(kt**2),      #6
                         (kt**3)*kl,      #7
                         kt*kl*c2,        #8
                         cg*kl*c2,        #9
                         c2*c2g,          #10
                         cg*kl*(kt**2),   #11
                         c2g*(kt**2),     #12
                         (kl**2)*cg*kt,   #13
                         c2g*kt*kl,       #14
                         cg*c2g*kl,       #15
                         (kt**3)*cg,      #16
                         kt*c2*cg,        #17
                         kt*(cg**2)*kl,   #18
                         kt*cg*c2g,       #19
                         (kt**2)*(cg**2), #20
                         c2*(cg**2),      #21
                         (cg**3)*kl,      #22
                         (cg**2)*c2g])    #23

v_parameters_SM = np.array([1.,0.,1.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.])
poly_signal = np.dot(coeffs_HEFT,v_parameters)
poly_SM = np.dot(coeffs_HEFT,v_parameters_SM)

poly_ratio = poly_signal / poly_SM

# Need to define Cnorm
c_norm = 0.
for ibin in range(nbins_mhh):
    c_norm += poly_ratio[ibin] * FracSM_HEFT[ibin]

weight_perbin = poly_ratio

print("Weights per m_hh bin:")
for ibin in range(nbins_mhh):
    print(f"Bin {ibin}: {weight_perbin[ibin]}")
