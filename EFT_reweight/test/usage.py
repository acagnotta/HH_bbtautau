import correctionlib
import time
import ROOT
import numpy as np

def read_couplings_from_string(str):
    """
    Extract coupling values from a string like:
    GluGlutoHHto*_kl_1p00_kt_0p00_c2_-0p50_*
    Returns a dict with 'kl', 'kt', 'c2' values as floats
    """
    couplings = {
        "kl": "0.0", "kt": "0.0", "c2": "0.0", "cg": "0.0", "c2g": "0.0"
    }
    
    for c in couplings.keys():
        if f"_{c}_" in str:
            c_str = str.split(f"_{c}_")[1].split("_")[0]
            couplings[c] = c_str.replace("p", ".").replace("m", "-")
    
    return couplings


if __name__ == "__main__":
    polys = "poly.json"
    cov = "cov.json"
    err = "err.json"
    
    input_sample = "GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00"
    coup_i = read_couplings_from_string(input_sample)
    target_sample = "GluGlutoHHto2B2Tau_kl_2p50_kt_1p00_c2_1p00"
    coup_o = read_couplings_from_string(target_sample)

    list_coup_i = ", ".join([f'"{coup_i[c]}"' for c in coup_i.keys()])
    list_coup_o = ", ".join([f'"{coup_o[c]}"' for c in coup_o.keys()])
    # print(list_coup_i)

    ROOT.gInterpreter.Declare('{}'.format(open("./usage.h", "r").read()))
    # Initialize RDataFrame with 10 entries and dummy variables
    np.random.seed(42)
    n_entries = 10
    
    data = {
        'LHEPart_pt': np.random.uniform(0, 100, n_entries),
        'LHEPart_eta': np.random.uniform(-2.5, 2.5, n_entries),
        'LHEPart_phi': np.random.uniform(-np.pi, np.pi, n_entries),
        'LHEPart_mass': np.random.uniform(250, 2000, n_entries),
        'LHEPart_pdgId': np.random.choice([-25, 25], n_entries),
    }
    print("Initializing RDataFrame with dummy data")
    t0 = time.time()
    df = ROOT.RDataFrame(n_entries)
    for col_name, col_values in data.items():
        df = df.Define(col_name, f"std::array<double, {n_entries}>{{{','.join(map(str, col_values))}}}[static_cast<int>(rdfentry_)]")

    df = df.Define("pthh_lhe", "GetPthhLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("mhh_lhe", "GetMhhLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("costhetastar_lhe", "GetCosthetaStarLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("weight", f'GetHEFTweight(mhh_lhe, pthh_lhe, costhetastar_lhe, {list_coup_i} , {list_coup_o}, "{polys}")')

    t1 = time.time()
    print(f"Time taken to define new columns: {t1 - t0:.2f} seconds")