import ROOT
import numpy as np

input_files = [
    "plots/anaTuples_werrors_kl_0p00_kt_1p00_c2_1p00input.root",
    "plots/anaTuples_werrors_kl_1p00_kt_1p00_c2_0p10input.root",
    "plots/anaTuples_werrors_kl_1p00_kt_1p00_c2_3p00input.root",
    "plots/anaTuples_werrors_kl_1p00_kt_1p00_c2_m2p00input.root",
    "plots/anaTuples_werrors_kl2p45_kt_1p00_c2_0p00input.root",
    "plots/anaTuples_werrors_kl_5p00_kt_1p00_c2_0p00input.root",
    "plots/anaTuples_werrors_sminput.root"
]

def weighted_avg(x, w):
    return np.sum(x * w) / np.sum(w)

def weighted_histo(histos, var): #stat only weighted average
    h_out = histos[0].Clone(var)
    h_out.Reset()
    print(h_out.GetBinContent(1), h_out.GetBinError(1))
    print("Calculating weighted average...")
    for iBin in range(1, h_out.GetNbinsX() + 1):
        x = np.array([h.GetBinContent(iBin) for h in histos])
        w = np.array([1/(h.GetBinError(iBin)**2) if h.GetBinError(iBin) != 0 else 0 for h in histos])
        avg = weighted_avg(x, w)
        h_out.SetBinContent(iBin, avg)
        h_out.SetBinError(iBin, np.sqrt(1/np.sum(w)))
    return h_out

def weighted_histo_withsyst(histos, histosUP, histosDOWN, var): #stat + variation weighted average
    h_out = histos[0].Clone(var+"_wsyst")
    h_out.Reset()
    print(h_out.GetBinContent(1), h_out.GetBinError(1))
    print("Calculating weighted average with systematics...")
    for iBin in range(1, h_out.GetNbinsX() + 1):
        x = np.array([h.GetBinContent(iBin) for h in histos])
        w = np.array([1/(h.GetBinError(iBin)**2 + ((h.GetBinContent(iBin) - hdown.GetBinContent(iBin))**2 + (hup.GetBinContent(iBin) - h.GetBinContent(iBin))**2)/2) if h.GetBinError(iBin) != 0 else 0 for (h, hup, hdown) in zip(histos, histosUP, histosDOWN)])
        avg = weighted_avg(x, w)
        h_out.SetBinContent(iBin, avg)
        h_out.SetBinError(iBin, np.sqrt(1/np.sum(w)))
    return h_out
    

tfiles = [ROOT.TFile.Open(f) for f in input_files]
histos_mhh = [f.Get("mhh_weighted") for f in tfiles]
histosUP_mhh = [f.Get("mhh_weighted_up") for f in tfiles]
histosDOWN_mhh = [f.Get("mhh_weighted_down") for f in tfiles]
histos_pthh = [f.Get("pthh_weighted") for f in tfiles]
histosUP_pthh = [f.Get("pthh_weighted_up") for f in tfiles]
histosDOWN_pthh = [f.Get("pthh_weighted_down") for f in tfiles]
histos_theta = [f.Get("costhetastar_weighted") for f in tfiles]
histosUP_theta = [f.Get("costhetastar_weighted_up") for f in tfiles]
histosDOWN_theta = [f.Get("costhetastar_weighted_down") for f in tfiles]

h_mhh = weighted_histo(histos_mhh, "mhh_weighted")
# h_mhh_ = weighted_histo_withsyst(histos_mhh, histosUP_mhh, histosDOWN_mhh, "mhh_weighted")
h_mhh_up = weighted_histo(histosUP_mhh, "mhh_weighted_up")
h_mhh_down = weighted_histo(histosDOWN_mhh, "mhh_weighted_down")

h_pthh = weighted_histo(histos_pthh, "pthh_weighted")
# h_pthh_ = weighted_histo_withsyst(histos_pthh, histosUP_pthh, histosDOWN_pthh, "pthh_weighted")
h_pthh_up = weighted_histo(histosUP_pthh, "pthh_weighted_up")
h_pthh_down = weighted_histo(histosDOWN_pthh, "pthh_weighted_down")

h_theta = weighted_histo(histos_theta, "costhetastar_weighted")
# h_theta_ = weighted_histo_withsyst(histos_theta, histosUP_theta, histosDOWN_theta, "costhetastar_weighted")
h_theta_up = weighted_histo(histosUP_theta, "costhetastar_weighted_up")
h_theta_down = weighted_histo(histosDOWN_theta, "costhetastar_weighted_down")



output_file = ROOT.TFile.Open("plots/anaTuples_werrors_weighted_avg.root", "RECREATE")
h_mhh.Write()
# h_mhh_.Write()
h_mhh_up.Write()
h_mhh_down.Write()
h_pthh.Write()
# h_pthh_.Write()
h_pthh_up.Write()
h_pthh_down.Write()
h_theta.Write()
# h_theta_.Write()
h_theta_up.Write()
h_theta_down.Write()
output_file.Close()