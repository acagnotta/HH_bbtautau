import ROOT
import os
import yaml

def DefineKinematicGenVariables(df):
    df = df.Define("mhh_gen", "GetMhhGen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId)")
    df = df.Define("pthh_gen", "GetPthhGen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId)")
    df = df.Define("costhetastar_gen", "GetCosThetaStarGen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId)")
    return df

with open("./config/input_samples.yaml", 'r') as f:
    input_samples = yaml.safe_load(f)
with open("./config/config.yaml", 'r') as f:
    config = yaml.safe_load(f)

target_sample_name = config['target_samples']
target_sample = input_samples['signals'][target_sample_name]['path']
target_xsec = input_samples['signals'][target_sample_name]['xs']

input_files_path = os.path.join(".",input_samples['SM_signal']['GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00']['path'])
input_xsec = input_samples['SM_signal']['GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00']['xs']
print("Input files path:", input_files_path)
print("Input files pattern:", os.path.join(input_files_path, "*.root"))
if 'HLepRare' in input_files_path:
    input_files_tchain = ROOT.TChain()
    print("Adding both Events and EventsNotSelected trees for HLepRare samples")
    input_files_tchain.Add(os.path.join(input_files_path, "*.root",'Events'))
    input_files_tchain.Add(os.path.join(input_files_path, "*.root",'EventsNotSelected'))
    target_file_tchain = ROOT.TChain()
    target_file_tchain.Add(os.path.join(target_sample, "*.root",'Events'))
    target_file_tchain.Add(os.path.join(target_sample, "*.root",'EventsNotSelected'))
else:
    input_files_tchain = ROOT.TChain('Events')
    input_files_tchain.Add(os.path.join(input_files_path, "*.root"))
    target_file_tchain = ROOT.TChain('Events')
    target_file_tchain.Add(os.path.join(target_sample, "*.root"))

file = f"HEFT_reweighting_{target_sample_name}.json"
ROOT.gInterpreter.Declare('{}'.format(open("./ApplyReweighting.h", "r").read()))
df = ROOT.RDataFrame(input_files_tchain)
df = df.Define("nloweight", "LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)")
NtotalEvents = df.Sum("nloweight").GetValue()

df = DefineKinematicGenVariables(df)
df = df.Define("w_nominal", f'nloweight * {input_xsec} / {NtotalEvents}')
df = df.Define("w_reweight", f'nloweight*GetWeight(mhh_gen, pthh_gen, costhetastar_gen, "{file}") * {input_xsec} / {NtotalEvents}')
# df = df.Define("w_reweight", f'GetWeight(mhh_gen, pthh_gen, costhetastar_gen, "{file}")')


output_file = ROOT.TFile(f"reweighted_sample_{target_sample_name}.root", "RECREATE")
df.Histo1D(("mhh_gen", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen").GetValue().Write()
df.Histo1D(("mhh_sm", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen", "w_nominal").GetValue().Write()
df.Histo1D(("pthh_sm", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_gen", "w_nominal").GetValue().Write()
df.Histo1D(("costhetastar_sm", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_gen", "w_nominal").GetValue().Write()
df.Histo1D(("mhh_weighted", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen", "w_reweight").GetValue().Write()
df.Histo1D(("pthh_weighted", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_gen", "w_reweight").GetValue().Write()
df.Histo1D(("costhetastar_weighted", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_gen", "w_reweight").GetValue().Write()

output_file.Close()
print(f"Histograms saved to reweighted_sample_{target_sample_name}.root")

df = ROOT.RDataFrame(target_file_tchain)
df = df.Define("nloweight", "LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)")
NtotalEvents = df.Sum("nloweight").GetValue()
df = DefineKinematicGenVariables(df)
df = df.Define("w_nominal", f'nloweight * {target_xsec}/ {NtotalEvents}')

output_file_target = ROOT.TFile(f"original_sample_{target_sample_name}.root", "RECREATE")
df.Histo1D(("mhh_gen_target", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen").GetValue().Write()
df.Histo1D(("pthh_gen_target", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_gen").GetValue().Write()
df.Histo1D(("costhetastar_gen_target", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_gen").GetValue().Write()
df.Histo1D(("mhh_target", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen", "w_nominal").GetValue().Write()
df.Histo1D(("pthh_target", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_gen", "w_nominal").GetValue().Write()
df.Histo1D(("costhetastar_target", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_gen", "w_nominal").GetValue().Write()
output_file_target.Close()
print(f"Histograms saved to original_sample_{target_sample_name}.root")