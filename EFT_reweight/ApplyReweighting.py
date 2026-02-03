import ROOT
import os
import yaml
ROOT.EnableImplicitMT()

doTarget = False

def DefineKinematicGenVariables(df):
    df = df.Define("mhh_gen", "GetMhhGen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId)")
    df = df.Define("pthh_gen", "GetPthhGen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId)")
    df = df.Define("costhetastar_gen", "GetCosThetaStarGen(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_statusFlags, GenPart_pdgId)")

    df = df.Define("mhh_LHE", "GetMhhLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("pthh_LHE", "GetPthhLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("costhetastar_LHE", "GetCosThetaStarLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    return df

with open("./config/input_samples.yaml", 'r') as f:
    input_samples = yaml.safe_load(f)
with open("./config/config.yaml", 'r') as f:
    config = yaml.safe_load(f)

target_sample_name = config['target_samples']
target_sample = input_samples['signals'][target_sample_name]['path']
target_xsec = input_samples['signals'][target_sample_name]['xs']

input_files_path = {}
input_xsec = {}
RSamples_inputs = ROOT.RDF.Experimental.RDatasetSpec()
for s in config['input_samples']:
    input_files_path[s] = os.path.join(".",input_samples['signals'][s]['path'])
    print(f"Input sample {s} path: {input_files_path[s]}")
    input_xsec[s] = input_samples['signals'][s]['xs']
    if "HLepRare" in input_files_path[s]:
        tmp_sample = ROOT.RDF.Experimental.RSample(s, "Events", os.path.join(input_files_path[s], "*.root"))
        RSamples_inputs = RSamples_inputs.AddSample(tmp_sample)
        tmp_sample = ROOT.RDF.Experimental.RSample(s, "EventsNotSelected", os.path.join(input_files_path[s], "*.root"))
        RSamples_inputs = RSamples_inputs.AddSample(tmp_sample)
    else:
        tmp_sample = ROOT.RDF.Experimental.RSample(s, "Events", os.path.join(input_files_path[s], "*.root"))
        RSamples_inputs = RSamples_inputs.AddSample(tmp_sample)

RSamples_target = ROOT.RDF.Experimental.RDatasetSpec()
if "HLepRare" in target_sample:
    tmp_sample = ROOT.RDF.Experimental.RSample(target_sample_name, "Events", os.path.join(target_sample, "*.root"))
    RSamples_target = RSamples_target.AddSample(tmp_sample)
    tmp_sample = ROOT.RDF.Experimental.RSample(target_sample_name, "EventsNotSelected", os.path.join(target_sample, "*.root"))
    RSamples_target = RSamples_target.AddSample(tmp_sample)
else:
    tmp_sample = ROOT.RDF.Experimental.RSample(target_sample_name, "Events", os.path.join(target_sample, "*.root"))
    RSamples_target = RSamples_target.AddSample(tmp_sample)

# if 'HLepRare' in input_files_path:
#     input_files_tchain = ROOT.TChain()
#     print("Adding both Events and EventsNotSelected trees for HLepRare samples")
#     input_files_tchain.Add(os.path.join(input_files_path, "*.root",'Events'))
#     input_files_tchain.Add(os.path.join(input_files_path, "*.root",'EventsNotSelected'))
#     target_file_tchain = ROOT.TChain()
#     target_file_tchain.Add(os.path.join(target_sample, "*.root",'Events'))
#     target_file_tchain.Add(os.path.join(target_sample, "*.root",'EventsNotSelected'))
# else:
#     input_files_tchain = ROOT.TChain('Events')
#     input_files_tchain.Add(os.path.join(input_files_path, "*.root"))
#     target_file_tchain = ROOT.TChain('Events')
#     target_file_tchain.Add(os.path.join(target_sample, "*.root"))

xsec_map_code = "std::map<std::string, double> xSec = {"
xsec_entries = []
for sample, xsec in input_xsec.items():
    xsec_entries.append(f'{{"{sample}", {xsec}}}')
xsec_map_code += ", ".join(xsec_entries) + "};"

ROOT.gInterpreter.Declare(xsec_map_code)


file = f"HEFT_reweighting.json"
ROOT.gInterpreter.Declare('{}'.format(open("./ApplyReweighting.h", "r").read()))

print("Creating RDataFrame for input samples...")
df = ROOT.RDataFrame(RSamples_inputs)
df = df.DefinePerSample("SampleName", "rdfsampleinfo_.GetSampleName()")

df = df.DefinePerSample("input_xsec", "GetInputXSec(rdfsampleinfo_.GetSampleName())")
df = df.Define("nloweight", "LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)")
NtotalEvents = df.Sum("nloweight").GetValue()
print(f"Total number of events (weighted) in input samples: {NtotalEvents}")
df = DefineKinematicGenVariables(df)
df = df.Define("w_nominal", f'nloweight * input_xsec / {NtotalEvents}')
df = df.Define("w_reweight", f'nloweight * GetWeight(mhh_gen, pthh_gen, costhetastar_gen, SampleName, "{file}", "{target_sample_name}") * input_xsec / {NtotalEvents}')
df = df.Define("w_reweightLHE", f'nloweight * GetWeight(mhh_LHE, pthh_LHE, costhetastar_LHE, SampleName, "{file}", "{target_sample_name}") * input_xsec / {NtotalEvents}')
# df = df.Define("w_reweight", f'GetWeight(mhh_gen, pthh_gen, costhetastar_gen, "{file}")')


output_file = ROOT.TFile(f"reweighted_sample_{target_sample_name}.root", "RECREATE")
df.Histo1D(("mhh_sm", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen", "w_nominal").GetValue().Write()
df.Histo1D(("pthh_sm", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_gen", "w_nominal").GetValue().Write()
df.Histo1D(("costhetastar_sm", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_gen", "w_nominal").GetValue().Write()
df.Histo1D(("mhh_weighted", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen", "w_reweight").GetValue().Write()
df.Histo1D(("pthh_weighted", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_gen", "w_reweight").GetValue().Write()
df.Histo1D(("costhetastar_weighted", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_gen", "w_reweight").GetValue().Write()
df.Histo1D(("mhh_weightedLHE", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_LHE", "w_reweightLHE").GetValue().Write()
df.Histo1D(("pthh_weightedLHE", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_LHE", "w_reweightLHE").GetValue().Write()
df.Histo1D(("costhetastar_weightedLHE", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_LHE", "w_reweightLHE").GetValue().Write()

output_file.Close()
print(f"Histograms saved to reweighted_sample_{target_sample_name}.root")

if doTarget:
    df_target = ROOT.RDataFrame(RSamples_target)
    df_target = df_target.Define("nloweight", "genWeight/abs(genWeight)")
    NtotalEvents_target = df_target.Sum("nloweight").GetValue()
    df_target = DefineKinematicGenVariables(df_target)
    df_target = df_target.Define("w_nominal", f'nloweight * {target_xsec}/ {NtotalEvents_target}')

    output_file_target = ROOT.TFile(f"original_sample_{target_sample_name}.root", "RECREATE")
    df_target.Histo1D(("mhh_gen_target", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen").GetValue().Write()
    df_target.Histo1D(("pthh_gen_target", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_gen").GetValue().Write()
    df_target.Histo1D(("costhetastar_gen_target", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_gen").GetValue().Write()
    df_target.Histo1D(("mhh_target", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_gen", "w_nominal").GetValue().Write()
    df_target.Histo1D(("pthh_target", ";p_{T,HH} [GeV]; Events", 200, 0, 2000), "pthh_gen", "w_nominal").GetValue().Write()
    df_target.Histo1D(("costhetastar_target", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar_gen", "w_nominal").GetValue().Write()
    df_target.Histo1D(("mhh_LHE_target", ";m_{HH} [GeV]; Events", 75, 250, 1000), "mhh_LHE", "w_nominal").GetValue().Write()
    output_file_target.Close()
    print(f"Histograms saved to original_sample_{target_sample_name}.root")