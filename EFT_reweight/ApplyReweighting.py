import ROOT
import os
import yaml
import array as array
import optparse
ROOT.EnableImplicitMT(4)

usage                   = 'python3 ApplyReweighting.py '
parser                  = optparse.OptionParser(usage)
parser.add_option('--doTarget', dest='doTarget', action='store_true', default=False, help='If true, the script will process the target sample. Enable this option only if doReweight is false')
parser.add_option('--doReweight', dest='doReweight', action='store_true', default=False, help='If true, the script will process the input samples and apply reweighting. Enable this option only if doTarget is false')
parser.add_option('--jsonFile', dest='jsonFile', type=str, default="./data/HEFT_reweighting.json", help='Path to the JSON file containing the reweighting information (output of HEFT_generalizedweight.py)')
parser.add_option('-c', '--config', dest='config', type=str, default="./config/config.yaml", help='Path to the config file')
parser.add_option('-o', '--output', dest='output', type=str, default="./plots/reweighted_sample.root", help='Path to the output ROOT file (of reweighted sample or target sample if --doTarget is set)')
(opt, args) = parser.parse_args()
doTarget = opt.doTarget
doReweight = opt.doReweight
output_name = opt.output
config_path = opt.config
json_file = opt.jsonFile
if doTarget and doReweight:
    print("Error: Both doTarget and doReweight options cannot be true at the same time. Please choose one of them.")
    exit(1)

def DefineKinematicGenVariables(df):
    df = df.Define("mhh", "GetMhhLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("pthh", "GetPthhLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("costhetastar", "GetCosThetaStarLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    return df

with open(config_path, 'r') as f:
    config = yaml.safe_load(f)
with open(config['input_samples_file'], 'r') as f:
    input_samples = yaml.safe_load(f)

target_sample_name = config['target_samples'][1]
target_sample = input_samples['signals'][target_sample_name]['path']
target_xsec = input_samples['signals'][target_sample_name]['xs']

input_anatuples_file = False
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
        if 'anaTuples' in input_files_path[s]:
            input_anatuples_file = True

RSamples_target = ROOT.RDF.Experimental.RDatasetSpec()
if "HLepRare" in target_sample:
    tmp_sample = ROOT.RDF.Experimental.RSample(target_sample_name, "Events", os.path.join(target_sample, "*.root"))
    RSamples_target = RSamples_target.AddSample(tmp_sample)
    tmp_sample = ROOT.RDF.Experimental.RSample(target_sample_name, "EventsNotSelected", os.path.join(target_sample, "*.root"))
    RSamples_target = RSamples_target.AddSample(tmp_sample)
else:
    tmp_sample = ROOT.RDF.Experimental.RSample(target_sample_name, "Events", os.path.join(target_sample, "*.root"))
    RSamples_target = RSamples_target.AddSample(tmp_sample)
    if 'anaTuples' in target_sample:
        input_anatuples_file = True

mhh_binning = array.array('d', [250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 650, 670, 690, 710, 730, 750, 800, 850, 900, 950, 1000, 1050, 1150, 1200, 1300, 1400, 1500, 1600, 1800, 2000])
pthh_binning = array.array('d', [0, 20, 40, 60, 80, 100, 140, 200, 290, 2500])
if doReweight:
    file = json_file

    print("Creating RDataFrame for input samples...")
    df = ROOT.RDataFrame(RSamples_inputs)
    if not input_anatuples_file:
        df = df.Define("nloweight", "genWeight/abs(genWeight)")
    df = df.DefinePerSample("SampleName", "rdfsampleinfo_.GetSampleName()")
    
    NtotalEvents = {}
    for s in config['input_samples']:
        df_sample = df.Filter(f'SampleName == "{s}"')
        if input_anatuples_file:
            Nsample = -1
        else:
            Nsample = df_sample.Sum("nloweight").GetValue()
        print(f"Total number of events (weighted) in input sample {s}: {Nsample}")
        NtotalEvents[s] = Nsample
    ntotEvents_map = "std::map<std::string, double> NtotalEvents = {"
    ntot_entries = []
    for sample, nEvents in NtotalEvents.items():
        ntot_entries.append(f'{{"{sample}", {nEvents}}}')
    ntotEvents_map += ", ".join(ntot_entries) + "};"

    xsec_map_code = "std::map<std::string, double> xSec = {"
    xsec_entries = []
    for sample, xsec in input_xsec.items():
        xsec_entries.append(f'{{"{sample}", {xsec}}}')
    xsec_map_code += ", ".join(xsec_entries) + "};"

    ROOT.gInterpreter.Declare(xsec_map_code)
    ROOT.gInterpreter.Declare(ntotEvents_map)
    ROOT.gInterpreter.Declare('{}'.format(open("./ApplyReweighting.h", "r").read()))

    df = df.DefinePerSample("input_xsec", "GetInputXSec(rdfsampleinfo_.GetSampleName())")
    df = DefineKinematicGenVariables(df)
    if input_anatuples_file:
        weight = "weight_MC_Lumi_pu * input_xsec"
    elif len(config['input_samples']) == 1:
        weight = f"nloweight * input_xsec / {NtotalEvents}"
        NtotalEvents = df.Sum("nloweight").GetValue()
        print(f"Total number of events (weighted) in input samples: {NtotalEvents}")
    else:
        weight = f'(1./{str(len(config["input_samples"]))}) * nloweight * input_xsec / NtotalEvents'
    
    if len(config['input_samples'])==1:
        df = df.Define("w_nominal", weight)
        df = df.Define("w_reweight", f'{weight} * GetWeightFromPoly(mhh, pthh, costhetastar, SampleName, "{file}", "{target_sample_name}")')
    else:
        df = df.DefinePerSample(f"NtotalEvents", f'GetNtotalEvents(rdfsampleinfo_.GetSampleName())')
        df = df.Define("w_nominal", f'nloweight * input_xsec / NtotalEvents')
        df = df.Define("w_reweight", f'(1./{str(len(config["input_samples"]))}) * nloweight * GetWeightFromPoly(mhh, pthh, costhetastar, SampleName, "{file}", "{target_sample_name}") * input_xsec / NtotalEvents')

    output_file = ROOT.TFile(f"{output_name}", "RECREATE")
    df.Histo1D(("mhh_nominal", ";m_{HH} [GeV]; Events", len(mhh_binning)-1, mhh_binning), "mhh", "w_nominal").GetValue().Write()
    df.Histo1D(("mhh_weighted", ";m_{HH} [GeV]; Events", len(mhh_binning)-1, mhh_binning), "mhh", "w_reweight").GetValue().Write()
    df.Histo1D(("pthh_weighted", ";p_{T,HH} [GeV]; Events", len(pthh_binning)-1, pthh_binning), "pthh", "w_reweight").GetValue().Write()
    df.Histo1D(("costhetastar_weighted", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar", "w_reweight").GetValue().Write()

    output_file.Close()
    print(f"Histograms saved to {output_name}")

if doTarget:
    df_target = ROOT.RDataFrame(RSamples_target)
    if not input_anatuples_file:
        df_target = df_target.Define("nloweight", "genWeight/abs(genWeight)")
    df_target = df_target.DefinePerSample("SampleName", "rdfsampleinfo_.GetSampleName()")

    NtotalEvents = {}
    for s in config['input_samples']:
        df_sample = df_target.Filter(f'SampleName == "{s}"')
        if input_anatuples_file:
            Nsample = -1
        else:
            Nsample = df_sample.Sum("nloweight").GetValue()
        print(f"Total number of events (weighted) in input sample {s}: {Nsample}")
        NtotalEvents[s] = Nsample
    ntotEvents_map = "std::map<std::string, double> NtotalEvents = {"
    ntot_entries = []
    for sample, nEvents in NtotalEvents.items():
        ntot_entries.append(f'{{"{sample}", {nEvents}}}')
    ntotEvents_map += ", ".join(ntot_entries) + "};"

    xsec_map_code = "std::map<std::string, double> xSec = {"
    xsec_entries = []
    for sample, xsec in input_xsec.items():
        xsec_entries.append(f'{{"{sample}", {xsec}}}')
    xsec_map_code += ", ".join(xsec_entries) + "};"

    ROOT.gInterpreter.Declare(xsec_map_code)
    ROOT.gInterpreter.Declare(ntotEvents_map)
    ROOT.gInterpreter.Declare('{}'.format(open("./ApplyReweighting.h", "r").read()))
    df_target = DefineKinematicGenVariables(df_target)
    if input_anatuples_file:
        weight = f"weight_MC_Lumi_pu * {target_xsec}"
    else:
        NtotalEvents_target = df_target.Sum("nloweight").GetValue()
        weight = f"nloweight * {target_xsec} / {NtotalEvents_target}"
    df_target = df_target.Define("w_nominal", f'{weight}')

    output_file_target = ROOT.TFile(f"{output_name}", "RECREATE")
    df_target.Histo1D(("mhh_target", ";m_{HH} [GeV]; Events", len(mhh_binning)-1, mhh_binning), "mhh", "w_nominal").GetValue().Write()
    df_target.Histo1D(("pthh_target", ";p_{T,HH} [GeV]; Events", len(pthh_binning)-1, pthh_binning), "pthh", "w_nominal").GetValue().Write()
    df_target.Histo1D(("costhetastar_target", ";cos(#theta^{*}); Events", 20, -1, 1), "costhetastar", "w_nominal").GetValue().Write()
    output_file_target.Close()
    print(f"Histograms saved to {output_name}")