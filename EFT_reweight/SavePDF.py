import ROOT, os
import yaml, json
import array
import numpy as np
import correctionlib
import correctionlib.schemav2 as cs

def DefineKinematicGenVariables(df):
    df = df.Define("mhh", "GetMhhLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("pthh", "GetPthhLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId)")
    df = df.Define("costhetastar", "abs(GetCosThetaStarLHE(LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId))")
    return df

config_path = "./config/config.yaml"
with open(config_path, 'r') as f:
    config = yaml.safe_load(f)
with open(config['input_samples_file'], 'r') as f:
    input_samples = yaml.safe_load(f)

input_files = config['target_samples'] 

RSamples_inputs = ROOT.RDF.Experimental.RDatasetSpec()
input_files_path = {}
for s in input_files:
    input_files_path[s] = os.path.join(".",input_samples['signals'][s]['path'])
    tmp_sample = ROOT.RDF.Experimental.RSample(s, "Events", os.path.join(input_files_path[s], "*.root"))
    RSamples_inputs = RSamples_inputs.AddSample(tmp_sample)
    tmp_sample = ROOT.RDF.Experimental.RSample(s, "EventsNotSelected", os.path.join(input_files_path[s], "*.root"))
    RSamples_inputs = RSamples_inputs.AddSample(tmp_sample)

ntotEvents_map = "std::map<std::string, double> NtotalEvents = {};"
xsec_map_code = "std::map<std::string, double> xSec = {};"
ROOT.gInterpreter.Declare(xsec_map_code)
ROOT.gInterpreter.Declare(ntotEvents_map)
ROOT.gInterpreter.Declare('{}'.format(open("./ApplyReweighting.h", "r").read()))
df = ROOT.RDataFrame(RSamples_inputs)
df = DefineKinematicGenVariables(df)
df = df.Define("nloweight", "genWeight/abs(genWeight)")
NtotalEventsW = df.Sum("nloweight").GetValue()
df = df.Define("NtotalEventsW", f"{NtotalEventsW}")
print(f"Total number of events (weighted) in input samples: {NtotalEventsW}")
NtotalEvents = df.Count().GetValue()
df = df.Define("NtotalEvents", f"{NtotalEvents}")
print(f"Total number of events (unweighted) in input samples: {NtotalEvents}")

binning = {
    "mhh": array.array('d', [250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 650, 670, 690, 710, 730, 750, 800, 850, 900, 950, 1000, 1050, 1150, 1200, 1300, 1400, 1500, 1600, 1800, 2000]),
    "pthh": array.array('d', [0,20,40,60,80,100,140,200,290,2500]),
    "costhetastar": array.array('d', [0.0, 0.25, 0.5, 0.75, 1.0])
}
df = df.Define("w", 'nloweight/ NtotalEventsW')
df = df.Define("w_", '1./ NtotalEvents')

h1 = df.Histo3D(("h3Dweighted", ";m_{HH} [GeV]; p_{T}^{HH} [GeV]; cos(#theta^{*})", len(binning["mhh"])-1, binning["mhh"], len(binning["pthh"])-1, binning["pthh"], len(binning["costhetastar"])-1, binning["costhetastar"]), "mhh", "pthh", "costhetastar", "w").GetValue()
h2 = df.Histo3D(("h3Dunweighted", ";m_{HH} [GeV]; p_{T}^{HH} [GeV]; cos(#theta^{*})", len(binning["mhh"])-1, binning["mhh"], len(binning["pthh"])-1, binning["pthh"], len(binning["costhetastar"])-1, binning["costhetastar"]), "mhh", "pthh", "costhetastar", "w_").GetValue()




count_flat = np.array([h1.GetBinContent(i+1, j+1, k+1) for i in range(len(binning["mhh"])-1) for j in range(len(binning["pthh"])-1) for k in range(len(binning["costhetastar"])-1)])
count_unweighted_flat = np.array([h2.GetBinContent(i+1, j+1, k+1) for i in range(len(binning["mhh"])-1) for j in range(len(binning["pthh"])-1) for k in range(len(binning["costhetastar"])-1)])

count = count_flat.reshape((len(binning["mhh"])-1, len(binning["pthh"])-1, len(binning["costhetastar"])-1))
count_unweighted = count_unweighted_flat.reshape((len(binning["mhh"])-1, len(binning["pthh"])-1, len(binning["costhetastar"])-1))

print(count.shape, count_unweighted.shape)
total_count = np.sum(count)
total_count_unweighted = np.sum(count_unweighted)
print(f"Total count (weighted) from histogram: {total_count}")
print(f"Total count (unweighted) from histogram: {total_count_unweighted}")

corr = []
corr.append(cs.Correction(
        name= "PDF_weighted_target",
        version=1,
        inputs=[
            # cs.Variable(name="input sample", type="string", description="input sample name"),
            cs.Variable(name="pthh", type="real", description="HH system transverse momentum pTHH"),
            cs.Variable(name="costhetastar", type="real", description="abs(cos(theta*))"),
            cs.Variable(name="mhh", type="real", description="HH invariant mass mHH"),
        ],
        output=cs.Variable(
            name="PDF value", type="real", description="value of PDF for the given sample and bin"
        ),
        data=cs.Binning(
                    nodetype="binning",
                    input="mhh",
                    edges=binning["mhh"],
                    content=[
                            cs.Binning(
                                nodetype="binning",
                                input="pthh",
                                edges=binning["pthh"],
                                content=[
                                    cs.Binning(
                                        nodetype="binning",
                                        input="costhetastar",
                                        edges=binning["costhetastar"],
                                        content=[count[i, j, k] for k in range(len(binning["costhetastar"]) - 1)],
                                        flow="clamp"
                                    ) for j in range(len(binning["pthh"]) - 1)
                                ],
                                flow="clamp"
                            ) for i in range(len(binning["mhh"]) - 1)
                        ],
                        flow="clamp"
                    ), 
            )                     
        )
corr.append(cs.Correction(
        name= "PDF_target",
        version=1,
        inputs=[
            # cs.Variable(name="input sample", type="string", description="input sample name"),
            cs.Variable(name="pthh", type="real", description="HH system transverse momentum pTHH"),
            cs.Variable(name="costhetastar", type="real", description="abs(cos(theta*))"),
            cs.Variable(name="mhh", type="real", description="HH invariant mass mHH"),
        ],
        output=cs.Variable(
            name="PDF value", type="real", description="value of PDF for the given sample and bin"
        ),
        data=cs.Binning(
                    nodetype="binning",
                    input="mhh",
                    edges=binning["mhh"],
                    content=[
                            cs.Binning(
                                nodetype="binning",
                                input="pthh",
                                edges=binning["pthh"],
                                content=[
                                    cs.Binning(
                                        nodetype="binning",
                                        input="costhetastar",
                                        edges=binning["costhetastar"],
                                        content=[count_unweighted[i, j, k] for k in range(len(binning["costhetastar"]) - 1)],
                                        flow="clamp"
                                    ) for j in range(len(binning["pthh"]) - 1)
                                ],
                                flow="clamp"
                            ) for i in range(len(binning["mhh"]) - 1)
                        ],
                        flow="clamp"
                    ), 
            )                     
        )

cset = cs.CorrectionSet(schema_version=2, corrections=corr)
with open(f"pdf_values_target.json", "w") as f:
    json.dump(cset.model_dump(), f, indent=2)

# for i, mhh_edge in enumerate(binning["mhh"][:-1]):
#     for j, pthh_edge in enumerate(binning["pthh"][:-1]):
#         for k, costhetastar_edge in enumerate(binning["costhetastar"][:-1]):
#             mhh_min, mhh_max = binning["mhh"][i], binning["mhh"][i+1]
#             pthh_min, pthh_max = binning["pthh"][j], binning["pthh"][j+1]
#             costhetastar_min, costhetastar_max = binning["costhetastar"][k], binning["costhetastar"][k+1]
#             idx = i * (len(binning["pthh"])-1) * (len(binning["costhetastar"])-1) + j * (len(binning["costhetastar"])-1) + k
#             f.write(f"Bin ({i+1},{j+1},{k+1}): mhh=[{mhh_min},{mhh_max}], pthh=[{pthh_min},{pthh_max}], costhetastar=[{costhetastar_min},{costhetastar_max}] -> PDF value (weighted) = {count[idx]:.6e}, PDF value (unweighted) = {count_unweighted[idx]:.6e}\n")
