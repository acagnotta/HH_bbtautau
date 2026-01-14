import ROOT

def read_files(dataset_path):
    os.listdir(dataset_path)
    files = []
    tchain = ROOT.TChain("Events")
    for file in os.listdir(dataset_path):
        if file.endswith(".root"):
            tchain.Add(os.path.join(dataset_path, file))
    return tchain


with open("./config/input_samples.yaml", 'r') as f:
    input_samples = yaml.safe_load(f)
