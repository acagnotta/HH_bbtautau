import ROOT
import correctionlib
import optparse
import array

usage                   = 'python3 PDFreweight.py '
parser                  = optparse.OptionParser(usage)
parser.add_option('-i', '--InputjsonFile', dest='inputjsonFile', type=str, default="./data/pdf_values.json", help='Path to the JSON file containing the pdf information input sample (output of SavePDF.py)')
parser.add_option('-t', '--TargetjsonFile', dest='targetjsonFile', type=str, default="./data/pdf_values.json", help='Path to the JSON file containing the pdf information target sample (output of SavePDF.py)')
parser.add_option('--weightsJsonFile', dest='weightsJsonFile', type=str, default="./data/HEFT_poly.json", help='Path to the JSON file containing the pdf information target sample (output of SavePDF.py)')
parser.add_option('-o', '--output', dest='output', type=str, default="./plots/pdfhisto.root", help='Path to the output ROOT file')
(opt, args) = parser.parse_args()


binning = {
    "mhh": [250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 650, 670, 690, 710, 730, 750, 800, 850, 900, 950, 1000, 1050, 1150, 1200, 1300, 1400, 1500, 1600, 1800, 2000],
    "pthh": [0, 20, 40, 60, 80, 100, 140, 200, 290, 2500],
    "costhetastar": [0.0, 0.25, 0.5, 0.75, 1.0]
}

corr_input = correctionlib.CorrectionSet.from_file(opt.inputjsonFile)
corr_target = correctionlib.CorrectionSet.from_file(opt.targetjsonFile)
weights_corr = correctionlib.CorrectionSet.from_file(opt.weightsJsonFile)
sigma_SM = 0.034170
sigma_target = 0.010448

output_file = ROOT.TFile(opt.output, "RECREATE")
h_mhh_output = ROOT.TH1D("h_mhh_output", "m_{HH} distribution after PDF reweighting; m_{HH} [GeV]; Events", len(binning["mhh"])-1, array.array('d', binning["mhh"]))
h_mhh_input = ROOT.TH1D("h_mhh_input", "m_{HH} distribution before PDF reweighting; m_{HH} [GeV]; Events", len(binning["mhh"])-1, array.array('d', binning["mhh"]))
h_mhh_target = ROOT.TH1D("h_mhh_target", "m_{HH} distribution of target sample; m_{HH} [GeV]; Events", len(binning["mhh"])-1, array.array('d', binning["mhh"]))
for i in range(len(binning["mhh"])-1):
    for j in range(len(binning["pthh"])-1):
        for k in range(len(binning["costhetastar"])-1):
            mhh_center = (binning["mhh"][i] + binning["mhh"][i+1]) / 2
            pthh_center = (binning["pthh"][j] + binning["pthh"][j+1]) / 2
            costhetastar_center = (binning["costhetastar"][k] + binning["costhetastar"][k+1]) / 2
            input_pdf_value = corr_input["PDF_weighted_SM"].evaluate(pthh_center, costhetastar_center, mhh_center)
            target_pdf_value = corr_target["PDF_weighted_target"].evaluate(pthh_center, costhetastar_center, mhh_center)
            poly_sm = weights_corr["HEFT_poly_GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p00"].evaluate(pthh_center, costhetastar_center, mhh_center)
            poly_i = weights_corr["HEFT_poly_GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p35"].evaluate(pthh_center, costhetastar_center, mhh_center)
            poly_ratio = poly_i / poly_sm if poly_sm != 0 else 1.0
            h_mhh_target.Fill(mhh_center, target_pdf_value * sigma_target)
            h_mhh_input.Fill(mhh_center, input_pdf_value)
            h_mhh_output.Fill(mhh_center, input_pdf_value * poly_ratio * sigma_SM)
output_file.cd()
h_mhh_input.Write()
h_mhh_target.Write()
h_mhh_output.Write()
output_file.Close()
            