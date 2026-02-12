import cmsstyle as CMS
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

def plot(h, folder, fillcolor, canv_name = "canv" ,extraTest="Preliminary", iPos=0, energy=13.6, lumi = 1,  addInfo="", ytitle = "Events"):
    if type(h)==list:
        h1 = h[0]
        # hist_dict = [k.GetName() for k in h]
    else:
        h1 = h
    CMS.SetExtraText(extraTest)
    iPos = iPos
    canv_name = canv_name
    CMS.SetLumi(lumi,"pb","",1)
    CMS.SetEnergy(energy)
    CMS.ResetAdditionalInfo()
    CMS.AppendAdditionalInfo(addInfo)

    x_min = h1.GetXaxis().GetXmin()
    x_max = h1.GetXaxis().GetXmax()
    y_min = h1.GetMinimum()
    if y_min !=0: y_min = 0
    y_max = h1.GetMaximum()
    y_max = y_max + 0.35 * (y_max - y_min)
    x_axis_name = h1.GetXaxis().GetTitle()
    canv = CMS.cmsCanvas(canv_name,x_min,x_max, y_min ,y_max,x_axis_name,ytitle,square=True, iPos=iPos)

    leg = CMS.cmsLeg(0.55, 0.79, 0.87, 0.89, textSize=0.04)
    if type(h)==list:
        CMS.cmsDraw(h[0], "P", marker= 8, lcolor =  ROOT.TColor.GetColor("#e42536"), mcolor =  ROOT.TColor.GetColor("#e42536"))
        leg.AddEntry(h[0], "reweighted sample", "p")
        CMS.cmsDraw(h[1], "P", marker= 8, lcolor = ROOT.kBlack, mcolor = ROOT.kBlack)
        leg.AddEntry(h[1], "Simulated sample", "p")
        leg.Draw("same")
    else:
        if ytitle=="SF":
            CMS.cmsDraw(h1, "P", marker= 8 , mcolor = fillcolor)
        else:
            CMS.cmsDraw(h1, "P", marker= 8 ,lcolor = fillcolor, mcolor = fillcolor)
    CMS.SaveCanvas(canv, folder+canv_name+".pdf")
    # CMS.SaveCanvas(canv, folder+"/png/"+canv_name+".png")
    # CMS.SaveCanvas(canv, folder+"/C/"+canv_name+".C")

def plot_withratio(h, folder, fillcolor, canv_name = "canv" ,extraTest="Preliminary", iPos=0, energy=13.6, lumi = 1,  addInfo="", ytitle = "Events"):
    CMS.SetExtraText(extraTest)
    iPos = iPos
    canv_name = canv_name
    CMS.SetLumi(lumi,"pb","",1)
    CMS.SetEnergy(energy)
    CMS.ResetAdditionalInfo()
    CMS.AppendAdditionalInfo(addInfo)
    if type(h)==list:
        h1 = h[0]
        # hist_dict = [k.GetName() for k in h]
    else:
        h1 = h
    x_min = h1.GetXaxis().GetXmin()
    x_max = h1.GetXaxis().GetXmax()
    y_min = h1.GetMinimum()
    if y_min !=0: y_min = 0
    y_max = h1.GetMaximum()
    y_max = y_max + 0.35 * (y_max - y_min)
    x_axis_name = h1.GetXaxis().GetTitle()

    dicanv = CMS.cmsDiCanvas(canv_name, x_min, x_max, y_min, y_max, 0, 2, x_axis_name, ytitle, "ratio", square =True, iPos=iPos)
    pad1 = dicanv.cd(1)
    leg = CMS.cmsLeg(0.55, 0.79, 0.87, 0.89, textSize=0.04)
    if type(h)==list:
        CMS.cmsDraw(h[0], "P", marker= 8, lcolor =  ROOT.TColor.GetColor("#e42536"), mcolor =  ROOT.TColor.GetColor("#e42536"))
        leg.AddEntry(h[0], "reweighted sample", "p")
        CMS.cmsDraw(h[1], "P", marker= 8, lcolor = ROOT.kBlack, mcolor = ROOT.kBlack)
        leg.AddEntry(h[1], "Simulated sample", "p")
        leg.Draw("same")
    else:
        pass
    pad2 = dicanv.cd(2)
    h_ratio = h[0].Clone()
    h_ratio.Divide(h[1])
    CMS.cmsDraw(h_ratio, "esame")
    ref_line = ROOT.TLine(x_min, 1, x_max, 1)
    CMS.cmsDrawLine(ref_line, lcolor=ROOT.kBlack, lstyle=ROOT.kDotted)
    CMS.SaveCanvas(dicanv, folder+canv_name+".pdf")

var = "mhh"
inputs = "1"
file_rw = ROOT.TFile.Open(f"reweighted_sample_GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p35_{inputs}.root")
h1 = file_rw.Get(f"{var}_weightedLHE")
file_or = ROOT.TFile.Open("original_sample_GluGlutoHHto2B2Tau_kl_1p00_kt_1p00_c2_0p35.root")
h2 = file_or.Get(f"{var}_LHE_target")
h1.Scale(1.0, "width")
h2.Scale(1.0, "width")


hist_ = [h1, h2]

plot_withratio(hist_, "./", ROOT.kBlack, f"distr_{var}_{inputs}", "Preliminary", 11, 13.6, 1, "(k_l=1.0, k_t=1.0, c_2=0.35)", ytitle = "Events / bin width [GeV]^{-1}")