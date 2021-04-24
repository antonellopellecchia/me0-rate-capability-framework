import ROOT as rt

rt.gROOT.SetBatch(rt.kTRUE)

def drawCmsPreliminaryLabel(x, y):
    latex = rt.TLatex()
    latex.SetTextSize(.055)
    latex.DrawLatexNDC(x, y, '#bf{CMS}')
    latex.SetTextSize(.04)
    latex.DrawLatexNDC(x, y, '          #it{Preliminary}')

def getStyle():
    style = rt.TStyle("CMS", "CMS Style")
    style.Reset()
    style.SetFillColor(1)
    style.SetFillStyle(1001)
    style.SetCanvasBorderMode(0)
    style.SetCanvasColor(0)
    style.SetCanvasPreferGL(rt.kTRUE)
    style.SetCanvasDefH(600)
    style.SetCanvasDefW(700)
    #style.SetPadBorderMode(1)
    style.SetPadColor(0)
    style.SetPadLeftMargin(0.15)
    style.SetPadBottomMargin(0.15)
    style.SetPadRightMargin(0.075)
    style.SetPadTopMargin(0.06)
    style.SetPadTickX(0)
    style.SetPadTickY(0)
    style.SetFrameFillColor(0)
    #style.SetFrameBorderMode(1)
    #style.SetDrawBorder(0)
    style.SetLegendBorderSize(0)

    style.SetGridColor(rt.kGray)
    style.SetGridStyle(3)
    style.SetGridWidth(1)
    style.SetPadGridX(rt.kFALSE)
    style.SetPadGridY(rt.kFALSE)
    
    font = 42
    tsize = 0.05
    style.SetTextFont(font)
    style.SetTextSize(tsize)
    style.SetTitleStyle(0)
    style.SetTitleBorderSize(0)
    style.SetTitleColor(1, "xyz")
    style.SetTitleColor(1, "t")
    style.SetTitleFillColor(0)
    style.SetTitleFont(font, "xyz")
    style.SetTitleFont(font, "t")
    style.SetTitleOffset(1.2, "xyz")
    style.SetTitleOffset(1.5, "y")
    style.SetTitleSize(tsize, "xyz")
    style.SetTitleSize(tsize, "t")

    style.SetLegendFont(font)
    style.SetStatStyle(0)
    style.SetStatBorderSize(0)
    style.SetStatColor(0)
    style.SetStatFont(font)
    style.SetStatFontSize(tsize)
    style.SetStatX(0.58)
    style.SetStatY(0.88)
    style.SetStatW(0.2)
    style.SetStatH(0.1)
    #style.SetOptStat(111110)
    style.SetOptStat(0)
    style.SetOptFit(1)
    style.SetStatFormat("6.3g")
    style.SetLabelFont(font, "xyz")
    style.SetLabelSize(0.04, "xyz")
    style.SetLabelOffset(0.01, "xyz")
    style.SetOptTitle(0)
    style.SetPaperSize(rt.TStyle.kA4)
    style.SetFuncWidth(2)
    style.SetHistLineColor(rt.kBlue - 3)
    style.SetPalette(1)
    style.SetAxisColor(rt.kBlack, "X")
    style.SetAxisColor(rt.kBlack, "Y")
    style.SetAxisColor(rt.kBlack, "Z")
    style.SetNdivisions(505, "x")
    style.SetNdivisions(510, "y")
    lw = 1
    style.SetLineWidth(lw)
    style.SetLineColor(rt.kBlack)
    #style.SetLineStyle(7)
    #style.SetLineStyleString(2, "[12 12]")
    style.SetFrameLineWidth(lw)
    style.SetHistLineWidth(lw)
    style.SetFuncWidth(lw)
    style.SetFuncColor(rt.kRed-2)
    style.SetGridWidth(lw) 
    style.SetMarkerSize(1.0)
    style.SetMarkerStyle(8)
    style.SetMarkerColor(rt.kBlue-1)
    style.cd()

    #rt.TGaxis.SetMaxDigits(3)
    
    return style

style = getStyle()
style.cd()
rt.gROOT.SetStyle("CMS")
rt.gROOT.ForceStyle()
