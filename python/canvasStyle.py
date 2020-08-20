import ROOT
def canv(c):
    #c = ROOT.TCanvas('c', '',1920,1080)
    
    #For the axis titles:
    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.045, "XYZ")
    
    ROOT.gStyle.SetTitleXOffset(1.2)
    ROOT.gStyle.SetTitleYOffset(0.9)
    
    #For the axis labels:
    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.005, "XYZ")
    ROOT.gStyle.SetLabelSize(0.035, "XYZ")
    
    #For the axis:
    ROOT.gStyle.SetAxisColor(1, "XYZ")
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.02, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetPadTickX(1)                          #To get tick marks on the opposite side of the frame
    ROOT.gStyle.SetPadTickY(1)
    #c.SetLeftMargin(0.20)
    #c.SetTopMargin(0.15)
    c.SetBottomMargin(0.15)
    c.SetRightMargin(0.20)
    c.SetGrid()
