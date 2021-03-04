import ROOT
import sys
import array
from DataFormats.FWLite import Events, Handle
def deltaPhi(phi1, phi2):
    pi=3.14;
    dphi = phi2-phi1
    if  dphi > pi:
        dphi -= 2.0*pi
    if dphi <= -pi:
        dphi += 2.0*pi
    return abs(dphi)
outHistFile=sys.argv[1]
#outHistFile="JetpTresponse_eta0pt5to1pt3"
#outHistFile="JetpTresponse_eta1pt3to2pt1"
#outHistFile="JetpTresponse_eta2pt1to2pt5"
#outHistFile="JetpTresponse_eta2pt5to3pt0"
if(outHistFile=="JetpTresponse_eta0pt5to1pt3") : histtiltle="Jet pT response, 0.5<|#eta|<1.3"
if(outHistFile=="JetpTresponse_eta1pt3to2pt1") : histtiltle="Jet pT response, 1.3<|#eta|<2.1"
if(outHistFile=="JetpTresponse_eta2pt1to2pt5") : histtiltle="Jet pT response, 2.1<|#eta|<2.5"
if(outHistFile=="JetpTresponse_eta2pt5to3pt0") : histtiltle="Jet pT response, 2.5<|#eta|<3.0"
#events = Events ("step3_unity_offset0.root")
events = Events ("step3_unity_nooffset.root")
events_ = Events ("step3_UL2017.root")

# create handle outside of loop
handle  = Handle ("std::vector<reco::PFJet>")
handle1  = Handle ("std::vector<reco::GenJet>")

# for now, label is just a tuple of strings that is initialized just
# like and edm::InputTag
label = ("ak4PFJets")
label1 = ("ak4GenJets")
pt_thresholds = [ 10**(x/10.) for x in range(7,31) ]
# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background
#zmassHist = ROOT.TH1F ("zmass", "Z Candidate Mass", 50, 20, 220)
PFJetPtHist = ROOT.TH1F ("PFJetPt", "AK4 PF Jets Pt", 500, 0, 500)
CaloJetPtHist = ROOT.TH1F ("CaloJetPt", "AK4 Gen Jets Pt", 500, 0, 500)
CaloPF_JetPtHist = ROOT.TProfile("CaloPF_JetPt",histtiltle,len(pt_thresholds)-1, array.array('d', pt_thresholds))
CaloPF_JetPtHist_ = ROOT.TProfile("CaloPF_JetPt_", "UL2017 response",len(pt_thresholds)-1, array.array('d', pt_thresholds))
line = ROOT.TLine(5,1 ,1000,1)
line1 = ROOT.TLine(5,1 ,1000,1)
#loop over events
for event in events:
    # use getByLabel, just like in cmsRun
    event.getByLabel (label, handle)
    event.getByLabel (label1, handle1)
    # get the product
    jets = handle.product()
    numjets = len (jets)
    for i in xrange (numjets - 1):
        ijets=jets[i]
        PFJetPtHist.Fill(ijets.pt())
    # print("no. of Ak4Jets  : ")
    # print(numjets)
    # use getByLabel, just like in cmsRun
    jet1s = handle1.product()
    numjet1s = len (jet1s)
    for i in xrange (numjet1s - 1):
        ijets=jet1s[i]
        CaloJetPtHist.Fill(ijets.pt())
    # print("no. of CaloJets  : ") 
    # print(numjet1s)
    # print("-=======-")

    for i in xrange (numjets - 1):
        for j in xrange (numjet1s - 1):
            ijets=jets[i]
            jjets=jet1s[j]
            deltaR2=deltaPhi(ijets.phi(), jjets.phi())**2 + (ijets.eta()- jjets.eta())**2
            if(outHistFile=="JetpTresponse_eta0pt5to1pt3" and deltaR2<0.04 and abs(ijets.eta())<1.3 and abs(jjets.eta())>0.5):
                CaloPF_JetPtHist.Fill(jjets.pt(),ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta1pt3to2pt1"and deltaR2<0.04 and abs(ijets.eta())<2.1 and abs(jjets.eta())>1.3):
                CaloPF_JetPtHist.Fill(jjets.pt(),ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta2pt1to2pt5"and deltaR2<0.04 and abs(ijets.eta())<2.5 and abs(jjets.eta())>2.1):
                CaloPF_JetPtHist.Fill(jjets.pt(),ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta2pt5to3pt0"and deltaR2<0.04 and abs(ijets.eta())<3.0 and abs(jjets.eta())>2.5):
                CaloPF_JetPtHist.Fill(jjets.pt(),ijets.pt()/jjets.pt())

for event in events_:
    # use getByLabel, just like in cmsRun
    event.getByLabel (label, handle)
    event.getByLabel (label1, handle1)
    # get the product
    jets = handle.product()
    numjets = len (jets)
    jet1s = handle1.product()
    numjet1s = len (jet1s)
    for i in xrange (numjets - 1):
        for j in xrange (numjet1s - 1):
            ijets=jets[i]
            jjets=jet1s[j]
            deltaR2=deltaPhi(ijets.phi(), jjets.phi())**2 + (ijets.eta()- jjets.eta())**2
            if(outHistFile=="JetpTresponse_eta0pt5to1pt3" and deltaR2<0.04 and abs(ijets.eta())<1.3 and abs(jjets.eta())>0.5):
                CaloPF_JetPtHist_.Fill(jjets.pt(),ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta1pt3to2pt1"and deltaR2<0.04 and abs(ijets.eta())<2.1 and abs(jjets.eta())>1.3):
                CaloPF_JetPtHist_.Fill(jjets.pt(),ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta2pt1to2pt5"and deltaR2<0.04 and abs(ijets.eta())<2.5 and abs(jjets.eta())>2.1):
                CaloPF_JetPtHist_.Fill(jjets.pt(),ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta2pt5to3pt0"and deltaR2<0.04 and abs(ijets.eta())<3.0 and abs(jjets.eta())>2.5):
                CaloPF_JetPtHist_.Fill(jjets.pt(),ijets.pt()/jjets.pt())

# make a canvas, draw, and save it


CaloPF_JetPtHist.SetStats(0)
CaloPF_JetPtHist_.SetStats(0)
CaloPF_JetPtHist_.GetYaxis().SetRangeUser(0,1.5)
CaloPF_JetPtHist.GetYaxis().SetRangeUser(0,1.5)
CaloPF_JetPtHist_.SetLineColor(ROOT.kBlue)
CaloPF_JetPtHist.SetLineColor(ROOT.kRed)
CaloPF_JetPtHist_.SetLineWidth(2)
CaloPF_JetPtHist.SetLineWidth(2)
CaloPF_JetPtHist.GetYaxis().SetTitle ("PFJet pT/GenJet pT")
CaloPF_JetPtHist.GetYaxis().SetLabelSize(0.05)
CaloPF_JetPtHist.GetYaxis().SetTitleSize (0.05)
CaloPF_JetPtHist.GetYaxis().SetTitleOffset (0.86)
CaloPF_JetPtHist_.GetYaxis().SetTitle ("PFJet pT/GenJet pT")
#CaloPF_JetPtHist.GetXaxis().SetTitleSize(0.05) 
line.SetLineColor(ROOT.kBlack)
line.SetLineWidth(2)
line1.SetLineColor(ROOT.kBlack)
line1.SetLineWidth(2)
# c1 = ROOT.TCanvas()
# PFJetPtHist.Draw()
# c1.Print ("PFJetPt.png")
# c2 = ROOT.TCanvas()
# CaloJetPtHist.Draw()
# c2.Print ("CaloJetPt.png")
# # c3 = ROOT.TCanvas()
# # CaloPF_JetPtHist.Draw()
# # c3.SetLogx ( True )
# # c3.Print ("CaloPF_JetPt.png")
# c4 = ROOT.TCanvas()
# CaloPF_JetPtHist.Draw()                                                                                                       
# CaloPF_JetPtHist_.Draw("same")
# line.Draw("same")
# c4.SetLogx ( True )
# c4.Print ("CaloPF_JetPt_UL.png")

legend = ROOT . TLegend (0.5 ,0.2 ,0.7 ,0.5)
legend.SetHeader("Ultralegacy 2017","C")
legend.AddEntry ( CaloPF_JetPtHist ,"Before PFCalibration")
legend.AddEntry ( CaloPF_JetPtHist_ ,"After PFCalibration ")
legend.SetLineWidth (0)
legend.SetTextSize(0.04)

latex = ROOT. TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.02)
latex . DrawText (0.7 ,0.83 , " Ultralegacy 2017 ")


ratio = CaloPF_JetPtHist.Clone()
ratio1 = CaloPF_JetPtHist_.Clone()
ratio.Divide(ratio1)
ratio.SetLineColor(ROOT.kWhite)
ratio.SetMarkerStyle(20)
ratio.SetMarkerSize(0.7)
ratio.SetMarkerColor(ROOT.kRed)
ratio.SetTitle("")
ratio.GetXaxis().SetLabelSize (0.12)
ratio.GetXaxis().SetTitleSize (0.12)
ratio.GetYaxis().SetRangeUser(0.8,1.2)
ratio.GetYaxis().SetLabelSize (0.13)
ratio.GetYaxis().SetTitleSize (0.15)
ratio.GetYaxis().SetTitle ("Ratio")
ratio.GetXaxis().SetTitle ("GenJet pT")
ratio.GetYaxis().SetTitleOffset (0.3)
ratio.GetYaxis().SetNdivisions(207)

c3= ROOT.TCanvas() 
pad1 = ROOT . TPad (" pad1 "," pad1 " ,0 ,0.3 ,1 ,1)
pad1 . SetBottomMargin (0)
pad1 . SetLogx ( True )
pad1 . Draw ()
pad1.cd()
CaloPF_JetPtHist.Draw("pe")
CaloPF_JetPtHist_.Draw("pe, same")
line.Draw("same")
latex.Draw(" same ")
legend.Draw (" same ")
c3.cd()
pad2 = ROOT . TPad (" pad2 "," pad2 " ,0 ,0.05 ,1 ,0.3)
pad2 . SetTopMargin (0)
pad2 . SetBottomMargin (0.25)
pad2 . SetLogx ( True )
pad2.Draw()
pad2.cd()
ratio.Draw("pe")
line.Draw("same")
c3.Print(outHistFile+".png")
c3.SaveAs(outHistFile+".root")
# ratio.Write()
# outHistFile.Close()
