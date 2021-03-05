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
#UL=2017
UL2017=sys.argv[2]
UL2018=sys.argv[3]

if(outHistFile=="JetpTresponse_eta0pt5to1pt3") : histtiltle="Jet pT response, 0.5<|#eta|<1.3"
if(outHistFile=="JetpTresponse_eta1pt3to2pt1") : histtiltle="Jet pT response, 1.3<|#eta|<2.1"
if(outHistFile=="JetpTresponse_eta2pt1to2pt5") : histtiltle="Jet pT response, 2.1<|#eta|<2.5"
if(outHistFile=="JetpTresponse_eta2pt5to3pt0") : histtiltle="Jet pT response, 2.5<|#eta|<3.0"
if(outHistFile=="JetpTresponse_eta") : histtiltle="Jet Response wrt Eta"
#events = Events ("step3_unity_offset0.root")
#events = Events ("../step3_unity_"+str(UL)+"_tmp.root")
#events = Events ("step3_UL2018.root")
events = Events ("../step3_"+str(UL2017)+".root")
events_ = Events ("../step3_"+str(UL2018)+".root")

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

CaloPF_JetPtHist2 = ROOT.TH1F ("CaloPF_JetPt2",  "1D "+histtiltle+ " ,for 100<jet pT<110", 100, 0, 2)
CaloPF_JetPtHist2_ = ROOT.TH1F ("CaloPF_JetPt2_", "1D "+histtiltle+" , for 100<jet pT<100", 100, 0, 2)

CaloPF_JetPt_EtaHist = ROOT.TProfile("CaloPF_JetPt_Eta",histtiltle, 30, 0, 3 )
CaloPF_JetPt_EtaHist_ = ROOT.TProfile("CaloPF_JetPt_Eta_",histtiltle, 30, 0, 3 )
line = ROOT.TLine(5,1 ,1000,1)
line1 = ROOT.TLine(5,1 ,1000,1)
line2 = ROOT.TLine(0,1 ,5,1)
line3 = ROOT.TLine(0,1 ,5,1)
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
            if(outHistFile=="JetpTresponse_eta0pt5to1pt3" and deltaR2<0.04 and abs(jjets.eta())<1.3 and abs(jjets.eta())>0.5):
                CaloPF_JetPtHist.Fill(jjets.pt(),ijets.pt()/jjets.pt())
                if(jjets.pt()>100 and jjets.pt()<110):
                    CaloPF_JetPtHist2.Fill(ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta1pt3to2pt1"and deltaR2<0.04 and abs(jjets.eta())<2.1 and abs(jjets.eta())>1.3):
                CaloPF_JetPtHist.Fill(jjets.pt(),ijets.pt()/jjets.pt())
                if(jjets.pt()>100 and jjets.pt()<110):
                    CaloPF_JetPtHist2.Fill(ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta2pt1to2pt5"and deltaR2<0.04 and abs(jjets.eta())<2.5 and abs(jjets.eta())>2.1):
                CaloPF_JetPtHist.Fill(jjets.pt(),ijets.pt()/jjets.pt())
                if(jjets.pt()>100 and jjets.pt()<110):
                    CaloPF_JetPtHist2.Fill(ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta2pt5to3pt0"and deltaR2<0.04 and abs(jjets.eta())<3.0 and abs(jjets.eta())>2.5):
                CaloPF_JetPtHist.Fill(jjets.pt(),ijets.pt()/jjets.pt())
                if(jjets.pt()>100 and jjets.pt()<110):
                    CaloPF_JetPtHist2.Fill(ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta" and deltaR2<0.04):
                CaloPF_JetPt_EtaHist.Fill(jjets.eta(),ijets.pt()/jjets.pt())

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
            if(outHistFile=="JetpTresponse_eta0pt5to1pt3" and deltaR2<0.04 and abs(jjets.eta())<1.3 and abs(jjets.eta())>0.5):
                CaloPF_JetPtHist_.Fill(jjets.pt(),ijets.pt()/jjets.pt())
                if(jjets.pt()>100 and jjets.pt()<110):
                    CaloPF_JetPtHist2_.Fill(ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta1pt3to2pt1"and deltaR2<0.04 and abs(jjets.eta())<2.1 and abs(jjets.eta())>1.3):
                CaloPF_JetPtHist_.Fill(jjets.pt(),ijets.pt()/jjets.pt())
                if(jjets.pt()>100 and jjets.pt()<110):
                    CaloPF_JetPtHist2_.Fill(ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta2pt1to2pt5"and deltaR2<0.04 and abs(jjets.eta())<2.5 and abs(jjets.eta())>2.1):
                CaloPF_JetPtHist_.Fill(jjets.pt(),ijets.pt()/jjets.pt())
                if(jjets.pt()>100 and jjets.pt()<110):
                    CaloPF_JetPtHist2_.Fill(ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta2pt5to3pt0"and deltaR2<0.04 and abs(jjets.eta())<3.0 and abs(jjets.eta())>2.5):
                CaloPF_JetPtHist_.Fill(jjets.pt(),ijets.pt()/jjets.pt())
                if(jjets.pt()>100 and jjets.pt()<110):
                    CaloPF_JetPtHist2_.Fill(ijets.pt()/jjets.pt())
            if(outHistFile=="JetpTresponse_eta" and deltaR2<0.04):
                CaloPF_JetPt_EtaHist_.Fill(jjets.eta(),ijets.pt()/jjets.pt())

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
CaloPF_JetPtHist.SetTitleSize (0.08)

CaloPF_JetPt_EtaHist.SetStats(0)
CaloPF_JetPt_EtaHist_.SetStats(0)
CaloPF_JetPt_EtaHist_.GetYaxis().SetRangeUser(0,1.5)
CaloPF_JetPt_EtaHist.GetYaxis().SetRangeUser(0,1.5)
CaloPF_JetPt_EtaHist_.SetLineColor(ROOT.kBlue)
CaloPF_JetPt_EtaHist.SetLineColor(ROOT.kRed)
CaloPF_JetPt_EtaHist_.SetLineWidth(2)
CaloPF_JetPt_EtaHist.SetLineWidth(2)
CaloPF_JetPt_EtaHist.GetYaxis().SetTitle ("PFJet pT/GenJet pT")
CaloPF_JetPt_EtaHist.GetYaxis().SetLabelSize(0.05)
CaloPF_JetPt_EtaHist.GetYaxis().SetTitleSize (0.05)
CaloPF_JetPt_EtaHist.GetYaxis().SetTitleOffset (0.86)
CaloPF_JetPt_EtaHist_.GetYaxis().SetTitle ("PFJet pT/GenJet pT")
CaloPF_JetPt_EtaHist.SetTitleSize(0.05) 
line.SetLineColor(ROOT.kBlack)
line.SetLineWidth(2)
line1.SetLineColor(ROOT.kBlack)
line1.SetLineWidth(2)

legend1 = ROOT . TLegend (0.6 ,0.3 ,0.8 ,0.5)
legend1.SetHeader("After Calibration","C")
legend1.AddEntry ( CaloPF_JetPtHist2 ,"Ultralegacy "+str(UL2017))
legend1.AddEntry ( CaloPF_JetPtHist2_ ,"Ultralegacy "+str(UL2018))
legend1.SetLineWidth (0)
legend1.SetTextSize(0.04)

CaloPF_JetPtHist2.SetStats(0)
CaloPF_JetPtHist2_.SetStats(0)
# CaloPF_JetPtHist2_.GetYaxis().SetRangeUser(0,1.5)
# CaloPF_JetPtHist2.GetYaxis().SetRangeUser(0,1.5)
CaloPF_JetPtHist2_.SetLineColor(ROOT.kBlue)
CaloPF_JetPtHist2.SetLineColor(ROOT.kRed)
CaloPF_JetPtHist2_.SetLineWidth(2)
CaloPF_JetPtHist2.SetLineWidth(2)
CaloPF_JetPtHist2.GetXaxis().SetTitle ("PFJet pT/GenJet pT")
CaloPF_JetPtHist2.SetTitleSize (0.05)
CaloPF_JetPtHist2.GetXaxis().SetLabelSize(0.05)
CaloPF_JetPtHist2.GetXaxis().SetTitleSize (0.05)
CaloPF_JetPtHist2.GetXaxis().SetTitleOffset (0.86)
CaloPF_JetPtHist2.GetYaxis().SetTitle ("Entries")

# c2 = ROOT.TCanvas()
# CaloPF_JetPtHist2_.Draw()
# c2.Print ("1DCaloPF_JetPtHist_.png")
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


legend = ROOT . TLegend (0.5 ,0.1 ,0.7 ,0.3)
legend.SetHeader("After Calibration","C")
legend.AddEntry ( CaloPF_JetPtHist ,"Ultralegacy "+str(UL2017))
legend.AddEntry ( CaloPF_JetPtHist_ ,"Ultralegacy "+str(UL2018))
legend.SetLineWidth (0)
legend.SetTextSize(0.06)

latex = ROOT. TLatex ()
latex . SetNDC ()
latex . SetTextSize (0.02)
latex . DrawText (0.7 ,0.83 , " Ultralegacy 2017 ")



if(outHistFile=="JetpTresponse_eta0pt5to1pt3" or outHistFile=="JetpTresponse_eta1pt3to2pt1" or outHistFile=="JetpTresponse_eta2pt1to2pt5" or outHistFile=="JetpTresponse_eta2pt5to3pt0"):
    c1 = ROOT.TCanvas()
    CaloPF_JetPtHist.SetStats(0)
    CaloPF_JetPtHist_.SetStats(0)
    CaloPF_JetPtHist2.Draw()
    CaloPF_JetPtHist2_.Draw("same")
    legend1.Draw("sames")
    c1.Print ("1DCaloPF_JetPtHist"+"_"+outHistFile+"_"+str(UL2017)+"wrt"+str(UL2018)+".png")
    c1.SaveAs ("1DCaloPF_JetPtHist"+"_"+outHistFile+"_"+str(UL2017)+"wrt"+str(UL2018)+".root")

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
    if(outHistFile=="JetpTresponse_eta2pt5to3pt0"): 
        ratio.GetYaxis().SetRangeUser(0.8,1.1)
    else : 
        ratio.GetYaxis().SetRangeUser(0.9,1.1)
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
    c3.Print(outHistFile+"_"+str(UL2017)+"wrt"+str(UL2018)+".png")
    c3.Print(outHistFile+"_"+str(UL2017)+"wrt"+str(UL2018)+".pdf")
    c3.SaveAs(outHistFile+"_"+str(UL2017)+"wrt"+str(UL2018)+".root")
# ratio.Write()
# outHistFile.Close()

if(outHistFile=="JetpTresponse_eta"):
    ratio = CaloPF_JetPt_EtaHist.Clone()
    ratio1 = CaloPF_JetPt_EtaHist_.Clone()
    ratio.Divide(ratio1)
    ratio.SetLineColor(ROOT.kWhite)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.7)
    ratio.SetMarkerColor(ROOT.kRed)
    ratio.SetTitle("")
    ratio.GetXaxis().SetLabelSize (0.12)
    ratio.GetXaxis().SetTitleSize (0.12)
    ratio.GetYaxis().SetRangeUser(0.9,1.1)
    ratio.GetYaxis().SetLabelSize (0.13)
    ratio.GetYaxis().SetTitleSize (0.15)
    ratio.GetYaxis().SetTitle ("Ratio")
    ratio.GetXaxis().SetTitle ("GenJet |#eta|")
    ratio.GetYaxis().SetTitleOffset (0.3)
    ratio.GetYaxis().SetNdivisions(207)
    
    c4= ROOT.TCanvas() 
    pad3 = ROOT . TPad (" pad3 "," pad3 " ,0 ,0.3 ,1 ,1)
    pad3 . SetBottomMargin (0)
    #pad3 . SetLogx ( True )
    pad3 . Draw ()
    pad3.cd()
    CaloPF_JetPt_EtaHist.Draw("pe")
    CaloPF_JetPt_EtaHist_.Draw("pe, same")
    line2.Draw("same")
    latex.Draw(" same ")
    legend.Draw (" same ")
    c4.cd()
    pad4 = ROOT . TPad (" pad4 "," pad4 " ,0 ,0.05 ,1 ,0.3)
    pad4 . SetTopMargin (0)
    pad4 . SetBottomMargin (0.25)
    #pad4 . SetLogx ( True )
    pad4.Draw()
    pad4.cd()
    ratio.Draw("pe")
    line3.Draw("same")
    c4.Print(outHistFile+"_"+str(UL2017)+"wrt"+str(UL2018)+".png")
    c4.Print(outHistFile+"_"+str(UL2017)+"wrt"+str(UL2018)+".pdf")
    c4.SaveAs(outHistFile+"_"+str(UL2017)+"wrt"+str(UL2018)+".root")
