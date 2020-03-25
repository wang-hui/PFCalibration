void bEtaEndcapH()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();
  TF1 *fbEtaEndcapH = new TF1("fbEtaEndcapH","[0]+[1]*exp(-x/[2])",1.,500.);
  fbEtaEndcapH->FixParameter(0,0.022428);
  fbEtaEndcapH->FixParameter(1,0.100791);
  fbEtaEndcapH->FixParameter(2,85.7951);
  TF1 *fbEtaEndcapH2 = new TF1("fbEtaEndcapH2","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",1.,1000.);
  fbEtaEndcapH2->FixParameter(0,0.0214894);
  fbEtaEndcapH2->FixParameter(1,-0.266704);
  fbEtaEndcapH2->FixParameter(2,5.2112);
  fbEtaEndcapH2->FixParameter(3,0.303578);
  fbEtaEndcapH2->FixParameter(4,-104.367);

  fbEtaEndcapH->GetYaxis()->SetRangeUser(-1, 1);
  fbEtaEndcapH->GetXaxis()->SetRangeUser(2., 500.);

  fbEtaEndcapH->SetLineColor(1);
   fbEtaEndcapH->Draw();
  fbEtaEndcapH2->Draw("sames");

  TLegend* legends = new TLegend(0.45, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(fbEtaEndcapH, "for UL2017 (chi2/Ndf = 15.3333)","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                                     
  legends->AddEntry(fbEtaEndcapH2,"for old param (chi2/Ndf = 33.0994)","l");                                  
legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("fbEtaEndcapH.gif");                                                     
}
