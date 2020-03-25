void bEtaBarrelEH()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();
  TF1 *fbEtaBarrelEH = new TF1("fbEtaBarrelEH","[0]+[1]*exp(-x/[2])",1.,500.);
  fbEtaBarrelEH->SetParameter(0,0.0557658);
  fbEtaBarrelEH->SetParameter(1,0.124518);
  fbEtaBarrelEH->SetParameter(2,204.727);
  TF1 *fbEtaBarrelEH2 = new TF1("fbEtaBarrelEH2","[0]+[1]*exp(-x/[2])",1.,500.);
  fbEtaBarrelEH2->SetParameter(0,0.0396458);
  fbEtaBarrelEH2->SetParameter(1,0.114128);
  fbEtaBarrelEH2->SetParameter(2,251.405);
  fbEtaBarrelEH->GetYaxis()->SetRangeUser(-1, 1);

  fbEtaBarrelEH->SetLineColor(1);
   fbEtaBarrelEH->Draw();
  fbEtaBarrelEH2->Draw("sames");

  TLegend* legends = new TLegend(0.45, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(fbEtaBarrelEH, "for UL2017 (chi2/Ndf = 12.3)","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                           \  
  legends->AddEntry(fbEtaBarrelEH2,"for old param (chi2/Ndf = 54.1033)","l");                                  
  legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("fbEtaBarrelEH.gif");                                                     
}
