void bEtaBarrelH()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();
  TF1 *fbEtaBarrelH = new TF1("fbEtaBarrelH","[0]+[1]*exp(-x/[2])",1.,500.);
  fbEtaBarrelH->SetParameter(0,0.000578309);
  fbEtaBarrelH->SetParameter(1,-0.491991);
  fbEtaBarrelH->SetParameter(2,1.445);

  TF1 *fbEtaBarrelH2 = new TF1("fbEtaBarrelH2","[0]+[1]*exp(-x/[2])",1.,500.);
  fbEtaBarrelH2->SetParameter(0,-0.0232604);
  fbEtaBarrelH2->SetParameter(1,0.0937525);
  fbEtaBarrelH2->SetParameter(2,34.9935);

  fbEtaBarrelH->GetYaxis()->SetRangeUser(-1, 1);

  fbEtaBarrelH->SetLineColor(1);
   fbEtaBarrelH->Draw();
  fbEtaBarrelH2->Draw("sames");

  TLegend* legends = new TLegend(0.45, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(fbEtaBarrelH, "for UL2017 (chi2/Ndf = 316.128)","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                                     
  legends->AddEntry(fbEtaBarrelH2,"for old param (chi2/Ndf = 46.5598)","l");                                  
legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("fbEtaBarrelH.gif");                                                     
}
