void aEtaEndcapH()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();
  TF1 *faEtaEndcapH = new TF1("faEtaEndcapH","[0]+[1]*exp(-x/[2])",1.,1000.);
  faEtaEndcapH->FixParameter(0,-0.00267139);
  faEtaEndcapH->FixParameter(1,0.0376755);
  faEtaEndcapH->FixParameter(2,42.5313);

  TF1 *faEtaEndcapH2 = new TF1("faEtaEndcapH2","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",1.,1000.);
  faEtaEndcapH2->FixParameter(0,-0.0106029);
  faEtaEndcapH2->FixParameter(1,-0.692207);
  faEtaEndcapH2->FixParameter(2,0.0542991);
  faEtaEndcapH2->FixParameter(3,-0.171435);
  faEtaEndcapH2->FixParameter(4,-61.2277);

  faEtaEndcapH->GetYaxis()->SetRangeUser(-1, 1);
  faEtaEndcapH->GetXaxis()->SetRangeUser(2., 500.);

  faEtaEndcapH->SetLineColor(1);
   faEtaEndcapH->Draw();
  faEtaEndcapH2->Draw("sames");
 

  TLegend* legends = new TLegend(0.45, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(faEtaEndcapH, "for UL2017 (chi2/Ndf = 3.32866)","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                                     
  legends->AddEntry(faEtaEndcapH2,"for old param (chi2/Ndf = 43.9929)","l");                                  
legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("faEtaEndcapH.gif");                                                     
}
