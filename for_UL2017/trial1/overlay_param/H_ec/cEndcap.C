void cEndcap()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();
  TF1 *  fcEndcap = new TF1("fcEndcap","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",2.,500.);
  fcEndcap->FixParameter(0,1.13985);
  fcEndcap->FixParameter(1,0.196757);
  fcEndcap->FixParameter(2,-2.56203);
  fcEndcap->FixParameter(3,22.1335);
  fcEndcap->FixParameter(4,0.231671);
  fcEndcap->FixParameter(5,0.0625535);
  fcEndcap->FixParameter(6,1.71307);
  fcEndcap->FixParameter(7,-0.987053);

  TF1 *fcEndcap2 = new TF1("fcEndcap2","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",2.,500.);
  fcEndcap2->FixParameter(0,1.01863);
  fcEndcap2->FixParameter(1,1.29787);
  fcEndcap2->FixParameter(2,-3.97293);
  fcEndcap2->FixParameter(3,21.7805);
  fcEndcap2->FixParameter(4,0.810195);
  fcEndcap2->FixParameter(5,0.234134);
  fcEndcap2->FixParameter(6,1.42226);
  fcEndcap2->FixParameter(7,-0.0997326);

  fcEndcap->SetLineColor(1);
  fcEndcap->Draw();
  fcEndcap2->Draw("sames");

  TLegend* legends = new TLegend(0.45, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(fcEndcap, "for UL2017 (chi2/Ndf = 1.32901 )","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                                     
  legends->AddEntry(fcEndcap2,"for old param (chi2/Ndf = 33.8286)","l");                                  
legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("fcEndcap.gif");                                                     
}
