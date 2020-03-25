void cBarrel()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();
  TF1 *fcBarrel = new TF1("fcBarrel","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,500.);
  fcBarrel->FixParameter(0,3.00232);
  fcBarrel->FixParameter(1,3.96519);
  fcBarrel->FixParameter(2,-10.7502);
  fcBarrel->FixParameter(3,3.25745);
  fcBarrel->FixParameter(4,4.68851);
  fcBarrel->FixParameter(5,0.335698);
  fcBarrel->FixParameter(6,0.0300156);
  fcBarrel->FixParameter(7,-0.6611);

  TF1 *fcBarrel2 = new TF1("fcBarrel2","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,500.);
  fcBarrel2->FixParameter(0,1.5126);
  fcBarrel2->FixParameter(1,0.855057);
  fcBarrel2->FixParameter(2,-6.04199);
  fcBarrel2->FixParameter(3,2.08229);
  fcBarrel2->FixParameter(4,0.592266);
  fcBarrel2->FixParameter(5,0.0291232);
  fcBarrel2->FixParameter(6,0.364802);
  fcBarrel2->FixParameter(7,-1.50142);

  fcBarrel->SetLineColor(1);
  fcBarrel->Draw();
  fcBarrel2->Draw("sames");

  TLegend* legends = new TLegend(0.45, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(fcBarrel, "for UL2017 (chi2/Ndf = 23.557","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                                     
  legends->AddEntry(fcBarrel2,"for old param (chi2/Ndf = 18.9596)","l");                                  
legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("fcBarrel.gif");                                                     
}
