void bBarrel()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();

  TF1 *fbBarrel = new TF1("fbBarrel","([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))",1.,1000.);
  fbBarrel->FixParameter(0,1.74951);
  fbBarrel->FixParameter(1,0.317424);
  fbBarrel->FixParameter(2,-7.78609);
  fbBarrel->FixParameter(3,1.44197);
  fbBarrel->FixParameter(4,0.818853);
  fbBarrel->FixParameter(5,0.0734628);
  fbBarrel->FixParameter(6,0.135386);
  fbBarrel->FixParameter(7,-0.836217);
  TF1 *fbBarrel2 = new TF1("fbBarrel2","([0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5])))))",1.,1000.);
  fbBarrel2->FixParameter(0,2.25366);
  fbBarrel2->FixParameter(1,0.537715);
  fbBarrel2->FixParameter(2,-4.81375);
  fbBarrel2->FixParameter(3,12.109);
  fbBarrel2->FixParameter(4,1.80577);
  fbBarrel2->FixParameter(5,0.187919);
  fbBarrel2->FixParameter(6,-6.26234);
  fbBarrel2->FixParameter(7,-0.607392);

  fbBarrel->SetLineColor(1);
  fbBarrel->Draw();
  fbBarrel2->Draw("sames");

  TLegend* legends = new TLegend(0.5, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(fbBarrel, "for UL2017 (chi2/Ndf = 11.6524)","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                                   
  legends->AddEntry(fbBarrel2,"for current param (chi2/Ndf = 125.033)","l");                                  
legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("fbBarrel.gif");                                                     
}
