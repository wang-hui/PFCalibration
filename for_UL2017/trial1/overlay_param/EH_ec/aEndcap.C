void aEndcap()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();
  // TF1 *  faEndcap = new TF1("faEndcap","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",2.,500.);
  // faEndcap->SetParameter(0,1.13418);
  // faEndcap->SetParameter(1,18.8699);
  // faEndcap->SetParameter(2,-41.7296);
  // faEndcap->SetParameter(3,0.814923);
  // faEndcap->SetParameter(4,0.000942355);
  // faEndcap->SetParameter(5,-0.0705279);
  // faEndcap->SetParameter(6,0.227606);
  // faEndcap->SetParameter(7,-0.573016);
  TF1 *faEndcap = new TF1("faEndcap","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[5]/[3])))-([4]))",2.,500.);
  faEndcap->SetParameter(0,28.9682);
  faEndcap->SetParameter(1,27.6333);
  faEndcap->SetParameter(2,-63.9996);
  faEndcap->SetParameter(3,0.670415);
  faEndcap->SetParameter(4,27.8404);
  faEndcap->SetParameter(5,0.20647);

  TF1 *faEndcap2 = new TF1("faEndcap2","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",2.,500.);
  faEndcap2->SetParameter(0,1.17227);
  faEndcap2->SetParameter(1,13.1489);
  faEndcap2->SetParameter(2,-29.1672);
  faEndcap2->SetParameter(3,0.604223);
  faEndcap2->SetParameter(4,0.0426363);
  faEndcap2->SetParameter(5,3.30898e-15);
  faEndcap2->SetParameter(6,0.165293);
  faEndcap2->SetParameter(7,-7.56786);


  faEndcap->SetLineColor(1);
  faEndcap->Draw();
  faEndcap2->Draw("sames");

  TLegend* legends = new TLegend(0.45, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(faEndcap, "for UL2017 (chi2/Ndf = 3.3819)","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                                     
  legends->AddEntry(faEndcap2,"for old param (chi2/Ndf = 2.85436)","l");                                  
legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("faEndcap.gif");                                                     
}
