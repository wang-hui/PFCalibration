void aBarrel()
{
  TCanvas *c2 = new TCanvas("c2", "c2",65,52,525,527);
  c2->Range(-60.25,-0.625,562.25,0.625);
  c2->SetLogx();
  TF1 *faBarrel = new TF1("faBarrel","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[4]/[3]))))",1.,500.);
  faBarrel->FixParameter(0,-1.28945);
  faBarrel->FixParameter(1,2.29041);
  faBarrel->FixParameter(2,4.93862);
  faBarrel->FixParameter(3,0.0744576);
  faBarrel->FixParameter(4,-1.80063);

  TF1  *faBarrel2 = new TF1("faBarrel2","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,500.);
  faBarrel2->FixParameter(0,-30.7141);
  faBarrel2->FixParameter(1,31.7583);
  faBarrel2->FixParameter(2,4.40594);
  faBarrel2->FixParameter(3,1.70914);
  faBarrel2->FixParameter(4,0.0613696);
  faBarrel2->FixParameter(5,0.000104857);
  faBarrel2->FixParameter(6,-1.38927);
  faBarrel2->FixParameter(7,-0.743082);
  faBarrel->SetLineColor(1);
  faBarrel->Draw();
  faBarrel2->Draw("sames");

  // cout<<"Ndf for  faBarrel ;  "<<faBarrel->GetNDF()<<endl;                                                                                                                                   
  // cout<<"Chi Square for  faBarrel ;  "<<faBarrel->GetChisquare()<<endl;                                                                                                                      
  // cout<<"Chi Square/Ndf for  faBarrel ;  "<<","<<faBarrel->GetChisquare()/faBarrel->GetNDF()<<endl<<endl;                                                                     


  TLegend* legends = new TLegend(0.45, 0.3, 0.9, 0.4,"","brNDC"); // the numbers determine the position of the box         
  
  legends->SetFillColor(0);                                                                       
  //legends->SetHeader(stat_name);                                                                  
  // //   legends->SetTextSize(40);                                                                  
  legends->AddEntry(faBarrel, "for UL2017 ( chi2/ndf = 89.9565 )","l");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing )                                                     
  legends->AddEntry(faBarrel2,"for old calib ( chi2/ndf = 9.47149 )","l");                                  
legends->Draw();
 c2->Modified();
 c2->cd();
 c2->SaveAs("faBarrel.gif");                                                     
}
