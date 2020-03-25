#include <stdio.h>
// #include<conio.h>
void generate_1Dplots()
{
  char* hname = new char[200];
  char* hist_name = new char[200];
  char* stat_name = new char[200];
  char* hname1 = new char[200];
  char* path = new char[200];
  char* plot_name = new char[200];
  char* old_rel = new char[200];
  char* new_rel = new char[200];
  char* old_release = new char[50];
  char* new_release = new char[50];
  char* hist = new char[200];
  sprintf(old_release, "corrEta");
  sprintf(new_release, "corrEta");
  sprintf(plot_name, "EH_ec_out");
  sprintf(hist, "histcorhybrid");
  double xmin= -1., xmax =1.;
  double ymin= 0, y1, y2;
  int n[69];
  n[0]=0;
  n[1]=2;
  for(int i=2; i<69; i++)
    {
      if (n[i-1]<12)
	n[i]=n[i-1]+2;
      else if (n[i-1]>=12 && n[i-1]<104)
	n[i]=n[i-1]+4;
      else
	n[i]=n[i-1]+10;
    }

  for(int i=2; i<68;i++)
    { 
      sprintf(stat_name,"Etrue bin %d-%d GeV",n[i],n[i+1]);
      // cout<<"\n i "<<i<<" = "<<n[i];}
      // cout<<stat_name<<"\n";


      sprintf(hname1,"correta_%s_10_6_0_pre2.root",plot_name);
      sprintf(hname,"../../../trial_oldparam/1D_response/correta_%s_10_6_0_pre2.root",plot_name);
      // sprintf(plot_name,"%s",plot);
      sprintf(hist_name,"%s%d",hist,n[i]);
      //     sprintf(path,"../overlay_response_plots/withbefore15may/%s/%s.gif",plot_name,hist_name);
      sprintf(path,"../../1D_response_plots_pdf/withraw/corr_%s/%s.pdf",plot_name,hist_name);
      sprintf(old_rel,"%s",old_release);
      sprintf(new_rel,"%s",new_release);



      // sprintf(hname,"%s.root",plot);    
      // cout<<"1 \n";

      // sprintf(hname1,"%s.root",plot);
      // cout<<hname1;
      // sprintf(plot_name,"%s",plot);
      // sprintf(hist_name,"%s%d",hist,n[i]);
      // sprintf(path,"../1D_response_plots/%s/%s.gif",plot,hist_name);
      // sprintf(old_rel,"%s",old_release);
      // sprintf(new_rel,"%s",new_release);


   TCanvas *Canvas_1_n2 = new TCanvas(hist_name, hist_name,65,52,525,527);
   Canvas_1_n2->Range(-60.25,-0.625,562.25,0.625);
   Canvas_1_n2->SetFillColor(0);
   Canvas_1_n2->SetBorderMode(0);
   Canvas_1_n2->SetBorderSize(2);
   Canvas_1_n2->SetGridx();
   Canvas_1_n2->SetGridy();
   Canvas_1_n2->SetFrameBorderMode(0);
   Canvas_1_n2->SetFrameBorderMode(0);
   TFile * inputfile1 = new TFile(hname,"READ");
   TFile *inputfile2 = new TFile(hname1,"READ");
   TH1D *resp1 = (TH1D*) inputfile1 -> Get(hist_name);
   TH1D *resp2 = (TH1D*) inputfile2 -> Get(hist_name);

   // respHisto__1->SetStats(0);

   // Int_t ci;      // for color index setting
   // TColor *color; // for color definition with alpha
   // ci = TColor::GetColor("#000099");
   // respHisto__1->SetLineColor(ci);
   // respHisto__1->Draw("");
   TPaveStats *ptstats = new TPaveStats(0.65,0.67,0.9,0.9,"brNDC");
   // ptstats = new TPaveStats(0.6260745,0.403397,0.8810888,0.6433121,"brNDC");
   //   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetLineColor(4);
   ptstats->SetTextAlign(12);
   ptstats->SetTextColor(4);
   ptstats->SetTextFont(42);
   ptstats->AddText(hist_name);
   ptstats->SetTextSize(0.0315317);
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(10001);
   ptstats->Draw();
   resp2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(resp2);

   TPaveStats *ptstats2 = new TPaveStats(0.65,0.47,0.9,0.67,"brNDC");
   //   TPaveStats *ptstats = new TPaveStats(0.6260745,0.6433121,0.8810888,0.8832272,"brNDC");
   //  ptstats->SetName(old_rel);
   ptstats2->SetBorderSize(1);
   ptstats2->SetFillColor(0);
   ptstats2->SetLineColor(kRed);
   ptstats2->SetTextAlign(12);
   ptstats2->SetTextColor(kRed);
   ptstats2->SetTextFont(42);
   TText *ptstats_LaTex = ptstats2->AddText(old_rel);
   ptstats_LaTex->SetTextSize(0.0315317);
   ptstats2->SetOptStat(1111);
   ptstats2->SetOptFit(10001);
   ptstats2->Draw();
   resp1->GetListOfFunctions()->Add(ptstats2);
   ptstats2->SetParent(resp1);

   // TAxis* ffy = resp1->GetXaxis();
   // TAxis* ffr = resp2->GetXaxis();   
   // y1=ffy->GetYmax();
   // y2=ffr->GetYmax();
   y1=resp1->GetMaximum();
   y2=resp2->GetMaximum();
   if (y1>y2)
     ymax=y1+100;
   else
     ymax=y2+100;

   resp2->SetTitle(stat_name);
   resp2->SetLineColor(4);
   resp2->SetLineWidth(2);
   // resp2->Scale(1.0/resp2->Integral());
   resp2->GetXaxis()->SetRangeUser(xmin , xmax);                                                                                                                                                          
   resp2->GetYaxis()->SetRangeUser(ymin , ymax);
   resp2->GetXaxis()->SetTitle("Energy response");
   resp2->GetXaxis()->SetTitleSize(22);
   resp2->GetXaxis()->SetTitleFont(43);
   resp2->GetXaxis()->SetTitleOffset(0.8);  
   resp2->SetName("New Calibration");
   //   resp2->Fit("gaus");
   // TF1* ga1 = (TF1*)resp1->GetFunction("gaus");
   // ga1->Delete();
   //    // cout<<"2 \n";
   // TF1* ga2 = (TF1*)resp2->GetFunction("gaus");
   // ga2->Delete();

   resp2->Draw();

   resp1->SetLineColor(kRed);
   resp1->SetLineWidth(2);
   // resp1->Scale(1.0/resp1->Integral());
   //   resp1->GetXaxis()->SetRangeUser(xmin,xmax); 
   resp1->GetYaxis()->SetRangeUser(ymin , ymax);                                                                                                                                                            
   resp1->SetName("Old Calibration");
   resp1->Draw("sames");



   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);

   TF1* f1 = (TF1*)resp1->GetFunction("gaus");
   TF1* f2 = (TF1*)resp2->GetFunction("gaus");
   f1->SetLineColor(2);
   f2->SetLineColor(4);
   // f1->Delete();
   // f2->Delete();
   
   // TLegend* legends = new TLegend(0.1, 0.85, 0.35, 0.9,"","brNDC"); // the numbers determine the position of the box 
   // legends->SetFillColor(0); 
   // legends->SetHeader(stat_name);
   // //   legends->SetTextSize(40);
   // // // legends->AddEntry(resp1, "for old (10_0_3)","lep");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
   // // // legends->AddEntry(resp2, "for new (10_6_0_pre2)","lep");
   // legends->Draw();
   Canvas_1_n2->SaveAs(path);
   if (n[i]==494) continue;
    }
}
