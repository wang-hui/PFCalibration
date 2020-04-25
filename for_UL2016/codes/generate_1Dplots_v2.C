#include <stdio.h>
// #include<conio.h>
void generate_1Dplots_v2(string plot_name)
{
  char* hname = new char[200];
  char* hist_name = new char[200];
  char* stat_name = new char[200];
  char* hname1 = new char[2000];
  char* path = new char[2000];
  char* full_path = new char[2000];
  char* path1 = new char[2000];
  char* full_path1 = new char[2000];
  char* path2 = new char[2000];
  char* full_path2 = new char[2000];
  //  char* plot_name = new char[200];
  char* old_rel = new char[200];
  char* new_rel = new char[200];
  char* old_release = new char[50];
  char* new_release = new char[50];
  char* hist = new char[200];
  // sprintf(old_release, "corrEta");
  // sprintf(new_release, "corrEta");
  //  sprintf(plot_name, "EH_barrel");
  //  sprintf(plot_name, "EH_ec_out");
  sprintf(hist, "histcorhybrid");

  //comparision for old calib and new calib
  sprintf(hname,"%s",plot_name.c_str());
  sprintf(path2,"%s",plot_name.c_str());

      
  int nlayer= 50;

  TFile * inputfile1 = new TFile(hname,"READ");
  
  for(int i=1; i<=nlayer;i++)
    {

      sprintf(hist_name,"Layer%d",i);
      sprintf(full_path2,"Layer_0%d.png",i);
 
   TCanvas *Canvas_1_n2 = new TCanvas(hist_name, hist_name,65,52,525,527);
   Canvas_1_n2->Range(-60.25,-0.625,562.25,0.625);
   Canvas_1_n2->SetFillColor(0);
   Canvas_1_n2->SetBorderMode(0);
   Canvas_1_n2->SetBorderSize(2);
   Canvas_1_n2->SetGridx();
   Canvas_1_n2->SetGridy();
   Canvas_1_n2->SetFrameBorderMode(0);
   Canvas_1_n2->SetFrameBorderMode(0);
   TH1D *resp1 = (TH1D*)inputfile1->Get(hist_name);


   TPaveStats *ptstats2 = new TPaveStats(0.65,0.47,0.9,0.67,"brNDC");
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
   //  ptstats2->Draw();
   resp1->GetListOfFunctions()->Add(ptstats2);
   ptstats2->SetParent(resp1);
   resp1->Draw();


   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);


   /* Canvas_1_n2->SaveAs(full_path); */
   /* Canvas_1_n2->SaveAs(full_path1); */
   Canvas_1_n2->SaveAs(full_path2);
    }

  
}
