#include <stdio.h>
// #include<conio.h>
void generate_1Dplots_overlay_v2(string file_name,string file_name1)
{
  char* hname = new char[200];
  char* hist_name = new char[200];
  char* hist_name1 = new char[200];
  char* hname1 = new char[2000];
  char* png_path = new char[2000];

  sprintf(hname,"%s",file_name.c_str());
  sprintf(hname1,"%s",file_name1.c_str());

     
  TFile * inputfile = new TFile(hname,"READ");
  TFile * inputfile1 = new TFile(hname1,"READ");
  
  sprintf(hist_name,"Layer1");
  sprintf(hist_name1,"Layer2");
  sprintf(png_path,"%s_vs_%s.png",hist_name,hist_name1);
 
  TCanvas *Canvas_1_n2 = new TCanvas(hist_name, hist_name,65,52,525,527);
  Canvas_1_n2->Range(-60.25,-0.625,562.25,0.625);
  Canvas_1_n2->SetFillColor(0);
  Canvas_1_n2->SetBorderMode(0);
  Canvas_1_n2->SetBorderSize(2);
  Canvas_1_n2->SetGridx();
  Canvas_1_n2->SetGridy();
  Canvas_1_n2->SetFrameBorderMode(0);
  Canvas_1_n2->SetFrameBorderMode(0);
  
  TH1D *resp = (TH1D*)inputfile->Get(hist_name);
  TH1D *resp1 = (TH1D*)inputfile1->Get(hist_name1);


  TPaveStats *ptstats2 = new TPaveStats(0.65,0.47,0.9,0.67,"brNDC");
  ptstats2->SetBorderSize(1);
  ptstats2->SetFillColor(0);
  ptstats2->SetLineColor(kRed);
  ptstats2->SetTextAlign(12);
  ptstats2->SetTextColor(kRed);
  ptstats2->SetTextFont(42);
  /* TText *ptstats_LaTex = ptstats2->AddText(old_rel); */
  /* ptstats_LaTex->SetTextSize(0.0315317); */
  ptstats2->SetOptStat(1111);
  ptstats2->SetOptFit(10001);
  resp1->GetListOfFunctions()->Add(ptstats2);
  ptstats2->SetParent(resp1);
  resp1->Draw();
  resp->Draw("sames");

  Canvas_1_n2->Modified();
  Canvas_1_n2->cd();
  Canvas_1_n2->SetSelected(Canvas_1_n2);


  Canvas_1_n2->SaveAs(png_path);
    

  
}
