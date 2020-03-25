#include <stdio.h>
// #include<conio.h>
void overlay_script_responsevsE()
{
//=========Macro generated from canvas: Canvas_1_n2/Canvas_1_n2
//=========  (Mon Apr  1 19:35:54 2019) by ROOT version 6.14/06
  char* hname = new char[200];
  char* hist_name = new char[200];
  char* stat_name = new char[200];
  char* hname1 = new char[200];
  char* path = new char[200];
  char* plot_name = new char[200];
  char* legendname = new char[200];
  char* old_rel = new char[200];
  char* new_rel = new char[200];
  char* old_release = new char[50];
  char* new_release = new char[50];
  char* hist = new char[200];
  char* path_root1 = new char[200];
  char* path_root2 = new char[200];
  char* rootfile1 = new char[200];
  char* rootfile2 = new char[200];


  sprintf(legendname,"H barrel");
  sprintf(plot_name, "H_barrel_10_6_0_pre2_pt");
  sprintf(path_root1,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2017/trial_oldparam/response_plots_root");
  sprintf(path_root2,"/home/work/bhumika/work/calibration/Ultra_legacy/for_UL2017/trial1/response_plots_root");
  sprintf(rootfile1,"correta");
  sprintf(rootfile2,"correta");

  sprintf(hname,"%s/%s_%s.root",path_root1,rootfile1,plot_name);
  sprintf(hname1,"%s/%s_%s.root",path_root2,rootfile2,plot_name);
  // sprintf(plot_name,"%s",plot);
  //  sprintf(hist_name,"%s",hist);
  sprintf(path,"../overlay_response_plots/withraw/%s_%s.gif",rootfile1,plot_name);
  // sprintf(old_rel,"%s",old_release);
  // sprintf(new_rel,"%s",new_release);


   TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",65,52,525,527);
   Canvas_1_n2->Range(-60.25,-0.625,562.25,0.625);
   Canvas_1_n2->SetFillColor(0);
   Canvas_1_n2->SetBorderMode(0);
   Canvas_1_n2->SetBorderSize(2);
   Canvas_1_n2->SetGridx();
   Canvas_1_n2->SetGridy();
   Canvas_1_n2->SetFrameBorderMode(0);
   Canvas_1_n2->SetFrameBorderMode(0);
   TFile * inputfile1 = new TFile(hname,"READ");
   TFile * inputfile2 = new TFile(hname1,"READ");
   TGraph* graph1 = (TGraph*) inputfile1 -> Get("Graph;1");
   TGraph* graph2 = (TGraph*) inputfile2 -> Get("Graph;1");
   TH2F *respHisto__1 = (TH2F*) inputfile1 -> Get("respHisto");

   // respHisto__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   respHisto__1->SetLineColor(ci);
   respHisto__1->Draw("");
   
   TPaveText *pt = new TPaveText(0.4076218,0.9338535,0.5923782,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Response");
   pt->Draw();

   
   graph1->SetName("Graph1");
   graph1->SetTitle("Graph");
   graph1->SetFillStyle(1000);
   graph1->SetMarkerColor(2);
   graph1->SetMarkerStyle(22);
   graph1->SetMarkerSize(1.1);

   graph2->SetName("Graph2");
   graph2->SetTitle("Graph");
   graph2->SetFillStyle(1000);
   graph2->SetMarkerColor(4);
   graph2->SetMarkerStyle(22);
   graph2->SetMarkerSize(1.1);
   
   TH1F *Graph_Graph01 = new TH1F("Graph_Graph01","Graph",100,0,554.4);
   Graph_Graph01->SetMinimum(-540.4823);
   Graph_Graph01->SetMaximum(5940.616);
   Graph_Graph01->SetDirectory(0);
   Graph_Graph01->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph01->SetLineColor(ci);
   Graph_Graph01->GetXaxis()->SetLabelFont(42);
   Graph_Graph01->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetXaxis()->SetTitleFont(42);
   Graph_Graph01->GetYaxis()->SetLabelFont(42);
   Graph_Graph01->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetYaxis()->SetTitleOffset(0);
   Graph_Graph01->GetYaxis()->SetTitleFont(42);
   Graph_Graph01->GetZaxis()->SetLabelFont(42);
   Graph_Graph01->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetZaxis()->SetTitleFont(42);
   graph1->SetHistogram(Graph_Graph01);
   

   TH1F *Graph_Graph12 = new TH1F("Graph_Graph12","Graph",100,0,554);
   Graph_Graph12->SetMinimum(0.0);
   Graph_Graph12->SetMaximum(0.0);
   Graph_Graph12->SetDirectory(0);
   Graph_Graph12->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph12->SetLineColor(ci);
   Graph_Graph12->GetXaxis()->SetLabelFont(42);
   Graph_Graph12->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleFont(42);
   Graph_Graph12->GetYaxis()->SetLabelFont(42);
   Graph_Graph12->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleOffset(0);
   Graph_Graph12->GetYaxis()->SetTitleFont(42);
   Graph_Graph12->GetZaxis()->SetLabelFont(42);
   Graph_Graph12->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleFont(42);
   graph2->SetHistogram(Graph_Graph12);
   
   graph2->Draw("p");
   graph1->Draw("p");

   // TLine *line = new TLine(1.2273,0,498.0851,0);
   TLine *line = new TLine(0,0,100,0);
   line->SetLineColor(2);
   line->SetLineWidth(2);
   line->Draw();
   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);


   TLegend* legends = new TLegend(0.44, 0.7, 0.9, 0.9,"","brNDC"); // the numbers determine the position of the box 
   legends->SetFillColor(0); 
   legends->SetHeader(legendname); 
   // legends->AddEntry(graph1,"Raw (no correction)","P");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
   // legends->AddEntry(graph2,"final correction","P");
   legends->AddEntry(graph1,"old calib (10_0_3)","P");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
   legends->AddEntry(graph2,"new calib","P");
   legends->SetTextSize(0.04);
   //   legends->SetMarkerStyle(1);
   legends->Draw();
   Canvas_1_n2->SaveAs(path);

}
