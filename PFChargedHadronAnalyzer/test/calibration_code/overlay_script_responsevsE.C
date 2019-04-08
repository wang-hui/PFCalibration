#include <stdio.h>
// #include<conio.h>
void overlay_script_responsevsE()
{
//=========Macro generated from canvas: Canvas_1_n2/Canvas_1_n2
//=========  (Mon Apr  1 19:35:54 2019) by ROOT version 6.14/06
  char* hname = new char[200];
  char* hname1 = new char[200];
  char* path = new char[200];
  char* plot_name = new char[200];
  char* old_rel = new char[200];
  char* new_rel = new char[200];

  char old_release[50] = "10_0_3_responseEta";
  char new_release[50] = "10_6_0_pre2_responseEta";
  char plot[50] = "Raw_EH";

  sprintf(hname,"%s_%s.root",plot,old_release);
  sprintf(hname1,"%s_%s.root",plot,new_release);
  sprintf(plot_name,"%s",plot);
  sprintf(path,"../overlay_response_plots/%s.gif",plot);
  sprintf(old_rel,"%s",old_release);
  sprintf(new_rel,"%s",new_release);


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
   TGraph* graph1 = (TGraph*) inputfile1 -> Get("Graph");
   TGraph* graph2 = (TGraph*) inputfile2 -> Get("Graph");
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
   graph1->SetMarkerSize(1.0);

   graph2->SetName("Graph2");
   graph2->SetTitle("Graph");
   graph2->SetFillStyle(1000);
   graph2->SetMarkerColor(4);
   graph2->SetMarkerStyle(22);
   graph2->SetMarkerSize(1.0);
   
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
   
   graph1->Draw("p");

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

   // TLine *line = new TLine(1.2273,0,498.0851,0);
   TLine *line = new TLine(0,0,3,0);
   line->SetLineColor(2);
   line->SetLineWidth(2);
   line->Draw();
   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);


   TLegend* legends = new TLegend(0.1, 0.7, 0.44, 0.9,"","brNDC"); // the numbers determine the position of the box 
   legends->SetFillColor(0); 
   legends->SetHeader(plot_name); 
   legends->AddEntry(graph1, "for old (10_0_3)","lep");//(name of hist,what you want it called in legend, l=line, p=polymarker, f=boxy thing ) 
   legends->AddEntry(graph2, "for new (10_6_0_pre2)","lep");
   legends->Draw();
   Canvas_1_n2->SaveAs(path);

}
