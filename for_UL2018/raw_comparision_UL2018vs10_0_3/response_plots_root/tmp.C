void tmp()
{
//=========Macro generated from canvas: Canvas_1_n2/Canvas_1_n2
//=========  (Tue Apr  2 22:04:30 2019) by ROOT version 6.14/06
   TCanvas *Canvas_1_n2 = new TCanvas("Canvas_1_n2", "Canvas_1_n2",65,52,700,500);
   Canvas_1_n2->Range(-0.375,-0.5,3.375,0.5);
   Canvas_1_n2->SetFillColor(0);
   Canvas_1_n2->SetBorderMode(0);
   Canvas_1_n2->SetBorderSize(2);
   Canvas_1_n2->SetFrameBorderMode(0);
   Canvas_1_n2->SetFrameBorderMode(0);
   
   TH2F *respHisto__1 = new TH2F("respHisto__1","Response",30,0,3,100,-1,1);
   respHisto__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   respHisto__1->SetLineColor(ci);
   respHisto__1->GetXaxis()->SetTitle("|#eta|");
   respHisto__1->GetXaxis()->SetRange(1,30);
   respHisto__1->GetXaxis()->SetLabelFont(42);
   respHisto__1->GetXaxis()->SetLabelSize(0.035);
   respHisto__1->GetXaxis()->SetTitleSize(0.035);
   respHisto__1->GetXaxis()->SetTitleFont(42);
   respHisto__1->GetYaxis()->SetTitle("(E_{cor}-E_{true})/E_{true}");
   respHisto__1->GetYaxis()->SetRange(31,70);
   respHisto__1->GetYaxis()->SetLabelFont(42);
   respHisto__1->GetYaxis()->SetLabelSize(0.035);
   respHisto__1->GetYaxis()->SetTitleSize(0.035);
   respHisto__1->GetYaxis()->SetTitleOffset(0);
   respHisto__1->GetYaxis()->SetTitleFont(42);
   respHisto__1->GetZaxis()->SetLabelFont(42);
   respHisto__1->GetZaxis()->SetLabelSize(0.035);
   respHisto__1->GetZaxis()->SetTitleSize(0.035);
   respHisto__1->GetZaxis()->SetTitleFont(42);
   respHisto__1->Draw("");
   
   TPaveText *pt = new TPaveText(0.4076218,0.9338535,0.5923782,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Response");
   pt->Draw();
   Canvas_1_n2->Modified();
   Canvas_1_n2->cd();
   Canvas_1_n2->SetSelected(Canvas_1_n2);
}
