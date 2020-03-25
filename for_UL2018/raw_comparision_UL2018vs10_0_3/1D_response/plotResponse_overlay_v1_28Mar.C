#include <iostream>
#include <iomanip>
#include <vector>
#include "TDirectory.h"

const int NEtaBins = 2;

TH1D* includeOverflow(TH1D *hin) {
  TH1D *hout = (TH1D*) hin->Clone();

  int nbins = hin->GetNbinsX();
  double overflow = hin->GetBinContent(nbins+1);
  double lastbin  = hin->GetBinContent(nbins);
  hout->SetBinContent(nbins, lastbin+overflow);
  return hout;		      
}


//void plotHists(TFile *file0, TFile *file1, char *histname, char* tag="", bool ifPrint=false){
void plotHists(TFile *file0, TFile *file1,TFile *file2, char const *histname, char const *tag="", bool ifPrint=false){

gStyle->SetStatW(0.22);
gStyle->SetStatH(0.12);
gStyle->SetOptFit(1111);
// gStyle->SetOptStat(1111111);
gStyle->SetTitleX(0.5);
gStyle->SetTitleY(0.955);
//gStyle->SetTitleW(0.8);

  int col1 = kBlack;
  int col2 = kRed;
  int col3 = kBlue;
  char hname[500], htit[500];
  bool rebin = true;
  
  file0->cd();

  sprintf(hname, "%s", histname);
  TH1D *h1 = (TH1D*)file0->FindObjectAny(hname);
  file1->cd();
  TH1D *h2 = (TH1D*)file1->FindObjectAny(hname);
  file2->cd();
  TH1D *h3 = (TH1D*)file2->FindObjectAny(hname);

  if( strcmp(histname,"h_HBHE_RecHitsEnergy") == 0)
    h1 = includeOverflow(h1);
  else if( strcmp(histname,"histcorhybrid2") == 0) {
    h1->SetTitle("Etrue 2-6 GeV"); h1->SetName("Etrue 2-6 GeV");
  }
  else if( strcmp(histname, "histcorhybrid10") == 0) {
    h1->SetTitle("Etrue 9-14 GeV"); h1->SetName("Etrue 9-14 GeV");
  }
  else if( strcmp(histname,"histcorhybrid52") == 0) {
    h1->SetTitle("Etrue 51-60 GeV"); h1->SetName("Etrue 51-60 GeV");
  }
  else if( strcmp(histname,"histcorhybrid44") == 0) {
    h1->SetTitle("Etrue 43-52 GeV"); 
  }
  else if( strcmp(histname,"histcorhybrid64") == 0) {
    h1->SetTitle("Etrue 43-72 GeV"); 
  }
  else if( strcmp(histname,"histcorhybrid80") == 0) {
    h1->SetTitle("Etrue 79-88 GeV"); 
  }

  
  else if( strcmp(histname,"histcorhybrid154") == 0) {
    h1->SetTitle("Etrue 153-174 GeV"); h1->SetName("Etrue 153-174 GeV");
  }
  else if( strcmp(histname,"histcorhybrid24") == 0){
    h1->SetTitle("Etrue 23-32 GeV"); h1->SetName("Etrue 23-32 GeV");
  }
  else if( strcmp(histname,"histcorhybrid60") == 0) {
    h1->SetTitle("Etrue 59-68 GeV"); h1->SetName("Etrue 59-68 GeV");
  }
  else if( strcmp(histname,"histcorhybrid124") == 0) {
    h1->SetTitle("Etrue 123-144 GeV"); h1->SetName("Etrue 123-144 GeV");
  }
  else if( strcmp(histname,"histcorhybrid274") == 0){
    h1->SetTitle("Etrue 273-294 GeV"); h1->SetName("Etrue 273-294 GeV");
  }
  
  else if( strcmp(histname,"histcorhybrid334") == 0){
    h1->SetTitle("Etrue 163-184 GeV"); h1->SetName("Etrue 163-184 GeV");
  }
  else if( strcmp(histname,"histcorhybrid334") == 0) {
    h1->SetTitle("Etrue 334-354 GeV"); h1->SetName("Etrue 333-354 GeV");
  }

  else if( strcmp(histname,"histcorhybrid444") == 0) {
    h1->SetTitle("Etrue 443-464 GeV"); h1->SetName("Etrue 443-464 GeV");
  }

  h1->SetXTitle("Response");
  //h1 = includeOverflow(h1);
  h1->SetLineColor(col1);
  h2->SetLineColor(col2);
  h3->SetLineColor(col3);
  //h2->SetLineColor(4);
  
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);

  h1->GetXaxis()->SetLabelSize(0.06);
  h2->GetXaxis()->SetLabelSize(0.06);
  h3->GetXaxis()->SetLabelSize(0.06);

  h1->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetLabelSize(0.05);
  h3->GetYaxis()->SetLabelSize(0.05);
   // h2->SetLineStyle(1);
   // h3->SetLineStyle(1);

  h1->Scale(1.0/h1->Integral());
  h2->Scale(1.0/h2->Integral());
  h3->Scale(1.0/h3->Integral());

  /* // ss
  double ymax=0.0;
  if (ymax < h1->GetBinContent(h1->GetMaximumBin())) ymax = h1->GetBinContent(h1->GetMaximumBin());  
  std::cout << "1 " << h1->GetBinContent(h1->GetMaximumBin()) << " " << ymax << std::endl;
  if (ymax < h2->GetBinContent(h2->GetMaximumBin())) ymax = h2->GetBinContent(h2->GetMaximumBin());
  std::cout << "2 " << h2->GetBinContent(h2->GetMaximumBin()) << " " << ymax << std::endl;
  if (ymax < h3->GetBinContent(h3->GetMaximumBin())) ymax = h3->GetBinContent(h3->GetMaximumBin());
  std::cout << "3 " << h3->GetBinContent(h3->GetMaximumBin()) << " " << ymax << std::endl;
  ymax = 1.1*ymax;
  h1->GetYaxis()->SetRangeUser(0, ymax);
  h2->GetYaxis()->SetRangeUser(0, ymax);
  h3->GetYaxis()->SetRangeUser(0, ymax);
  */
  // double ymax;
  // double y1 = h1->GetBinContent(h1->GetMaximumBin());
  // double y2 = h2->GetBinContent(h2->GetMaximumBin());
  // double y3 = h3->GetBinContent(h3->GetMaximumBin());
  // if (y1 > y2 ) y2=y1;
  // if (y2 > y3 ) ymax = y2;
  // else ymax = y3;
  // h1->GetYaxis()->SetRangeUser(0, ymax);                                                                                                      
  // h2->GetYaxis()->SetRangeUser(0, ymax);                                                                                                      
  // h3->GetYaxis()->SetRangeUser(0, ymax);    
  
  sprintf(htit, "c_%s", hname);
  TCanvas *c = new TCanvas(htit, htit, 400, 250);
  //TH1D* h1 = (TH1D*)file0->Draw();
  if (strcmp (histname, "histcorhybrid8")||strcmp (histname, "histcorhybrid12")||strcmp (histname, "histcorhybrid24")||strcmp (histname, "histcorhybrid44")) {
    h3->SetStats(1);
    h3->GetXaxis()->SetRangeUser(-1, 10);
    h3->Draw("hist");
    h2->Draw("hist sames");
    h2->GetXaxis()->SetRangeUser(-1, 10);
    h2->SetStats(1);
    h1->SetStats(1);
    h1->GetXaxis()->SetRangeUser(-1, 10);
    h1->Draw("hist sames");
  } 
  
  else { 
    h2->Draw("hist");
    h2->SetStats(1);
    h1->SetStats(1);
    h1->Draw("hist sames");
    h3->SetStats(1);
    h3->Draw("hist sames");
  }
  
  // h1->GetXaxis()->SetRange(0, 10);
  // h2->GetXaxis()->SetRange(0, 10);
  // h3->GetXaxis()->SetRange(0, 10);

  cout<<"histname:"<<histname<<endl;
  //if (tag == "EC_EH_beyond_tracker_v900" || tag == "EC_H_beyond_tracker_v900" || tag == "EC_EH_beyond_tracker_v900" ||
  //  || tag == "EC_EH_beyond_tracker_v910pre1" || tag == "EC_H_beyond_tracker_v910pre1") {
  if (true) {
  //if (histname == "histcorhybrid24") { 
    TF1 *f1 = (TF1*) h1->FindObject("gaus");
    delete f1;
    TF1 *f2 = (TF1*) h2->FindObject("gaus");
    delete f2;
    TF1 *f3 = (TF1*) h3->FindObject("gaus");
    delete f3;

  }

  
  if((strcmp(histname,"histcorhybrid8") == 0) || (strcmp(histname,"histcorhybrid12") == 0) || (strcmp(histname,"histcorhybrid24") == 0) || (strcmp(histname,"histcorhybrid44") == 0)) {
    if(rebin) {
      
      h1->Rebin(6);
      h2->Rebin(6);
      h3->Rebin(6);
    }
  }
  else {
    
    h1->GetXaxis()->SetRangeUser(-1.5,1.5);
    h2->GetXaxis()->SetRangeUser(-1.5,1.5);
    h3->GetXaxis()->SetRangeUser(-1.5,1.5);
  }
  
  gPad->Update();

  h1->SetName("10_0_2");
  h2->SetName("10_0_3");
  h3->SetName("10_0_3_aternate");
  TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats");
  //st1->AddText("yolo");
  st1->SetTextColor(h1->GetLineColor());
  st1->SetLineColor(h1->GetLineColor());
  //st1->SetText("hello");
  // st1->SetX1NDC(0.6);
  // st1->SetY1NDC(0.9);
  // st1->SetX2NDC(0.9);
  // st1->SetY2NDC(0.5);
  st1->SetX1NDC(0.7);  //x1 left start
  st1->SetY1NDC(0.9);  //
  st1->SetX2NDC(0.9);  //x2  right end
  st1->SetY2NDC(0.6);



  TPaveStats *st2 = (TPaveStats*)h2->FindObject("stats");
  st2->SetTextColor(h2->GetLineColor());
  st2->SetLineColor(h2->GetLineColor());
  //st2->SetText("hello");
  // st2->SetX1NDC(0.6);
  // st2->SetY1NDC(0.9);
  // st2->SetX2NDC(0.9);
  // st2->SetY2NDC(0.5);
  st2->SetX1NDC(0.7);  //x1 left start
  st2->SetY1NDC(0.6);  //
  st2->SetX2NDC(0.9);  //x2  right end
  st2->SetY2NDC(0.3);

  TPaveStats *st3 = (TPaveStats*)h3->FindObject("stats");
  st3->SetTextColor(h3->GetLineColor());
  st3->SetLineColor(h3->GetLineColor()); 
  st3->SetX1NDC(0.7);  //x1 left start                                                                                                        
  st3->SetY1NDC(0.3);  //                                                                                                                     
  st3->SetX2NDC(0.9);  //x2  right end                                                                                                        
  st3->SetY2NDC(0);
  // TLegend *l1 = new TLegend(0.7,0.6,0.9,0.9);
  // l1->AddEntry(h1,"E Corr","lp");
  // l1->AddEntry(h2,"Eta Corr","lp");
  // //l1->Draw();

                   
  gPad->Modified();
  //c->Modified();
  
  if(ifPrint){
    char cname[200];
    sprintf(cname, "%s_%s.gif", tag,htit);
    //std::cout << cname << std::endl;
    c->Print(cname);
  }

}


//void plotHists(char* fname1,char* fname2, char* tag="", bool ifPrint = false){
void plotHists(char const *fname1,char const *fname2, char const *fname3, char const *tag="", bool ifPrint = false){ 

  TFile *file0 = TFile::Open(fname1);
  TFile *file1 = TFile::Open(fname2);
  TFile *file2 = TFile::Open(fname3);
  // plotHists(file0, "histcorhybrid2",  tag, ifPrint);
  // plotHists(file0, "histcorhybrid10",  tag, ifPrint);
  // plotHists(file0, "histcorhybrid52",  tag, ifPrint);
  // plotHists(file0, "histcorhybrid154", tag, ifPrint);


  // plotHists(file0, "histcorhybrid24",  tag, ifPrint);
  // plotHists(file0, "histcorhybrid60",  tag, ifPrint);
  // plotHists(file0, "histcorhybrid124",  tag, ifPrint);
  // plotHists(file0, "histcorhybrid164", tag, ifPrint);


  // plotHists(file0, file1,  "histEta;1",  tag, ifPrint);
  // plotHists(file0, file1,  "histEta;5",  tag, ifPrint);
  // plotHists(file0, file1,  "histEta;9",  tag, ifPrint);
  // plotHists(file0, file1,  "histEta;13",  tag, ifPrint);


  // plotHists(file0, file1,  "histcorhybrid2",  tag, ifPrint);
  // plotHists(file0, file1,  "histcorhybrid4",  tag, ifPrint);
  // plotHists(file0, file1,  "histcorhybrid6",  tag, ifPrint);
  plotHists(file0, file1,file2 ,  "histcorhybrid8",  tag, ifPrint);
  plotHists(file0, file1,file2,  "histcorhybrid12",  tag, ifPrint);

  plotHists(file0, file1,file2, "histcorhybrid24",  tag, ifPrint);

  plotHists(file0, file1,file2, "histcorhybrid44",  tag, ifPrint);
  // plotHists(file0, file1,file2, "histcorhybrid60",  tag, ifPrint);
  plotHists(file0, file1,file2, "histcorhybrid124",  tag, ifPrint);
  plotHists(file0, file1,file2, "histcorhybrid164", tag, ifPrint);
    plotHists(file0, file1,file2, "histcorhybrid244", tag, ifPrint);
  
  // plotHists(file0, file1,  "histEta;5",  tag, ifPrint);
  // plotHists(file0, file1,  "histEta;9",  tag, ifPrint);
  // plotHists(file0, file1,  "histEta;13",  tag, ifPrint);

}


/*

// =========================================================================

void plotEtaVsPhiPlots(char* fname1, char *beam, bool ifPrint=false){

  TFile *file0 = TFile::Open(fname1);
  plotEtaVsPhiPlots(file0, "h_HBHE_EtavsPhiEnergy_1p5GeV", beam, ifPrint);
  plotEtaVsPhiPlots(file0, "h_HBHE_EtavsPhiEnergy_1p5GeV_Dep1", beam, ifPrint);
  plotEtaVsPhiPlots(file0, "h_HBHE_EtavsPhiEnergy_1p5GeV_Dep2", beam, ifPrint);
  plotEtaVsPhiPlots(file0, "h_HBHE_EtavsPhiEnergy_1p5GeV_Dep3", beam, ifPrint);
  plotEtaVsPhiPlots(file0, "h_HBHE_EtavsPhiEnergy_30GeV", beam, ifPrint);
  plotEtaVsPhiPlots(file0, "h_HBHE_EtavsPhiEnergy_30GeV_Dep1", beam, ifPrint);
  plotEtaVsPhiPlots(file0, "h_HBHE_EtavsPhiEnergy_30GeV_Dep2", beam, ifPrint);
  plotEtaVsPhiPlots(file0, "h_HBHE_EtavsPhiEnergy_30GeV_Dep3", beam, ifPrint);
}

void plotEtaVsPhiPlots(TFile *file0, char *histname, char *beam, bool ifPrint=false){

  
  Int_t palette[50] = {19,18,17,16,15,14,13,12,11,20,
		       21,22,23,24,25,26,27,28,29,30, 8,
		       31,32,33,34,35,36,37,38,39,40, 9,
		       41,42,43,44,45,47,48,49,46,50, 2,
		       7, 6, 5, 4, 3, 2,1};                                                                                                                                      
   for(int i=0; i<50; i++) {
     palette[i] = 51+i;
     //if(i==20 || i==23 || i==22 ) palette[i]=10;
     //if(i>20 && i<30) palette[i]=10;
     //if(i>23 && i<28) palette[i]=10;
     //if(i==24 || i==27) palette[i]=10;
     //std::cout << i <<" "<<palette[i]<<std::endl;
   }
   gStyle->SetPalette(50, palette);
  

  //gStyle->SetOptStat(11111);
  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleY(0.955);

  file0->cd();

  char hname[200], htit[200];
  sprintf(hname, "%s", histname);
  TH2F *h1 = (TH2F*)file0->FindObjectAny(hname);
  //h1->GetZaxis()->SetRangeUser(-600,600);

  sprintf(htit, "c_%s_%s", hname, beam);
  TCanvas *c = new TCanvas(htit, htit, 500, 500);
  h1->Draw("colz");

  if(ifPrint){
    char cname[200];
    sprintf(cname, "%s.eps", htit);
    //std::cout << cname << std::endl;
    c->Print(cname);
  }

}

void plotHBHEEtaRings(TFile *file0, bool ifPrint=false){

  //gStyle->SetOptTitle(0);
  gStyle->SetOptStat(11111);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleY(0.955);
  
  char hname[200];

  //const double xx[7] = {-0.1, 0.1, 0.3, 0.5, 0.7, 0.9, -0.1};
  const double xx[7] = {-0.1, 0.1, 0.3, 0.5, 0.7, 0.1, 0.3};
  const double yy[7] = {0.9, 0.9, 0.9, 0.9, 0.9, 0.7, 0.7};

  const int ietaCol[7] = {600, 632, 432, 616, 418, 1, 500};
  const int ietaTh[11] = {-28, -22, -16, -10, -5,-1, 5, 10, 16, 22, 28};
  const char* ietaThName[11] = {"m28", "m22", "m16", "m10", "m5","m1", "p5", "p10", "p16", "p22", "p28"};
  TCanvas *c[11];
  for(int ii=0; ii<10; ii++){
    sprintf(hname, "c_HBHEEtaRings_eta_%sTo%s", ietaThName[ii], ietaThName[ii+1]);
    c[ii] = new TCanvas(hname, hname, 500, 500);
  }
  
  int icnt=0;
  for (int ieta=0; ieta<60; ieta++) {
    if(ieta-30 <= 0) {
      sprintf(hname, "h_HBHEEtaRings_RecHitTime30GeV_etam%i", -1*(ieta-30) );
    } else if(ieta-30 > 0) {
      sprintf(hname, "h_HBHEEtaRings_RecHitTime30GeV_etap%i", ieta-30);
    }

    TH1F *h1 = (TH1F*)file0->FindObjectAny(hname);
    if(ieta-30 <= 0) {
      sprintf(hname, "i#eta = %i", (ieta-30) );
    } else if(ieta-30 > 0) {
      sprintf(hname, "i#eta = %i", ieta-30);
    }
    h1->SetName(hname);
    
    for(int ii=0; ii<10; ii++) {
      if( ieta-30 >= ietaTh[ii] && ieta-30 < ietaTh[ii+1]) {
	if(ieta-30 == ietaTh[ii]) icnt=0;
	//cout << ieta-30 << " "<< icnt << endl;

	h1->SetLineColor(ietaCol[icnt]);
	c[ii]->cd();
	if(ieta==ietaTh[ii]) 
	  h1->Draw();
	else
	  h1->Draw("sames");
	icnt++;

	gPad->Update();
	TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats");
	st1->SetTextColor(h1->GetLineColor());
	st1->SetLineColor(h1->GetLineColor());
	st1->SetX1NDC(xx[icnt]);
	st1->SetY1NDC(yy[icnt]);
	st1->SetX2NDC(xx[icnt]+0.2);
	st1->SetY2NDC(yy[icnt]-0.2);

	gPad->Modified();
	
      }

    }
    
  }

}





*************************************************************8

  faBarrel = new TF1("faBarrel","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);

  funtionBarrelEcalHcalB->FixParameter(0,-13.9219);
  funtionBarrelEcalHcalB->FixParameter(1,14.9124);
  funtionBarrelEcalHcalB->FixParameter(2,5.38578);
  funtionBarrelEcalHcalB->FixParameter(3,0.861981);
  funtionBarrelEcalHcalB->FixParameter(4,-0.00759275);
  funtionBarrelEcalHcalB->FixParameter(5,3.73563e-23);
  funtionBarrelEcalHcalB->FixParameter(6,-1.17946);
  funtionBarrelEcalHcalB->FixParameter(7,-13.3644);



  fbBarrel = new TF1("fbBarrel","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  functionBarrelEcalHcalC->FixParameter(0,3.35133);
  functionBarrelEcalHcalC->FixParameter(1,2.3293);
  functionBarrelEcalHcalC->FixParameter(2,-7.83787);
  functionBarrelEcalHcalC->FixParameter(3,1.2109e+06);
  functionBarrelEcalHcalC->FixParameter(4,4.742);
  functionBarrelEcalHcalC->FixParameter(5,0.427321);
  functionBarrelEcalHcalC->FixParameter(6,-0.383789);
  functionBarrelEcalHcalC->FixParameter(7,-0.545161);



  fcBarrel = new TF1("fcBarrel","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  functionBarrelHcalC->FixParameter(0,1.50693);
  functionBarrelHcalC->FixParameter(1,0.58792);
  functionBarrelHcalC->FixParameter(2,-4.41621);
  functionBarrelHcalC->FixParameter(3,4.24264);
  functionBarrelHcalC->FixParameter(4,0.603079);
  functionBarrelHcalC->FixParameter(5,0.0235189);
  functionBarrelHcalC->FixParameter(6,0.524996);
  functionBarrelHcalC->FixParameter(7,-1.58837);



  faEtaBarrelEH = new TF1("faEtaBarrelEH","[0]+[1]*exp(-x/[2])",1.,1000.);
  functionBarrelAlphaEcalHcal->FixParameter(0,0.0273791);
  functionBarrelAlphaEcalHcal->FixParameter(1,-0.0538802);
  functionBarrelAlphaEcalHcal->FixParameter(2,548.713);


  fbEtaBarrelEH = new TF1("fbEtaBarrelEH","[0]+[1]*exp(-x/[2])",1.,1000.);
  functionBarrelBetaEcalHcal->FixParameter(0,0.0437171);
  functionBarrelBetaEcalHcal->FixParameter(1,0.109082);
  functionBarrelBetaEcalHcal->FixParameter(2,234.472);

  faEtaBarrelH = new TF1("faEtaBarrelH","[0]+[1]*x",1.,1000.);
  functionBarrelAlphaHcal->FixParameter(0,0.00343543);
  functionBarrelAlphaHcal->FixParameter(1,-1.02499e-06);

  fbEtaBarrelH = new TF1("fbEtaBarrelH","[0]+[1]*exp(-x/[2])",1.,1000.);
  functionBarrelBetaHcal->FixParameter(0,-0.0285311);
  functionBarrelBetaHcal->FixParameter(1,0.0761233);
  functionBarrelBetaHcal->FixParameter(2,54.1946);

  faEndcap = new TF1("faEndcap","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  functionEndcapEcalHcalB->FixParameter(0,0.962468);
  functionEndcapEcalHcalB->FixParameter(1,11.9536);
  functionEndcapEcalHcalB->FixParameter(2,-27.7088);
  functionEndcapEcalHcalB->FixParameter(3,0.755474);
  functionEndcapEcalHcalB->FixParameter(4,0.0791012);
  functionEndcapEcalHcalB->FixParameter(5,2.6901e-11);
  functionEndcapEcalHcalB->FixParameter(6,0.158734);
  functionEndcapEcalHcalB->FixParameter(7,-6.92163);

  fbEndcap = new TF1("fbEndcap","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  functionEndcapEcalHcalC->FixParameter(0,-0.43671);
  functionEndcapEcalHcalC->FixParameter(1,2.90096);
  functionEndcapEcalHcalC->FixParameter(2,-5.10099);
  functionEndcapEcalHcalC->FixParameter(3,1.20771);
  functionEndcapEcalHcalC->FixParameter(4,-1.30656);
  functionEndcapEcalHcalC->FixParameter(5,0.0189607);
  functionEndcapEcalHcalC->FixParameter(6,0.270027);
  functionEndcapEcalHcalC->FixParameter(7,-2.30372);

  fcEndcap = new TF1("fcEndcap","[0]+((([1]+([2]/sqrt(x)))*exp(-(x^[6]/[3])))-([4]*exp(-(x^[7]/[5]))))",1.,1000.);
  functionEndcapHcalC->FixParameter(0,1.13795);
  functionEndcapHcalC->FixParameter(1,1.21698);
  functionEndcapHcalC->FixParameter(2,-3.81192);
  functionEndcapHcalC->FixParameter(3,60.0406);
  functionEndcapHcalC->FixParameter(4,0.673456);
  functionEndcapHcalC->FixParameter(5,0.217077);
  functionEndcapHcalC->FixParameter(6,1.95596);
  functionEndcapHcalC->FixParameter(7,-0.252215);

  faEtaEndcapEH = new TF1("faEtaEndcapEH","[0]+[1]*exp(-x/[2])",1.,1000.);
  functionEndcapAlphaEcalHcal->FixParameter(0,-0.0426257);
  functionEndcapAlphaEcalHcal->FixParameter(1,0.0406977);
  functionEndcapAlphaEcalHcal->FixParameter(2,-1143.46);

  fbEtaEndcapEH = new TF1("fbEtaEndcapEH","[0]+[1]*exp(-x/[2])",1.,1000.);
  functionEndcapBetaEcalHcal->FixParameter(0,0.0312776);
  functionEndcapBetaEcalHcal->FixParameter(1,-1.31529);
  functionEndcapBetaEcalHcal->FixParameter(2,11.5701);

  faEtaEndcapH = new TF1("faEtaEndcapH","[0]+[1]*x",1.,1000.);
  functionEndcapAlphaHcal->FixParameter(0,-0.0825626);
  functionEndcapAlphaHcal->FixParameter(1,8.68807e-05);

  fbEtaEndcapH = new TF1("fbEtaEndcapH","[0]+[1]*exp(-x/[2])+[3]*[3]*exp(-x*x/([4]*[4]))",1.,1000.);
  functionEndcapBetaHcal->FixParameter(0,0.0719133);
  functionEndcapBetaHcal->FixParameter(1,-0.27956);
  functionEndcapBetaHcal->FixParameter(2,5.32799);
  functionEndcapBetaHcal->FixParameter(3,0.272516);
  functionEndcapBetaHcal->FixParameter(4,161.719);


*/
