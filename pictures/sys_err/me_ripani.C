TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;


TH1D *m_pip_p_bin_corr,*m_pip_pim_bin_corr,*m_pim_p_bin_corr;
TH1D *theta_p_bin_corr,*theta_pim_bin_corr,*theta_pip_bin_corr;
TH1D *alpha_p_bin_corr,*alpha_pim_bin_corr,*alpha_pip_bin_corr;

TH1D *m_pip_p_model,*m_pip_pim_model,*m_pim_p_model;
TH1D *theta_p_model,*theta_pim_model,*theta_pip_model;
TH1D *alpha_p_model,*alpha_pim_model,*alpha_pip_model;


TH1D *h_w[12];
TH1D *h_w_model[12];
TH1D *h_w_err[12];
Float_t full_err[22][12];

TGraphAsymmErrors* gae;

TCanvas *c = new TCanvas("c","c",800,600);


Float_t Q2_bin,W_bin[30];

Double_t w_err;



TSpline5 *spline;
ostringstream qqq;


 

Double_t f_fit(Double_t *x, Double_t *par) {

return spline->Eval(x[0]);

}


TH1D *h_bin_corr_theta(TH1D *h) {

TH1D *h_out;

Double_t *x = new double[(h->GetNbinsX())-1];
Double_t *y = new double[(h->GetNbinsX())-1];

for (Int_t aa = 1; aa <=((h->GetNbinsX())-1); aa++) {
x[aa-1] = (h->GetBinCenter(aa)+ h->GetBinCenter(aa+1))/2.;
y[aa-1] = (h->GetBinContent(aa) + h->GetBinContent(aa+1))/2.;
};


TGraph *gr = new TGraphErrors(((h->GetNbinsX())-1),x,y);
spline = new TSpline5("spline",gr);
spline->SetLineWidth(2);

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(1);
h_out->SetMarkerColor(1);


qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out->SetName(qqq.str().c_str());
qqq.str("");


TF1 *f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

TH1D *tmp = new TH1D("tmp","tmp",100,h->GetBinLowEdge(0),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1));
for (Int_t subbin=1; subbin<101;subbin++) {

tmp->SetBinContent(subbin,f->Eval(tmp->GetBinCenter(subbin)));
};
tmp->SetLineWidth(2);
//tmp->Draw("same c");

for (Int_t bin=1; bin<=(h_out->GetNbinsX());bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

/*
for (Int_t subbin=0; subbin<100;subbin++) {

position = left + width*subbin/99.;
avrg = avrg + spline->Eval(position);

};
avrg = avrg/100.;
*/

avrg = (f->Integral(left,left+width))/width;



if ((h->GetBinContent(bin)) == 0.) avrg = 0;

//cout << m_pip_p->GetBinContent(bin) << "  " << avrg << endl;

if ((h->GetBinContent(bin)) != 0.) {
h_out->SetBinContent(bin,(h->GetBinContent(bin))*(f->Eval(h->GetBinCenter(bin)))/avrg);
h_out->SetBinError(bin,(h->GetBinError(bin))/(h->GetBinContent(bin))*(h_out->GetBinContent(bin)));
} else {
h_out->SetBinContent(bin,0.);
h_out->SetBinError(bin,0.000001);
};

};




return h_out;


}; 





TH1D *h_alpha_sym(TH1D *h) {

TH1D *h_out;




Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName();
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(1);
h_out->SetMarkerColor(1);


qqq.str("");
qqq << h->GetName();
h_out->SetName(qqq.str().c_str());
qqq.str("");




for (Int_t bin=1; bin<=Int_t((h_out->GetNbinsX())/2.);bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

h_out->SetBinContent(bin,((h->GetBinContent(bin))+(h->GetBinContent((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinContent(((h->GetNbinsX())-bin+1),((h->GetBinContent(bin))+(h->GetBinContent((h->GetNbinsX())-bin+1)))/2.);


h_out->SetBinError(bin,sqrt((h->GetBinError(bin))*(h->GetBinError(bin))+(h->GetBinError((h->GetNbinsX())-bin+1))*(h->GetBinError((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinError(((h->GetNbinsX())-bin+1),sqrt((h->GetBinError(bin))*(h->GetBinError(bin))+(h->GetBinError((h->GetNbinsX())-bin+1))*(h->GetBinError((h->GetNbinsX())-bin+1)))/2.);


};

return h_out;


}; 



TH1D *h_bin_corr(TH1D *h) {

TH1D *h_out;


TGraph *gr = new TGraph(h);
spline = new TSpline5("spline",gr);
spline->SetLineWidth(2);

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(1);
h_out->SetMarkerColor(1);


qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out->SetName(qqq.str().c_str());
qqq.str("");




for (Int_t bin=1; bin<=(h_out->GetNbinsX());bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

for (Int_t subbin=0; subbin<100;subbin++) {

position = left + width*subbin/99.;
avrg = avrg + spline->Eval(position);

};
avrg = avrg/100.;

if ((h->GetBinContent(bin)) == 0.) avrg = 0;

//cout << m_pip_p->GetBinContent(bin) << "  " << avrg << endl;

if ((h->GetBinContent(bin)) != 0.) {
h_out->SetBinContent(bin,(h->GetBinContent(bin))*(spline->Eval(h->GetBinCenter(bin)))/avrg);
h_out->SetBinError(bin,(h->GetBinError(bin))/(h->GetBinContent(bin))*(h_out->GetBinContent(bin)));
} else {
h_out->SetBinContent(bin,0.);
h_out->SetBinError(bin,0.000001);
};

};

return h_out;


}; 


void draw_1d_hist (Int_t canvas, TH1D *h, string title, string name, string ytitle, string xtitle, Int_t color, string draw_options, string distr_flag,Int_t i) {

c->cd(canvas);
c->cd(canvas)->SetBottomMargin(0.2);
c->cd(canvas)->SetTopMargin(0.1);
c->cd(canvas)->SetLeftMargin(0.23);
c->cd(canvas)->SetRightMargin(0.01);
gPad->SetFillStyle(0);

h->SetMarkerStyle(20);
h->SetMarkerColor(color);
h->SetLineColor(color);
h->SetOption("pX0");
h->SetTitle(title.c_str());
h->SetTitleSize(0.1);


h->SetName(name.c_str());


 h->GetYaxis()->SetTitle(ytitle.c_str());
 h ->GetXaxis()->SetTitle(xtitle.c_str());
 h ->GetXaxis()->SetTitleSize(0.08);
 h->GetYaxis()->SetTitleSize(0.08);
 h->GetYaxis()->SetTitleOffset(1.3);
 h->GetXaxis()->SetLabelSize(0.08);
 h->GetXaxis()->SetNdivisions(6);
 h->GetYaxis()->SetLabelSize(0.07);
 h->GetYaxis()->SetNdivisions(5);


h->SetAxisRange(0, h-> GetMaximum()+(h-> GetMaximum())/2., "Y");
if ((name == "model_thetaP_")&&(W_bin[i]>1.65))h->SetAxisRange(0, h-> GetMaximum()+(h-> GetMaximum()), "Y");


};



Int_t get_max_w (Float_t Q2_bin) {
Int_t get_max_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.5))get_max_w = 21;
if ((Q2_bin> 0.5)&&(Q2_bin< 0.6))get_max_w = 20;
if ((Q2_bin> 0.6)&&(Q2_bin< 0.65))get_max_w = 19;
if ((Q2_bin> 0.65)&&(Q2_bin< 0.7))get_max_w = 18;
if ((Q2_bin> 0.7)&&(Q2_bin< 0.75))get_max_w = 17;
if ((Q2_bin> 0.75)&&(Q2_bin< 0.8))get_max_w = 16;
if ((Q2_bin> 0.8)&&(Q2_bin< 0.85))get_max_w = 14;
if ((Q2_bin> 0.85)&&(Q2_bin< 0.9))get_max_w = 13;
if ((Q2_bin> 0.9)&&(Q2_bin< 0.95))get_max_w = 12;
if ((Q2_bin> 0.95)&&(Q2_bin< 1.0))get_max_w = 10;
if ((Q2_bin> 1.0)&&(Q2_bin< 1.05))get_max_w = 9;
if ((Q2_bin> 1.05)&&(Q2_bin< 1.1))get_max_w = 7;
if ((Q2_bin> 1.1)&&(Q2_bin< 1.15))get_max_w = 6;
if ((Q2_bin> 1.15)&&(Q2_bin< 1.2))get_max_w = 4;
if ((Q2_bin> 1.2)&&(Q2_bin< 1.25))get_max_w = 2;
if ((Q2_bin> 1.25)&&(Q2_bin< 1.3))get_max_w = 1;
return get_max_w;
};

Int_t get_min_w (Float_t Q2_bin) {
Int_t get_min_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 12;
return get_min_w;
};



void read_data_rec(TFile *file_eff,Int_t i) {

file_eff->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p);
qqq.str("");




qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_pim_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p);
qqq.str("");





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_p);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pim);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pip);
qqq.str("");






qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_p);
qqq.str("");





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim);
qqq.str("");




qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip);
qqq.str("");



};

void read_data_eff(TFile *file_eff,Int_t i) {

file_eff->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p_bin_corr);
qqq.str("");




qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_pim_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim_bin_corr);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p_bin_corr);
qqq.str("");





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_p_bin_corr);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pim_bin_corr);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pip_bin_corr);
qqq.str("");






qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_p_bin_corr);
qqq.str("");





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim_bin_corr);
qqq.str("");




qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip_bin_corr);
qqq.str("");



};


void draw_1d_canvas( Int_t i, Int_t qq2 ) {

c = new TCanvas("c","c",700,700);

//c->Divide(3,3);



qqq << "Q^{2} = " << Q2_bin << " GeV^{2}" <<"; W = " << W_bin[i] <<" GeV" ;

Float_t max;

max = m_pip_p->GetMaximum();

if ((m_pip_pim->GetMaximum()) > max) max = m_pip_pim->GetMaximum();

if ((m_pim_p->GetMaximum()) > max) max = m_pim_p->GetMaximum();


for (Int_t kk=1; kk<4; kk++) { 
c->cd(kk)->SetBottomMargin(0.2);
c->cd(kk)->SetTopMargin(0.15);
c->cd(kk)->SetLeftMargin(0.23);
c->cd(kk)->SetRightMargin(0.01);
c->cd(kk)->SetFillStyle(4000);
};

for (Int_t kk=4; kk<10; kk++) { 
c->cd(kk)->SetBottomMargin(0.2);
c->cd(kk)->SetTopMargin(0.005);
c->cd(kk)->SetLeftMargin(0.23);
c->cd(kk)->SetRightMargin(0.01);
c->cd(kk)->SetFillStyle(4000);
};
gPad->SetFillStyle(4000);

c->cd(1);
m_pip_p->SetLineColor(1);
m_pip_p->SetMarkerColor(1);
m_pip_p->SetMinimum(0);
m_pip_p->SetMaximum(1.2*max);
m_pip_p->Draw("e1");
m_pip_p->Draw("e0 same");


c->cd(2);
m_pip_pim->SetLineColor(1);
m_pip_pim->SetMarkerColor(1);
m_pip_pim->SetMinimum(0);
m_pip_pim->SetMaximum(1.2*max);
m_pip_pim->Draw("e1");
m_pip_pim->Draw("e0 same");


c->cd(3);
m_pim_p->SetLineColor(1);
m_pim_p->SetMarkerColor(1);
m_pim_p->SetMinimum(0);
m_pim_p->SetMaximum(1.2*max);
m_pim_p->Draw("e1");
m_pim_p->Draw("e0 same");


max = theta_p->GetMaximum();

if ((theta_pip->GetMaximum()) > max) max = theta_pip->GetMaximum();

if ((theta_pim->GetMaximum()) > max) max = theta_pim->GetMaximum();

c->cd(4);
theta_p->SetLineColor(1);
theta_p->SetMarkerColor(1);
theta_p->SetMinimum(0);
theta_p->SetMaximum(1.2*max);
theta_p->GetXaxis()->SetTitle("#theta_{p'} (deg)");
theta_p->Draw("e1");


c->cd(5);
theta_pim->SetLineColor(1);
theta_pim->SetMarkerColor(1);
theta_pim->SetMinimum(0);
theta_pim->SetMaximum(1.2*max);
theta_pim->GetXaxis()->SetTitle("#theta_{#pi^{-}} (deg)");
theta_pim->Draw("e1");

c->cd(6);
theta_pip->SetLineColor(1);
theta_pip->SetMarkerColor(1);
theta_pip->SetMinimum(0);
theta_pip->SetMaximum(1.2*max);
theta_pip->GetXaxis()->SetTitle("#theta_{#pi^{+}} (deg)");
theta_pip->Draw("e1");

max = alpha_p->GetMaximum();

if ((alpha_pip->GetMaximum()) > max) max = alpha_pip->GetMaximum();

if ((alpha_pim->GetMaximum()) > max) max = alpha_pim->GetMaximum();


c->cd(7);
alpha_p->SetLineColor(1);
alpha_p->SetMarkerColor(1);
alpha_p->SetMinimum(0);
alpha_p->SetMaximum(1.2*max);
alpha_p->Draw("e1");


c->cd(8);
alpha_pim->SetLineColor(1);
alpha_pim->SetMarkerColor(1);
alpha_pim->SetMinimum(0);
alpha_pim->SetMaximum(1.2*max);
alpha_pim->GetXaxis()->SetTitle("#alpha_{#pi^{-}} (deg)");
alpha_pim->Draw("e1");


c->cd(9);
alpha_pip->SetLineColor(1);
alpha_pip->SetMarkerColor(1);
alpha_pip->SetMinimum(0);
alpha_pip->SetMaximum(1.2*max);
alpha_pip->Draw("e1");


c->cd();
TLatex latex;
latex.SetTextAlign(12);  //centered
latex.DrawLatex(.2,.9745,qqq.str().c_str());
qqq.str("");
c->Update();

qqq << "Q_" << Q2_bin << "_W_" << W_bin[i] <<".pdf";
c->Print(qqq.str().c_str());
qqq.str("");
//c->Delete();

};

void read_EG(TFile *file_eff,Int_t i) {

file_eff->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_pipP_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p_model);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_pippim_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim_model);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_pimP_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p_model);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_thetaP_";
gDirectory->GetObject(qqq.str().c_str(),theta_p_model);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_thetapim_";
gDirectory->GetObject(qqq.str().c_str(),theta_pim_model);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_thetapip_";
gDirectory->GetObject(qqq.str().c_str(),theta_pip_model);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_alpha_proton_";
gDirectory->GetObject(qqq.str().c_str(),alpha_p_model);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_alpha_pim_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim_model);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/model_alpha_pip_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip_model);
qqq.str("");






};




void me_ripani() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
gStyle->SetTitleSize(0.07,"t");
gStyle->SetTitleY(1.01);
gStyle->SetOptStat(0);
gStyle->SetErrorX(0);
gErrorIgnoreLevel = kError;
gStyle->SetStatY(0.88); 


 gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   Float_t small = 1e-5;
//   c1.Divide(2,2,small,small);

//c->Divide(4,3,small,small);

c->cd();
TLegend *leg = new TLegend(0.2,0.81,0.99,0.95);
leg->SetNColumns(4);

//TLegend *leg_rip = new TLegend(0.35,0.2,0.99,0.4,"#splitline{M. Ripani et al.,}{#splitline{Phys.Rev. Lett.}{91, 022002 (2003)}}");
TLegend *leg_rip = new TLegend(0.45,0.2,0.99,0.33,"M. Ripani et al., Phys.Rev. Lett. 91, 022002 (2003)");
leg_rip->SetNColumns(2);

/*
  Double_t W  = 0.24;          // Pad Width
   Int_t    Nx = 4;            // Number of pads along X
   Double_t Xm = (1-(Nx*W))/2; // X Margin
   Double_t dw = (W*0.1)/4;
   Double_t rshift = 0.02;

   TPad *p1 = new TPad("p1", "p1", Xm, 0.65, Xm+W+dw+rshift, 0.95, 0, 0, 0);
   p1->SetRightMargin(0);
   p1->Draw();

   TPad *p2 = new TPad("p2", "p2", Xm+W+dw+rshift, 0.65, Xm+2*W-dw+rshift, 0.95, 0, 0, 0);
   p2->SetRightMargin(0);
   p2->SetLeftMargin(0);
   p2->Draw();

   TPad *p3 = new TPad("p3", "p3", Xm+2*W-dw+rshift, 0.65, Xm+3*W+rshift, 0.95, 0, 0, 0);  
   p3->SetRightMargin(0);
   p3->SetLeftMargin(0);
   p3->Draw();
   
   TPad *p4 = new TPad("p4", "p4", Xm+3*W-dw+rshift, 0.65, Xm+4*W+rshift, 0.95, 0, 0, 0);  
   p4->SetLeftMargin(0);
   p4->Draw();   
   
   TPad *p5 = new TPad("p5", "p5", Xm, 0.35, Xm+W+dw+rshift, 0.65, 0, 0, 0);
   p5->SetTopMargin(0);
   p5->SetRightMargin(0);
   p5->SetLeftMargin(0);
   p5->Draw();  
   
   TPad *p6 = new TPad("p6", "p6", Xm+W+dw+rshift, 0.35, Xm+2*W-dw+rshift, 0.65, 0, 0, 0);
   p6->SetRightMargin(0);
   p6->SetLeftMargin(0);
   p6->Draw();

   TPad *p7 = new TPad("p7", "p7", Xm+2*W-dw+rshift, 0.35, Xm+3*W+rshift, 0.65, 0, 0, 0);  
   p7->SetRightMargin(0);
   p7->SetLeftMargin(0);
   p7->Draw();
   
   TPad *p8 = new TPad("p8", "p8", Xm+3*W-dw+rshift, 0.35, Xm+4*W+rshift, 0.65, 0, 0, 0);  
   p8->SetLeftMargin(0);
   p8->Draw();       


   TPad *p9 = new TPad("p9", "p9", Xm, 0.0, Xm+W+dw+rshift, 0.35, 0, 0, 0);
   p9->SetTopMargin(0);
   p9->SetRightMargin(0);
   p9->SetLeftMargin(0);
   p9->Draw();  
   
   TPad *p10 = new TPad("p10", "p10", Xm+W+dw+rshift, 0.0, Xm+2*W-dw+rshift, 0.35, 0, 0, 0);
   p10->SetRightMargin(0);
   p10->SetLeftMargin(0);
   p10->Draw();

   TPad *p11 = new TPad("p11", "p11", Xm+2*W-dw+rshift, 0.0, Xm+3*W+rshift, 0.35, 0, 0, 0);  
   p11->SetRightMargin(0);
   p11->SetLeftMargin(0);
   p11->Draw();
   
   TPad *p12 = new TPad("p12", "p12", Xm+3*W-dw+rshift, 0.0, Xm+4*W+rshift, 0.35, 0, 0, 0);  
   p12->SetLeftMargin(0);
   p12->Draw();   

*/
//Define input files

TFile *file_cr_sec_pim = new TFile("out_cr_sec_all_top_final_bin_corr.root","READ");

TFile *file_max_eff = new TFile("out_cr_sec_all_top_final_bin_corr_eff.root","READ");

TFile *file_EG = new TFile("EG.root","READ");


 for (Int_t qq2=0; qq2<12;qq2++) {
 Q2_bin = 0.425 + 0.05*qq2;
   double ax[21];
   double ay[21];
   double aexl[21];
   double aexh[21];
   double aeyl[21];
   double aeyh[21];
for (Int_t i=0; i<21;i++) {
ax[i] = 1.3125+0.025*i;
ay[i]=0;
aexl[i]=0;
aexh[i]=0;
aeyl[i]=0;
aeyh[i]=0;

};  

qqq.str("");
qqq << "Q^{2} = " << Q2_bin << " GeV^{2}";
//h_w[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),(get_max_w(Q2_bin)-get_min_w(Q2_bin)),1.3+0.025*get_min_w(Q2_bin),1.3+0.025*get_max_w(Q2_bin));

h_w[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21,1.3,1.825);
qqq.str("");

qqq.str("");
qqq << "Q^{2}_{model} = " << Q2_bin << " GeV^{2}";


h_w_model[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21,1.3,1.825);
qqq.str("");

qqq << "Q^{2} = " << Q2_bin << " GeV^{2} err";
h_w_err[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21,1.3,1.825);

qqq.str("");
//for (Int_t i=13; i<14;i++) {
for (Int_t i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {
 W_bin[i] = 1.3125+0.025*i; 

read_data_rec(file_cr_sec_pim,i);
read_data_eff(file_max_eff,i);
read_EG(file_EG,i);

m_pip_pim->Scale(m_pip_pim->GetBinWidth(1));
m_pip_p->Scale(m_pip_p->GetBinWidth(1));
m_pim_p->Scale(m_pim_p->GetBinWidth(1));
alpha_pim->Scale((alpha_pim->GetBinWidth(1))*3.1415/180.);
alpha_p->Scale((alpha_p->GetBinWidth(1))*3.1415/180.);
alpha_pip->Scale((alpha_pip->GetBinWidth(1))*3.1415/180.);

m_pip_pim_bin_corr->Scale(m_pip_pim_bin_corr->GetBinWidth(1));
m_pip_p_bin_corr->Scale(m_pip_p_bin_corr->GetBinWidth(1));
m_pim_p_bin_corr->Scale(m_pim_p_bin_corr->GetBinWidth(1));
alpha_pim_bin_corr->Scale((alpha_pim_bin_corr->GetBinWidth(1))*3.1415/180.);
alpha_p_bin_corr->Scale((alpha_p_bin_corr->GetBinWidth(1))*3.1415/180.);
alpha_pip_bin_corr->Scale((alpha_pip_bin_corr->GetBinWidth(1))*3.1415/180.);

TH1D *h_cos_th;
//FILLING OUT -d(cos) HISTOGRAM
Int_t n_theta_bins;
if ((i==0)||(i==1)) n_theta_bins = 6;
if ((i==2)||(i==3)) n_theta_bins = 8;
if (i>=4) n_theta_bins = 10; 
h_cos_th = new TH1D("h_cos_th","h_cos_th",n_theta_bins,0.,180.);  
Double_t temp;   
for (Int_t j=1; j<=n_theta_bins; j++) {
temp = cos((h_cos_th->GetBinLowEdge(j))*3.1415/180.)-cos(3.1415/180.*(h_cos_th->GetBinLowEdge(j)+h_cos_th->GetBinWidth(j)));
temp = 1./temp;
h_cos_th->SetBinContent(j,temp);
h_cos_th->SetBinError(j,0.);
};
theta_p->Divide(h_cos_th);
theta_pip->Divide(h_cos_th);
theta_pim->Divide(h_cos_th);

theta_p_bin_corr->Divide(h_cos_th);
theta_pip_bin_corr->Divide(h_cos_th);
theta_pim_bin_corr->Divide(h_cos_th);

Double_t w_int;
//w_int = (m_pip_pim->Integral())+(m_pip_p->Integral())+(m_pim_p->Integral());
w_int = (alpha_pim->Integral())+(alpha_p->Integral())+(alpha_pip->Integral());
//w_int = w_int + (theta_pim->Integral())+(theta_p->Integral())+(theta_pip->Integral());

w_int = w_int/3.;

Double_t w_int_eff;
//w_int = (m_pip_pim->Integral())+(m_pip_p->Integral())+(m_pim_p->Integral());
w_int_eff = (alpha_pim_bin_corr->Integral())+(alpha_p_bin_corr->Integral())+(alpha_pip_bin_corr->Integral());
//w_int = w_int + (theta_pim->Integral())+(theta_p->Integral())+(theta_pip->Integral());

w_int_eff = w_int_eff/3.;


//cout << w_int << "  " << w_int_eff << endl;

h_w[qq2]->SetBinContent(i+1,w_int);

m_pip_pim_model->Scale(m_pip_pim_model->GetBinWidth(1));
m_pip_p_model->Scale(m_pip_p_model->GetBinWidth(1));
m_pim_p_model->Scale(m_pim_p_model->GetBinWidth(1));
alpha_pim_model->Scale((alpha_pim_model->GetBinWidth(1))*3.1415/180.);
alpha_p_model->Scale((alpha_p_model->GetBinWidth(1))*3.1415/180.);
alpha_pip_model->Scale((alpha_pip_model->GetBinWidth(1))*3.1415/180.);
Double_t w_int_model;
//w_int_model = (m_pip_pim_model->Integral())+(m_pip_p_model->Integral())+(m_pim_p_model->Integral());
w_int_model = (alpha_pim_model->Integral())+(alpha_p_model->Integral())+(alpha_pip_model->Integral());

//cout << alpha_pim_model->Integral(1,alpha_pim_model->GetNbinsX()) << "   " << alpha_pim->Integral() << endl;

w_int_model = w_int_model/3.;

h_w_model[qq2]->SetBinContent(i+1,w_int_model);




//cout << w_int << "  " << w_int_eff << endl;

h_w[qq2]->SetBinContent(i+1,w_int);
//if (qq2 == 3) cout << m_pip_pim->Integral() << " " << w_int  << " " << w_int_eff << endl;
Double_t err_total;
//m_pip_pim->IntegralAndError(1,m_pip_pim->GetNbinsX(),w_err);
//err_total = w_err*w_err;
//err_total = w_err;
//m_pip_p->IntegralAndError(1,m_pip_p->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;
//m_pim_p->IntegralAndError(1,m_pim_p->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;

alpha_pim->IntegralAndError(1,alpha_pim->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
err_total = w_err;
alpha_p->IntegralAndError(1,alpha_p->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
err_total = err_total + w_err;
alpha_pip->IntegralAndError(1,alpha_pip->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
err_total = err_total + w_err;

//theta_pim->IntegralAndError(1,theta_pim->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;
//theta_p->IntegralAndError(1,theta_p->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;
//theta_pip->IntegralAndError(1,theta_pip->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;

//err_total = sqrt(err_total)/9.;
err_total = err_total/3.;
//cout << err_total << endl;

h_w[qq2]->SetBinError(i+1,err_total);
//draw_1d_canvas(i,qq2);
ax[i] = W_bin[i];
ay[i] = w_int;
aexl[i] = 0.0125;
aexh[i] = 0.0125;
aeyl[i] = w_int*0.05;
aeyh[i] = w_int*0.05;

Double_t sys_err;


h_w_model[qq2]->SetBinError(i+1,0.);


/*
sys_err = ((w_int - m_pip_pim->Integral())*(w_int - m_pip_pim->Integral()));
sys_err = sys_err + ((w_int - m_pip_p->Integral())*(w_int - m_pip_p->Integral()));
sys_err = sys_err + ((w_int - m_pim_p->Integral())*(w_int - m_pim_p->Integral()));
*/
sys_err = ((w_int - alpha_p->Integral())*(w_int - alpha_p->Integral()));
sys_err = sys_err + ((w_int - alpha_pip->Integral())*(w_int - alpha_pip->Integral()));
sys_err = sys_err + ((w_int - alpha_pim->Integral())*(w_int - alpha_pim->Integral()));
sys_err = sys_err + ((w_int - w_int_eff)*(w_int - w_int_eff));

//cout << sqrt(sys_err)/w_int*100. << endl;

/*
sys_err = sys_err + ((w_int - theta_p->Integral())*(w_int - theta_p->Integral()));
sys_err = sys_err + ((w_int - theta_pip->Integral())*(w_int - theta_pip->Integral()));
sys_err = sys_err + ((w_int - theta_pim->Integral())*(w_int - theta_pim->Integral()));
*/
sys_err = sys_err + (w_int*0.03)*(w_int*0.03);

aeyl[i] = aeyl[i]*aeyl[i]+sys_err;
aeyh[i] = aeyh[i]*aeyh[i]+sys_err;

h_w_err[qq2]->SetBinContent(i+1,sqrt(aeyl[i]));

aeyl[i] = aeyl[i]+err_total*err_total;
aeyh[i] = aeyh[i]+err_total*err_total;

aeyl[i] = sqrt(aeyl[i]);
aeyh[i] = sqrt(aeyh[i]);

full_err[i][qq2] = aeyh[i];

};
//c->cd();
//c->cd(qq2+1); 

/*
if (qq2 == 0)p1->cd();
if (qq2 == 1)p2->cd();
if (qq2 == 2)p3->cd();
if (qq2 == 3)p4->cd();
if (qq2 == 4)p5->cd();
if (qq2 == 5)p6->cd();
if (qq2 == 6)p7->cd();
if (qq2 == 7)p8->cd();
if (qq2 == 8)p9->cd();
if (qq2 == 9)p10->cd();
if (qq2 == 10)p11->cd();
if (qq2 == 11)p12->cd();
*/

/*
//if (qq2 < 4)gStyle->SetTitleY(0.);
 TStyle *st1 = new TStyle("st1","my style");
if (qq2 >= 4) {
st1->SetTitleY(0.9);
} else {
st1->SetTitleY(0.5);
};
    gROOT->SetStyle("st1");
*/

 
//c->cd(qq2+1)->SetBottomMargin(0.16);
//c->cd(qq2+1)->SetTopMargin(0.08);
//c->cd(qq2+1)->SetLeftMargin(0.23);
//c->cd(qq2+1)->SetRightMargin(0.01);
gStyle->SetTitleY(0.97);

/*
if ((qq2 == 0) || (qq2 == 4)) {
gPad->SetBottomMargin(small);
gPad->SetRightMargin(small);
gPad->SetLeftMargin(0.17);
gPad->SetTopMargin(small);

};
if (((qq2 >0) && (qq2 < 3)) ||((qq2 >4) && (qq2 < 7)))  {   
   gPad->SetBottomMargin(small);
   gPad->SetRightMargin(small);
   gPad->SetLeftMargin(small);
   gPad->SetTopMargin(small);
};  

if ((qq2 == 3) ||(qq2 == 7))  {   
   gPad->SetBottomMargin(small);
//   gPad->SetRightMargin(small);
   gPad->SetLeftMargin(small);
   gPad->SetTopMargin(small);
}; 

if (qq2 == 8) {
gPad->SetBottomMargin(0.17);
gPad->SetRightMargin(small);
gPad->SetLeftMargin(0.17);
gPad->SetTopMargin(small);

};

if ((qq2 >8) && (qq2 < 11))  {   
   gPad->SetBottomMargin(0.17);
   gPad->SetRightMargin(small);
   gPad->SetLeftMargin(small);
   gPad->SetTopMargin(small);
};  


if (qq2 == 11) {
   gPad->SetBottomMargin(0.17);
//   gPad->SetRightMargin(small);
   gPad->SetLeftMargin(small);
   gPad->SetTopMargin(small);
};
*/
   gPad->SetBottomMargin(0.17);
   gPad->SetLeftMargin(0.17);
   gPad->SetRightMargin(small);
   gPad->SetTopMargin(0.05);

//c->cd(qq2+1)->SetFillStyle(4000);
//gPad->SetBottomMargin(0.001);
//gPad->DrawFrame(1.2,0,1.9,40);
h_w[qq2] ->GetXaxis()->SetTitleSize(0.08);
h_w[qq2] ->GetXaxis()->SetRangeUser(1.2875,1.85);
 h_w[qq2]->GetYaxis()->SetTitleSize(0.08);
 h_w[qq2]->GetYaxis()->SetTitleOffset(1.);
 h_w[qq2]->GetXaxis()->SetLabelSize(0.08);
// h_w[qq2]->GetXaxis()->SetNdivisions(6);
 h_w[qq2]->GetYaxis()->SetLabelSize(0.07);
// h_w[qq2]->GetYaxis()->SetNdivisions(5);
h_w[qq2]->SetMinimum(0.);
h_w[qq2]->SetMaximum(40.);
h_w[qq2]->GetXaxis()->SetTitle("W (GeV)");
h_w[qq2]->GetYaxis()->SetTitle("#sigma (#mub)");
h_w[qq2]->SetMarkerStyle(20);
h_w[qq2]->SetMarkerSize(0.5);
h_w[qq2]->SetMarkerColor(1);
//h_w[qq2]->Draw("e1P");





gae = new TGraphAsymmErrors(21, ax, ay, aexl, aexh, aeyl, aeyh);
qqq.str("");
qqq << "Q^{2} = " << Q2_bin << " GeV^{2}";
gae->SetTitle(qqq.str().c_str());

gae->GetXaxis()->SetTitleSize(0.08);
 gae->GetYaxis()->SetTitleSize(0.08);
 gae->GetYaxis()->SetTitleOffset(1.);
 gae->GetXaxis()->SetLabelSize(0.08);
 gae->GetYaxis()->SetLabelSize(0.08); 
gae->GetXaxis()->SetTitle("W (GeV) ");
gae->GetYaxis()->SetTitle("#sigma (#mub) ");   
   gae->SetFillColor(21);
//   gae->SetFillStyle(1); 
   gae->SetMinimum(0.2);
   gae->SetMaximum(40.);
   gae->GetXaxis()->SetRangeUser(1.2875,1.85);

//gae->SetMarkerStyle(20);
//gae->SetMarkerSize(0.7);
//gae->GetXaxis()->ChangeLabel(1,-1.,-1.,1,-1,-1,"1.3");
//gae->GetXaxis()->ChangeLabel(2,-1.,-1.,1,-1,-1,"1.4");
//gae->GetXaxis()->ChangeLabel(3,-1.,-1.,1,-1,-1,"1.5");
//gae->GetXaxis()->ChangeLabel(4,-1.,-1.,1,-1,-1,"1.6");
//gae->GetXaxis()->ChangeLabel(5,-1.,-1.,1,-1,-1,"1.7");
if (qq2 != 8) {
gae->GetXaxis()->ChangeLabel(1,-1.,-1.,-1,-1,-1," ");
} else{
gae->GetXaxis()->ChangeLabel(1,-1.,-1.,-1,-1,-1,"1.3");
};
gae->GetYaxis()->ChangeLabel(8,-1.,-1.,-1,-1,-1," ");
//gae->Draw("a2");
//gae->Draw("AP0");
h_w[qq2]->SetMarkerStyle(20);
h_w[qq2]->SetMarkerSize(0.5);
h_w[qq2]->SetLineColor(1);
h_w[qq2]->SetTitle("");

if(qq2==0)h_w[qq2]->SetMarkerStyle(20); 
if(qq2==0)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==0)leg->AddEntry(h_w[qq2],"Q^{2} = 0.425 GeV^{2}","p");

if(qq2==1)h_w[qq2]->SetMarkerStyle(21); 
if(qq2==1)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==1)leg->AddEntry(h_w[qq2],"Q^{2} = 0.475 GeV^{2}","p");

if(qq2==2)h_w[qq2]->SetMarkerStyle(22); 
if(qq2==2)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==2)leg->AddEntry(h_w[qq2],"Q^{2} = 0.525 GeV^{2}","p");

if(qq2==3)h_w[qq2]->SetMarkerStyle(23); 
if(qq2==3)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==3)leg->AddEntry(h_w[qq2],"Q^{2} = 0.575 GeV^{2}","p");

if(qq2==4)h_w[qq2]->SetMarkerStyle(24); 
if(qq2==4)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==4)leg->AddEntry(h_w[qq2],"Q^{2} = 0.625 GeV^{2}","p");

if(qq2==5)h_w[qq2]->SetMarkerStyle(25); 
if(qq2==5)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==5)leg->AddEntry(h_w[qq2],"Q^{2} = 0.675 GeV^{2}","p");

if(qq2==6)h_w[qq2]->SetMarkerStyle(26); 
if(qq2==6)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==6)leg->AddEntry(h_w[qq2],"Q^{2} = 0.725 GeV^{2}","p");

if(qq2==7)h_w[qq2]->SetMarkerStyle(27); 
if(qq2==7)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==7)leg->AddEntry(h_w[qq2],"Q^{2} = 0.775 GeV^{2}","p");

if(qq2==8)h_w[qq2]->SetMarkerStyle(28); 
if(qq2==8)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==8)leg->AddEntry(h_w[qq2],"Q^{2} = 0.825 GeV^{2}","p");

if(qq2==9)h_w[qq2]->SetMarkerStyle(29); 
if(qq2==9)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==9)leg->AddEntry(h_w[qq2],"Q^{2} = 0.875 GeV^{2}","p");

if(qq2==10)h_w[qq2]->SetMarkerStyle(30); 
if(qq2==10)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==10)leg->AddEntry(h_w[qq2],"Q^{2} = 0.925 GeV^{2}","p");

if(qq2==11)h_w[qq2]->SetMarkerStyle(31); 
if(qq2==11)h_w[qq2]->SetMarkerSize(1.2);
if(qq2==11)leg->AddEntry(h_w[qq2],"Q^{2} = 0.975 GeV^{2}","p");


if (qq2 == 0) {
h_w[qq2]->Draw("e1P");
} else {
h_w[qq2]->Draw("e1P same");
};

//h_w_err[qq2]->SetMarkerStyle(20);
//h_w_err[qq2]->SetMarkerSize(2);
h_w_err[qq2]->SetLineColor(1);
h_w_err[qq2]->SetFillColor(1);
//h_w_err[qq2]->Draw("same ][");


h_w_model[qq2]->SetLineWidth(2.);
h_w_model[qq2]->SetLineColor(kBlack);
Double_t factor = ((h_w[qq2]->Integral())/(h_w_model[qq2]->Integral()));
cout << "factor = " << factor << " Q2 = " << Q2_bin << endl;
factor = 1.;
h_w_model[qq2]->Scale((h_w[qq2]->Integral())/(h_w_model[qq2]->Integral()));
//cout << "factor[" << qq2 << "]=" <<  (h_w[qq2]->Integral())/(h_w_model[qq2]->Integral()) << ";" << endl;

if(qq2==0)h_w_model[0]->GetXaxis()->SetRangeUser(1.6,1.8125);
if(qq2==2)h_w_model[2]->GetXaxis()->SetRangeUser(1.3,1.7875);
if(qq2==3)h_w_model[3]->GetXaxis()->SetRangeUser(1.3,1.7875);
if(qq2==4)h_w_model[4]->GetXaxis()->SetRangeUser(1.3,1.7625);
if(qq2==5)h_w_model[5]->GetXaxis()->SetRangeUser(1.3,1.7375);
if(qq2==6)h_w_model[6]->GetXaxis()->SetRangeUser(1.3,1.7125);
if(qq2==7)h_w_model[7]->GetXaxis()->SetRangeUser(1.3,1.6875);
if(qq2==8)h_w_model[8]->GetXaxis()->SetRangeUser(1.3,1.6375);
if(qq2==9)h_w_model[9]->GetXaxis()->SetRangeUser(1.3,1.6125);
if(qq2==10)h_w_model[10]->GetXaxis()->SetRangeUser(1.3,1.5875);
if(qq2==11)h_w_model[11]->GetXaxis()->SetRangeUser(1.3,1.5375);

//h_w_model[qq2]->Draw("same C");












if (qq2==0) {
Double_t w_values1[9]={1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={3.07158, 5.02006, 7.2034, 11.2869, 14.1014, 17.9653, 20.5538, 27.2129, 30.0085, 32.6551, 32.1707, 30.9483, 29.3869, 25.357, 19.3566, 13.8919, 10.9587};

//Double_t data_int_mod_res[17] ={1.9264, 3.09148, 4.47598, 7.02632, 8.37201, 9.02445, 8.46081, 8.62866, 9.19142, 10.7838, 12.3823, 14.8718, 15.8047, 14.764, 12.3983, 9.8233, 8.34427};

Double_t data_int_mod_res1[9] ={9.19142, 10.7838, 12.3823, 14.8718, 15.8047, 14.764, 12.3983, 9.8233, 8.34427};

Double_t data_int_mod_res1_trans[9] ={6.90167, 8.53918, 10.1375, 12.6247, 13.6339, 13.027, 11.2075, 9.06667, 7.84885};




TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(9,w_values1 ,data_int_mod_res1);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");



TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(9,w_values1 ,data_int_mod_res1_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res_trans->Draw("same");




};






if (qq2==1) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={2.64723, 4.29793, 6.17725, 9.69122, 12.0081, 14.8764, 16.5312, 21.2103, 23.3065, 25.696, 25.9832, 26.2029, 25.6492, 22.4934, 17.3137, 12.5073, 9.97956};

Double_t data_int_mod_res[17] ={1.84033, 2.93754, 4.25271, 6.68189, 7.95812, 8.54494, 7.95527, 8.02193, 8.50501, 10.0102, 11.5962, 14.18, 15.2281, 14.2865, 11.9904, 9.43985, 8.03662};

Double_t data_int_mod_res_trans[17] ={0.633473, 1.27918, 2.05641, 3.66853, 4.68618, 5.49442, 5.37161, 5.67707, 6.21956, 7.75531, 9.33405, 11.8792, 12.9781, 12.4686, 10.7408, 8.65969, 7.52828};


TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(17,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res->Draw("same");


TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(17,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");


};








if (qq2==2) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={2.34137, 3.77058, 5.42355, 8.51528, 10.4781, 12.6646, 13.6885, 17.0021, 18.5783, 20.6855, 21.4031, 22.5324, 22.6738, 20.1758, 15.6592, 11.3803, 9.18491};

//Double_t data_int_mod_res[16] ={2.34137, 3.77058, 5.42355, 8.51528, 10.4781, 12.6646, 13.6885, 17.0021, 18.5783, 20.6855, 21.4031, 22.5324, 22.6738, 20.1758, 15.6592, 11.3803};

Double_t data_int_mod_res[16] ={1.76479, 2.79763, 4.04651, 6.35905, 7.5726, 8.11025, 7.50851, 7.4954, 7.91057, 9.34207, 10.9165, 13.5849, 14.7241, 13.8644, 11.6367, 9.10399};

Double_t data_int_mod_res_trans[16] ={0.583919, 1.19727, 1.94301, 3.49912, 4.47155, 5.19842, 5.00572, 5.17408, 5.63386, 7.07674, 8.62905, 11.2023, 12.3602, 11.9327, 10.3004, 8.28606};




TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(16,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");


TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(17,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");



};






if (qq2==3) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={2.1135, 3.37078, 4.84748, 7.61084, 9.31088, 11.0154, 11.6025, 13.946, 15.1302, 16.9752, 17.9396, 19.6721, 20.3047, 18.3081, 14.3306, 10.4695, 8.54533};

//Double_t data_int_mod_res[16] ={2.1135, 3.37078, 4.84748, 7.61084, 9.31088, 11.0154, 11.6025, 13.946, 15.1302, 16.9752, 17.9396, 19.6721, 20.3047, 18.3081, 14.3306, 10.4695};


Double_t data_int_mod_res[16] ={1.78712, 2.79597, 4.01087, 6.25888, 7.45211, 8.02306, 7.48197, 7.53435, 8.01827, 9.68565, 11.5884, 15.2591, 17.4502, 16.4844, 13.0421, 9.64917};

Double_t data_int_mod_res_trans[16] ={0.544801, 1.1274, 1.84308, 3.3442, 4.27507, 4.93225, 4.68099, 4.72755, 5.10866, 6.44912, 7.94996, 10.4812, 11.633, 11.2686, 9.77547, 7.87635};


TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(16,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");


TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(16,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");


};








if (qq2==4) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={1.94107, 3.06088, 4.39583, 6.89446, 8.39278, 9.74936, 10.031, 11.6736, 12.5616, 14.1845, 15.2957, 17.4439, 18.4216, 16.8059, 13.2689, 9.73721, 8.0334};

//Double_t data_int_mod_res[15] ={1.94107, 3.06088, 4.39583, 6.89446, 8.39278, 9.74936, 10.031, 11.6736, 12.5616, 14.1845, 15.2957, 17.4439, 18.4216, 16.8059, 13.2689};


Double_t data_int_mod_res[15] ={1.70933, 2.65116, 3.79774, 5.92327, 7.05253, 7.57866, 7.03202, 7.00709, 7.41049, 8.94891, 10.763, 14.3281, 16.4234, 15.5425, 12.3873};

Double_t data_int_mod_res_trans[15] ={0.518026, 1.07273, 1.76092, 3.20986, 4.10363, 4.70453, 4.40902, 4.35745, 4.67321, 5.92713, 7.38066, 9.86781, 11.0074, 10.696, 9.32607};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(15,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");



TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(15,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");




};




if (qq2==5) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={1.81058, 2.8181, 4.03632, 6.31531, 7.65434, 8.75588, 8.82375, 9.95462, 10.6192, 12.0608, 13.2584, 15.6832, 16.8864, 15.5582, 12.3907, 9.13265, 7.61308};

//Double_t data_int_mod_res[17] ={1.66996, 2.5727, 3.67779, 5.72491, 6.82847, 7.38529, 6.9075, 6.97216, 7.4556, 9.20648, 11.129, 14.4894, 16.1412, 15.0591, 12.0162, 8.88832, 7.43579};


Double_t data_int_mod_res[14] ={1.66996, 2.5727, 3.67779, 5.72491, 6.82847, 7.38529, 6.9075, 6.97216, 7.4556, 9.20648, 11.129, 14.4894, 16.1412, 15.0591};

Double_t data_int_mod_res_trans[14] ={0.503173, 1.03257, 1.69561, 3.09445, 3.95424, 4.50948, 4.18216, 4.05331, 4.31489, 5.49522, 6.90494, 9.34703, 10.4698, 10.2029};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(14,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");



TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(14,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");




};



if (qq2==6) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={1.71295, 2.62694, 3.74691, 5.83872, 7.04814, 7.95922, 7.87786, 8.63096, 9.12473, 10.4133, 11.6484, 14.2148, 15.5352, 14.4269, 11.5888, 8.59326, 7.24104};

Double_t data_int_mod_res[13] ={1.7127, 2.6265, 3.74627, 5.83767, 7.04671, 7.95698, 7.87484, 8.62605, 9.11793, 10.4016, 11.6315, 14.1883, 15.5081};

Double_t data_int_mod_res_trans[13] ={0.499818, 1.00618, 1.64618, 2.99639, 3.82426, 4.34203, 3.99351, 3.8053, 4.02156, 5.13808, 6.50654, 8.90314, 10.0059};


//Double_t data_int_mod_res[17] ={1.7127, 2.6265, 3.74627, 5.83767, 7.04671, 7.95698, 7.87484, 8.62605, 9.11793, 10.4016, 11.6315, 14.1883, 15.5081, 14.4008, 11.5623, 8.57487, 7.22601};

TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(13,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");



TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(13,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");






};



if (qq2==7) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={1.6413, 2.47586, 3.5114, 5.43957, 6.5401, 7.30539, 7.12002, 7.59014, 7.94959, 9.09923, 10.3263, 12.9011, 14.246, 13.313, 10.7877, 8.07464, 6.88746};

Double_t data_int_mod_res[12] ={1.6413, 2.47586, 3.5114, 5.43957, 6.5401, 7.30539, 7.12002, 7.59014, 7.94959, 9.09923, 10.3263, 12.9011};

Double_t data_int_mod_res_trans[12] ={0.506814, 0.991575, 1.60985, 2.91125, 3.70733, 4.19145, 3.82829, 3.59035, 3.76246, 4.80754, 6.12133, 8.44408};

//Double_t data_int_mod_res[17] ={1.6413, 2.47586, 3.5114, 5.43957, 6.5401, 7.30539, 7.12002, 7.59014, 7.94959, 9.09923, 10.3263, 12.9011, 14.246, 13.313, 10.7877, 8.07464, 6.88746};

TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(12,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");






TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(12,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");








};



if (qq2==8) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={1.59047, 2.35652, 3.31844, 5.10099, 6.10773, 6.75907, 6.50226, 6.75916, 7.01152, 8.03277, 9.21851, 11.7036, 13.0081, 12.2204, 9.99384, 7.57735, 6.55313};

Double_t data_int_mod_res[10] ={1.59047, 2.35652, 3.31844, 5.10099, 6.10773, 6.75907, 6.50226, 6.75916, 7.01152, 8.03277};

Double_t data_int_mod_res_trans[10] ={0.524378, 0.989176, 1.58735, 2.84039, 3.60512, 4.05966, 3.6887, 3.41249, 3.54444, 4.51984};


//Double_t data_int_mod_res[17] ={1.6413, 2.47586, 3.5114, 5.43957, 6.5401, 7.30539, 7.12002, 7.59014, 7.94959, 9.09923, 10.3263, 12.9011, 14.246, 13.313, 10.7877, 8.07464, 6.88746};


TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(10,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");





TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(10,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");









};


if (qq2==9) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={1.55666, 2.26314, 3.16047, 4.81228, 5.73703, 6.29845, 5.9948, 6.09279, 6.26104, 7.16944, 8.29923, 10.6505, 11.8823, 11.2146, 9.26174, 7.12728, 6.25477};

Double_t data_int_mod_res[9] ={1.55666, 2.26314, 3.16047, 4.81228, 5.73703, 6.29845, 5.9948, 6.09279, 6.26104};

Double_t data_int_mod_res_trans[9] ={0.552325, 0.998773, 1.57851, 2.78385, 3.51759, 3.94601, 3.57332, 3.26961, 3.36596};


//Double_t data_int_mod_res[17] ={1.55666, 2.26314, 3.16047, 4.81228, 5.73703, 6.29845, 5.9948, 6.09279, 6.26104, 7.16944, 8.29923, 10.6505, 11.8823, 11.2146, 9.26174, 7.12728, 6.25477};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(9,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");







TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(9,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");


















};



if (qq2==10) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={1.53701, 2.19138, 3.03188, 4.56564, 5.41792, 5.90783, 5.57599, 5.55767, 5.66121, 6.47539, 7.54795, 9.76091, 10.9108, 10.3405, 8.62711, 6.73985, 6.00117};

Double_t data_int_mod_res[8] ={1.53701, 2.19138, 3.03188, 4.56564, 5.41792, 5.90783, 5.57599, 5.55767};

Double_t data_int_mod_res_trans[8] ={0.590147, 1.0197, 1.58259, 2.74074, 3.44355, 3.84815, 3.47834, 3.15556};


//Double_t data_int_mod_res[17] ={1.53701, 2.19138, 3.03188, 4.56564, 5.41792, 5.90783, 5.57599, 5.55767, 5.66121, 6.47539, 7.54795, 9.76091, 10.9108, 10.3405, 8.62711, 6.73985, 6.00117};

TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(8,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");



TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(8,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");





};



if (qq2==11) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};

//Double_t data_int_mod_res[17] ={1.52924, 2.13774, 2.92806, 4.35458, 5.142, 5.57416, 5.22772, 5.1251, 5.17893, 5.91548, 6.93444, 9.02027, 10.09, 9.59839, 8.08973, 6.41172, 5.78866};

//Double_t data_int_mod_res[6] ={1.52924, 2.13774, 2.92806, 4.35458, 5.142, 5.57416};
Double_t data_int_mod_res[6] ={1.52924, 2.13774, 2.92806, 4.35458, 5.142, 5.57416};

Double_t data_int_mod_res_trans[6] ={0.63721, 1.05119, 1.59873, 2.71021, 3.38199, 3.76413};

//Double_t data_int_mod_res[17] ={1.52924, 2.13774, 2.92806, 4.35458, 5.142, 5.57416, 5.22772, 5.1251, 5.17893, 5.91548, 6.93444, 9.02027, 10.09, 9.59839, 8.08973, 6.41172, 5.78866};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(6,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");


TGraph *gr_int_mod_res_trans;
gr_int_mod_res_trans = new TGraph(6,w_values ,data_int_mod_res_trans);
gr_int_mod_res_trans->SetMarkerStyle(20);
gr_int_mod_res_trans->SetLineWidth(2);
gr_int_mod_res_trans->SetLineStyle(4);
gr_int_mod_res_trans->SetMarkerColor(kBlack);


//  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res_trans->GetN(); iiii++){     gr_int_mod_res_trans->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res_trans->SetPoint(iiii, x, y);

  } 
    
//gr_int_mod_res_trans->Draw("same");




};















};




std::ostringstream file2;   
file2 << "ripani065.dat";
 cout << file2.str() << "\n";




 float sigma3[30],error3[30],theta3[30],errorx3[30];
 
 

std::ifstream sect3(file2.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (Int_t i= 0; i <= 19; i++)
{
  sect3 >> theta3[i] >> sigma3[i] >> error3[i] >> errorx3[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta3[i] << "  " << sigma3[i] << "  " << error3[i] << "\n";

}
 sect3.close();   
 
   TGraphErrors *Tg2 = new TGraphErrors(20,theta3,sigma3,errorx3,error3);
   
   Tg2->SetMarkerStyle(20);
   Tg2->SetMarkerSize(1.5);
   Tg2->SetMarkerColor(2);
   Tg2->SetLineColor(2); 
   Tg2->SetMaximum(33.);
   Tg2->GetXaxis()->SetLimits(1.3,1.9);
   Tg2->SetMinimum(0.);
   

   Tg2->Draw("P same");
   
   


std::ostringstream file1;   
file1 << "ripani095.dat";
 cout << file1.str() << "\n";




 float sigma2[30],error2[30],theta2[30],errorx2[30];
 
 

std::ifstream sect2(file1.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (Int_t i= 0; i <= 19; i++)
{
  sect2 >> theta2[i] >> sigma2[i] >> error2[i] >> errorx2[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta2[i] << "  " << sigma2[i] << "  " << error2[i] << "\n";

}
 sect2.close();   
 
   TGraphErrors *Tg1 = new TGraphErrors(20,theta2,sigma2,errorx2,error2);
   
   Tg1->SetMarkerStyle(20);
   Tg1->SetMarkerSize(1.5);
   Tg1->SetMarkerColor(3);
   Tg1->SetLineColor(3); 
   Tg1->SetMaximum(33.);
   Tg1->GetXaxis()->SetLimits(1.3,1.9);
   Tg1->SetMinimum(0.);
   

   Tg1->Draw("P same");





    

leg->Draw();

leg_rip->AddEntry(Tg2,"Q^{2} = 0.65 GeV^{2}","p");
leg_rip->AddEntry(Tg1,"Q^{2} = 0.95 GeV^{2}","p");

leg_rip->Draw();



}; //end of main program
