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

TLegend *leg;

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


f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

tmp = new TH1D("tmp","tmp",100,h->GetBinLowEdge(0),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1));
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

c->Divide(3,3);

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
//m_pip_p->Draw("e1");
//m_pip_p->Draw("e0 same");


c->cd(2);
m_pip_pim->SetLineColor(1);
m_pip_pim->SetMarkerColor(1);
m_pip_pim->SetMinimum(0);
m_pip_pim->SetMaximum(1.2*max);
//m_pip_pim->Draw("e1");
//m_pip_pim->Draw("e0 same");


c->cd(3);
m_pim_p->SetLineColor(1);
m_pim_p->SetMarkerColor(1);
m_pim_p->SetMinimum(0);
m_pim_p->SetMaximum(1.2*max);
//m_pim_p->Draw("e1");
//m_pim_p->Draw("e0 same");


max = theta_p->GetMaximum();

if ((theta_pip->GetMaximum()) > max) max = theta_pip->GetMaximum();

if ((theta_pim->GetMaximum()) > max) max = theta_pim->GetMaximum();

c->cd(4);
theta_p->SetLineColor(1);
theta_p->SetMarkerColor(1);
theta_p->SetMinimum(0);
theta_p->SetMaximum(1.2*max);
//theta_p->GetXaxis()->SetTitle("#theta_{p'} (deg)");
//theta_p->Draw("e1");


c->cd(5);
theta_pim->SetLineColor(1);
theta_pim->SetMarkerColor(1);
theta_pim->SetMinimum(0);
theta_pim->SetMaximum(1.2*max);
//theta_pim->GetXaxis()->SetTitle("#theta_{#pi^{-}} (deg)");
//theta_pim->Draw("e1");

c->cd(6);
theta_pip->SetLineColor(1);
theta_pip->SetMarkerColor(1);
theta_pip->SetMinimum(0);
theta_pip->SetMaximum(1.2*max);
theta_pip->GetXaxis()->SetTitle("#theta_{#pi^{+}} (deg)");
//theta_pip->Draw("e1");

max = alpha_p->GetMaximum();

if ((alpha_pip->GetMaximum()) > max) max = alpha_pip->GetMaximum();

if ((alpha_pim->GetMaximum()) > max) max = alpha_pim->GetMaximum();


c->cd(7);
alpha_p->SetLineColor(1);
alpha_p->SetMarkerColor(1);
alpha_p->SetMinimum(0);
alpha_p->SetMaximum(1.2*max);
//alpha_p->Draw("e1");


c->cd(8);
alpha_pim->SetLineColor(1);
alpha_pim->SetMarkerColor(1);
alpha_pim->SetMinimum(0);
alpha_pim->SetMaximum(1.2*max);
alpha_pim->GetXaxis()->SetTitle("#alpha_{#pi^{-}} (deg)");
//alpha_pim->Draw("e1");


c->cd(9);
alpha_pip->SetLineColor(1);
alpha_pip->SetMarkerColor(1);
alpha_pip->SetMinimum(0);
alpha_pip->SetMaximum(1.2*max);
//alpha_pip->Draw("e1");


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




void res_cont() {
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




c->Divide(4,3);


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
c->cd();
c->cd(qq2+1); 
c->cd(qq2+1)->SetBottomMargin(0.16);
c->cd(qq2+1)->SetTopMargin(0.08);
c->cd(qq2+1)->SetLeftMargin(0.23);
c->cd(qq2+1)->SetRightMargin(0.01);
c->cd(qq2+1)->SetFillStyle(4000);
h_w[qq2] ->GetXaxis()->SetTitleSize(0.08);
h_w[qq2] ->GetXaxis()->SetRangeUser(1.2875,1.85);
 h_w[qq2]->GetYaxis()->SetTitleSize(0.08);
 h_w[qq2]->GetYaxis()->SetTitleOffset(1.3);
 h_w[qq2]->GetXaxis()->SetLabelSize(0.08);
// h_w[qq2]->GetXaxis()->SetNdivisions(6);
 h_w[qq2]->GetYaxis()->SetLabelSize(0.07);
// h_w[qq2]->GetYaxis()->SetNdivisions(5);
h_w[qq2]->SetMinimum(0.);
h_w[qq2]->SetMaximum(40.);
h_w[qq2]->GetXaxis()->SetTitle("W (GeV)");
h_w[qq2]->GetYaxis()->SetTitle("#sigma (#mubn)");
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
 gae->GetYaxis()->SetTitleOffset(1.3);
 gae->GetXaxis()->SetLabelSize(0.08);
 gae->GetYaxis()->SetLabelSize(0.08); 
gae->GetXaxis()->SetTitle("W (GeV)");
gae->GetYaxis()->SetTitle("#sigma (#mubn)");   
//   gae->SetFillColor(1);
   gae->SetFillStyle(3005); 
   gae->SetMinimum(0.2);
   gae->SetMaximum(40.);
   gae->GetXaxis()->SetRangeUser(1.2875,1.85);

gae->SetMarkerStyle(20);
gae->SetMarkerSize(0.6);
//gae->Draw("a2");
//gae->Draw("AP0");
//h_w[qq2]->Draw("e1P same");

//h_w_err[qq2]->SetMarkerStyle(20);
//h_w_err[qq2]->SetMarkerSize(2);
h_w_err[qq2]->SetLineColor(1);
h_w_err[qq2]->SetFillColor(1);
//h_w_err[qq2]->Draw("][");


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

//h_w_model[qq2]->Draw("][");












if (qq2==0) {
Double_t q2gr1617[8];
q2gr1617[0] = 0.425;
Double_t res1617[8];
Double_t res1718[4];
Double_t w_values1[9]={1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};


//Double_t data_int_mod_res1[9] ={9.19142, 10.7838, 12.3823, 14.8718, 15.8047, 14.764, 12.3983, 9.8233, 8.34427};

Double_t data_int_mod_res1[9] ={15.0997, 17.3058, 18.7788, 21.2375, 22.1427, 20.2333, 16.1943, 12.0668, 9.8289};

TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(9,w_values1 ,data_int_mod_res1);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+13,y);
  } 
  
  
  
  
  

gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825); 
res1617[0]  = resvl->GetBinContent(4);
res1718[0]  = resvl->GetBinContent(5);
cout << "Q2 = 0.425; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
cout << "Q2 = 0.425; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
cout << "Q2 = 0.425; W = 1.6-1.7; ratio = " <<  resvl->GetBinContent(4) << endl;
cout << "Q2 = 0.425; W = 1.7-1.8; ratio = " <<  resvl->GetBinContent(5) << endl;

 resvl->Draw("][");

};






if (qq2==1) {
Double_t q2gr[11];
q2gr[0] = 0.475;
Double_t res1415[11];
Double_t res1516[10];


Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};


//Double_t data_int_mod_res[17] ={1.84033, 2.93754, 4.25271, 6.68189, 7.95812, 8.54494, 7.95527, 8.02193, 8.50501, 10.0102, 11.5962, 14.18, 15.2281, 14.2865, 11.9904, 9.43985, 8.03662};

Double_t data_int_mod_res[17] ={2.16086, 3.46886, 4.99628, 7.83298, 9.49867, 10.943, 11.2003, 13.0055, 14.1315, 16.1338, 17.4996, 19.8728, 20.8261, 19.0478, 15.2064, 11.3073, 9.24795};


TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(17,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);
      

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825); 
cout << "Q2 = 0.475; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
res1415[0] = resvl->GetBinContent(2);
res1516[0] = resvl->GetBinContent(3);
q2gr1617[1] = 0.475;
res1617[1] = resvl->GetBinContent(4);
res1718[1] = resvl->GetBinContent(5);
cout << "Q2 = 0.475; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
cout << "Q2 = 0.475; W = 1.6-1.7; ratio = " <<  resvl->GetBinContent(4) << endl;
cout << "Q2 = 0.475; W = 1.7-1.8; ratio = " <<  resvl->GetBinContent(5) << endl;
 resvl->Draw("][");

};








if (qq2==2) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};





//Double_t data_int_mod_res[16] ={1.76479, 2.79763, 4.04651, 6.35905, 7.5726, 8.11025, 7.50851, 7.4954, 7.91057, 9.34207, 10.9165, 13.5849, 14.7241, 13.8644, 11.6367, 9.10399};

Double_t data_int_mod_res[16] ={2.06286, 3.29445, 4.74447, 7.44534, 9.03217, 10.3956, 10.6107, 12.2636, 13.2648, 15.0873, 16.3617, 18.6756, 19.6784, 18.0178, 14.3528, 10.6465};




TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(16,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);    
      resvl->SetBinContent(iiii+5,y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");

resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825);  
cout << "Q2 = 0.525; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[1] =  0.525;
res1415[1] = resvl->GetBinContent(2);
res1516[1] = resvl->GetBinContent(3);
q2gr1617[2] = 0.525;
res1617[2] = resvl->GetBinContent(4);
res1718[2] = resvl->GetBinContent(5);
cout << "Q2 = 0.525; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
cout << "Q2 = 0.525; W = 1.6-1.7; ratio = " <<  resvl->GetBinContent(4) << endl;
cout << "Q2 = 0.525; W = 1.7-1.8; ratio = " <<  resvl->GetBinContent(5) << endl;
 resvl->Draw("][");

};






if (qq2==3) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};



//Double_t data_int_mod_res[16] ={1.78712, 2.79597, 4.01087, 6.25888, 7.45211, 8.02306, 7.48197, 7.53435, 8.01827, 9.68565, 11.5884, 15.2591, 17.4502, 16.4844, 13.0421, 9.64917};

Double_t data_int_mod_res[16] ={1.97661, 3.13607, 4.51238, 7.08245, 8.59653, 9.89394, 10.0807, 11.6039, 12.4953, 14.164, 15.3653, 17.6502, 18.7038, 17.1459, 13.6355, 10.0852};


TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(16,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();
  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);

  } 
  
  
  
  
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
 resvl->Rebin(4);
 resvl->Scale(0.25); 
resvl->GetXaxis()->SetRangeUser(1.41,1.825); 
cout << "Q2 = 0.575; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[2] =  0.575;
res1415[2] = resvl->GetBinContent(2);
res1516[2] = resvl->GetBinContent(3);
q2gr1617[3] = 0.575;
res1617[3] = resvl->GetBinContent(4);
res1718[3] = resvl->GetBinContent(5);
cout << "Q2 = 0.575; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
cout << "Q2 = 0.575; W = 1.6-1.7; ratio = " <<  resvl->GetBinContent(4) << endl;
cout << "Q2 = 0.575; W = 1.7-1.8; ratio = " <<  resvl->GetBinContent(5) << endl;
 resvl->Draw("][");  

//gr_int_mod_res->Draw("same");

};








if (qq2==4) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};


//Double_t data_int_mod_res[15] ={1.70933, 2.65116, 3.79774, 5.92327, 7.05253, 7.57866, 7.03202, 7.00709, 7.41049, 8.94891, 10.763, 14.3281, 16.4234, 15.5425, 12.3873};

Double_t data_int_mod_res[15] ={1.90223, 2.99408, 4.30056, 6.74485, 8.19127, 9.43507, 9.60588, 11.0206, 11.8172, 13.3567, 14.5019, 16.7822, 17.8814, 16.4101, 13.0363};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(15,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();


  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);
  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825);  
cout << "Q2 = 0.625; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[3] =  0.625;
res1415[3] = resvl->GetBinContent(2);
res1516[3] = resvl->GetBinContent(3);
q2gr1617[4] = 0.625;
res1617[4] = resvl->GetBinContent(4);
cout << "Q2 = 0.625; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
cout << "Q2 = 0.625; W = 1.6-1.7; ratio = " <<  resvl->GetBinContent(4) << endl;
 resvl->Draw("][");
};




if (qq2==5) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};


//Double_t data_int_mod_res[14] ={1.66996, 2.5727, 3.67779, 5.72491, 6.82847, 7.38529, 6.9075, 6.97216, 7.4556, 9.20648, 11.129, 14.4894, 16.1412, 15.0591};

Double_t data_int_mod_res[14] ={1.81058, 2.8181, 4.03632, 6.31531, 7.65434, 8.75588, 8.82375, 9.95462, 10.6192, 12.0608, 13.2584, 15.6832, 16.8864, 15.5582};


TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(14,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);
  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825);  
cout << "Q2 = 0.675; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[4] =  0.675;
res1415[4] = resvl->GetBinContent(2);
res1516[4] = resvl->GetBinContent(3);
q2gr1617[5] = 0.675;
res1617[5] = resvl->GetBinContent(4);
cout << "Q2 = 0.675; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
cout << "Q2 = 0.675; W = 1.6-1.7; ratio = " <<  resvl->GetBinContent(4) << endl;
 resvl->Draw("][");
};



if (qq2==6) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};



//Double_t data_int_mod_res[13] ={1.7127, 2.6265, 3.74627, 5.83767, 7.04671, 7.95698, 7.87484, 8.62605, 9.11793, 10.4016, 11.6315, 14.1883, 15.5081};

Double_t data_int_mod_res[13] ={1.71295, 2.62694, 3.74691, 5.83872, 7.04814, 7.95922, 7.87786, 8.63096, 9.12473, 10.4133, 11.6484, 14.2148, 15.5352};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(13,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);
  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825);  
cout << "Q2 = 0.725; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[5] =  0.725;
res1415[5] = resvl->GetBinContent(2);
res1516[5] = resvl->GetBinContent(3);
q2gr1617[6] = 0.725;
res1617[6] = resvl->GetBinContent(4);
cout << "Q2 = 0.725; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
cout << "Q2 = 0.725; W = 1.6-1.7; ratio = " <<  resvl->GetBinContent(4) << endl;
 resvl->Draw("][");
};



if (qq2==7) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};



//Double_t data_int_mod_res[12] ={1.6413, 2.47586, 3.5114, 5.43957, 6.5401, 7.30539, 7.12002, 7.59014, 7.94959, 9.09923, 10.3263, 12.9011};

Double_t data_int_mod_res[12] ={1.6413, 2.47586, 3.5114, 5.43957, 6.5401, 7.30539, 7.12002, 7.59014, 7.94959, 9.09923, 10.3263, 12.9011};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(12,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825);  
cout << "Q2 = 0.775; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[6] =  0.775;
res1415[6] = resvl->GetBinContent(2);
res1516[6] = resvl->GetBinContent(3);
q2gr1617[7] = 0.775;
res1617[7] = resvl->GetBinContent(4);
cout << "Q2 = 0.775; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
cout << "Q2 = 0.775; W = 1.6-1.7; ratio = " <<  resvl->GetBinContent(4) << endl;
 resvl->Draw("][");
};



if (qq2==8) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};



//Double_t data_int_mod_res[10] ={1.59047, 2.35652, 3.31844, 5.10099, 6.10773, 6.75907, 6.50226, 6.75916, 7.01152, 8.03277};

Double_t data_int_mod_res[10] ={1.59047, 2.35652, 3.31844, 5.10099, 6.10773, 6.75907, 6.50226, 6.75916, 7.01152, 8.03277};


TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(10,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);

  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825);  
cout << "Q2 = 0.825; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[7] =  0.825;
res1415[7] = resvl->GetBinContent(2);
res1516[7] = resvl->GetBinContent(3);
cout << "Q2 = 0.825; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
 resvl->Draw("][");
};


if (qq2==9) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};



//Double_t data_int_mod_res[9] ={1.55666, 2.26314, 3.16047, 4.81228, 5.73703, 6.29845, 5.9948, 6.09279, 6.26104};

Double_t data_int_mod_res[9] ={1.55666, 2.26314, 3.16047, 4.81228, 5.73703, 6.29845, 5.9948, 6.09279, 6.26104};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(9,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);
  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825);  
cout << "Q2 = 0.875; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[8] =  0.875;
res1415[8] = resvl->GetBinContent(2);
res1516[8] = resvl->GetBinContent(3);
cout << "Q2 = 0.875; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
 resvl->Draw("][");
};



if (qq2==10) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};


//Double_t data_int_mod_res[8] ={1.53701, 2.19138, 3.03188, 4.56564, 5.41792, 5.90783, 5.57599, 5.55767};

Double_t data_int_mod_res[8] ={1.53701, 2.19138, 3.03188, 4.56564, 5.41792, 5.90783, 5.57599, 5.55767};


TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(8,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);
  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825);  
cout << "Q2 = 0.925; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl;
q2gr[9] =  0.925;
res1415[9] = resvl->GetBinContent(2);
res1516[9] = resvl->GetBinContent(3);
cout << "Q2 = 0.925; W = 1.5-1.6; ratio = " <<  resvl->GetBinContent(3) << endl;
 resvl->Draw("][");
};



if (qq2==11) {
Double_t w_values[17]={1.4125,1.4375,1.4625,1.4875,1.5125,1.5375,1.5625,1.5875,1.6125,1.6375,1.6625,1.6875,1.7125,1.7375,1.7625,1.7875,1.8125};


//Double_t data_int_mod_res[6] ={1.52924, 2.13774, 2.92806, 4.35458, 5.142, 5.57416};

Double_t data_int_mod_res[6] ={1.52924, 2.13774, 2.92806, 4.35458, 5.142, 5.57416};



TGraph *gr_int_mod_res;
gr_int_mod_res = new TGraph(6,w_values ,data_int_mod_res);
gr_int_mod_res->SetMarkerStyle(20);
gr_int_mod_res->SetLineWidth(2);
gr_int_mod_res->SetLineStyle(2);
gr_int_mod_res->SetMarkerColor(kBlack);

TH1F *resvl = (TH1F*)h_w_model[qq2]->Clone();

  Double_t x,y; 
  for (Int_t iiii = 0; iiii < gr_int_mod_res->GetN(); iiii++){     gr_int_mod_res->GetPoint(iiii, x, y); 

      y = y * factor;
      gr_int_mod_res->SetPoint(iiii, x, y);
      resvl->SetBinContent(iiii+5,y);
  } 
  
  
  
  
  

//gr_int_mod_res->Draw("same");
resvl->Divide(h_w_model[qq2]);
resvl->SetMinimum(0.);
resvl->SetMaximum(0.8);
  resvl->Rebin(4);
  resvl->Scale(0.25);
resvl->GetXaxis()->SetRangeUser(1.41,1.825); 
cout << "Q2 = 0.975; W = 1.4-1.5; ratio = " <<  resvl->GetBinContent(2) << endl; 
q2gr[10] =  0.975;
res1415[10] = resvl->GetBinContent(2);
 resvl->Draw("][");
};











};



TGraphErrors *gr_res1415;
Double_t res1415_err[11];
for(Int_t i = 0; i < 11; i++) res1415_err[i] = 0.1*res1415[i];
gr_res1415 = new TGraphErrors(11,q2gr,res1415,0,res1415_err);

TCanvas *c_res1415 = new TCanvas("c_res1415","c_res1415",400,300);

c_res1415->cd();
c_res1415->SetBottomMargin(0.15);
c_res1415->SetTopMargin(0.13);
c_res1415->SetLeftMargin(0.15);
c_res1415->SetRightMargin(0.05);
gr_res1415->SetTitle("");
gr_res1415->SetMarkerStyle(20);
for (Int_t i=0;i<gr_res1415->GetN();i++) gr_res1415->GetY()[i] *= 100.;
gr_res1415->GetYaxis()->SetRangeUser(0.,80.);
gr_res1415->GetXaxis()->SetLimits(0.35,1.05);
gr_res1415->GetXaxis()->SetLabelSize(0.07);
gr_res1415->GetYaxis()->SetLabelSize(0.07);
gr_res1415->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
gr_res1415->GetYaxis()->SetTitle("Resonant contribution (%)");
gr_res1415->GetXaxis()->SetTitleSize(0.08);
gr_res1415->GetXaxis()->SetTitleOffset(0.85);
gr_res1415->GetYaxis()->SetTitleSize(0.08);
gr_res1415->GetYaxis()->SetTitleOffset(0.85);
gr_res1415->Draw("APL");



TGraphErrors *gr_res1516;
Double_t res1516_err[10];
for(Int_t i = 0; i < 10; i++) res1516_err[i] = 0.1*res1516[i];
gr_res1516 = new TGraphErrors(10,q2gr,res1516,0,res1516_err);
for (Int_t i=0;i<gr_res1516->GetN();i++) gr_res1516->GetY()[i] *= 100.;
gr_res1516->SetMarkerStyle(21);
gr_res1516->SetMarkerSize(0.8);
//gr_res1516->SetMarkerColor(kBlue);
//gr_res1516->SetLineColor(kBlue);
gr_res1516->Draw("PL");


TGraphErrors *gr_res1617;
Double_t res1617_err[8];
for(Int_t i = 0; i < 8; i++) res1617_err[i] = 0.1*res1617[i];
gr_res1617 = new TGraphErrors(8,q2gr1617,res1617,0,res1617_err);
for (Int_t i=0;i<gr_res1617->GetN();i++) gr_res1617->GetY()[i] *= 100.;
gr_res1617->SetMarkerStyle(22);
gr_res1617->SetMarkerSize(1.1);
//gr_res1617->SetMarkerColor(kRed);
//gr_res1617->SetLineColor(kRed);
gr_res1617->Draw("PL");


TGraphErrors *gr_res1718;
Double_t res1718_err[8];
for(Int_t i = 0; i < 4; i++) res1718_err[i] = 0.1*res1718[i];
gr_res1718 = new TGraphErrors(4,q2gr1617,res1718,0,res1718_err);
for (Int_t i=0;i<gr_res1718->GetN();i++) gr_res1718->GetY()[i] *= 100.;
gr_res1718->SetMarkerStyle(34);
gr_res1718->SetMarkerSize(1.1);
//gr_res1718->SetMarkerColor(kGreen);
//gr_res1718->SetLineColor(kGreen);
gr_res1718->Draw("PL");

leg = new TLegend(0.62,0.73,0.99,0.99);
leg->AddEntry(gr_res1718,"1.7 < W < 1.8 GeV","p");
leg->AddEntry(gr_res1617,"1.6 < W < 1.7 GeV","p");
leg->AddEntry(gr_res1516,"1.5 < W < 1.6 GeV","p");
leg->AddEntry(gr_res1415,"1.4 < W < 1.5 GeV","p");

//leg->AddEntry(gr_int_mod,"JM model, total cr","p");
//leg->AddEntry(gr_int_mod_res,"JM model, res contr","p");
leg->SetTextSize(0.05);
leg->SetFillColor(kWhite);
leg->Draw();


}; //end of main program
