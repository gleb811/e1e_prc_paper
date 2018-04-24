//TH1D *h_cos_th;
# define M_PI           3.14159265358979323846  /* pi */

THnSparseD *h_data_1[21],*h_data_2[21],*h_data_3[21],*h_data_4[21];
THnSparseD *h_data_1_tmp[21],*h_data_2_tmp[21],*h_data_3_tmp[21];
THnSparseD *h_empty_1[21],*h_empty_2[21],*h_empty_3[21];
THnSparseD *h_empty_fill_1[21],*h_empty_fill_2[21],*h_empty_fill_3[21];
THnSparseD *h_rec_1[21],*h_rec_2[21],*h_rec_3[21];

THnSparseD *h_gen_1[21],*h_gen_2[21],*h_gen_3[21];
THnSparseD *h_eff_1[21],*h_eff_2[21],*h_eff_3[21];
THnSparseD *h_cr_sect_1[21],*h_cr_sect_2[21],*h_cr_sect_3[21];


THnSparseD *h_cr_sect_noemptcells_1[21],*h_cr_sect_noemptcells_2[21],*h_cr_sect_noemptcells_3[21];
THnSparseD *h_rec_unit_1[21],*h_rec_unit_2[21],*h_rec_unit_3[21];
THnSparseD *h_gen_unit_1[21],*h_gen_unit_2[21],*h_gen_unit_3[21];
THnSparseD *model_1[21],*model_2[21],*model_3[21];

Float_t rad_corr[25];
ostringstream qqq;
 Float_t W_bin[21];
 Float_t Q2_bin = 0.475;
Float_t m_proton = 0.938272;
 TLegend *leg;
 TLegend *leg_w_int;
 
 
 

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
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 8;
return get_min_w;
};







 
 
 
 
void avrg_eff() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
 
gErrorIgnoreLevel = kError;

Int_t k,i,j;





Double_t w_err;
Float_t Int[21];


Int_t sum_pim = 0;
Int_t sum_pip = 0;
Int_t sum_pr = 0;
Int_t sum_0 = 0;





//Define input files


TFile *file_data = new TFile("combine_small_bin_all_top.root","READ");


Float_t tmp_eff = 0.;
Float_t num_of_cells = 0.;
Float_t sum_eff = 0.;

 for (Int_t qq2=0; qq2<12;qq2++) {
 Q2_bin = 0.425 + 0.05*qq2;




for (i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {

 W_bin[i] = 1.3125+0.025*i; 


//cout << i << "\n";

TH1D *h_prj_crs;

//Reading input files

//read_data(file_data,i);


file_data->cd();
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_1[i]);
qqq.str("");

//cout << Q2_bin << " " << W_bin[i] << " " << " " << h_rec_1[i]->GetEntries() <<"\n";

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_1[i]);
qqq.str("");

//cout << Q2_bin << " " << W_bin[i] << " " << " " << h_gen_1[i]->GetEntries() <<"\n";


h_rec_1[i]->Divide(h_gen_1[i]);

Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t y_max = 8;
Int_t t_max = 8;

if ((i==0)||(i==1)) {
o_max = p_max = 8;
r_max = 6;
y_max = 5; 
t_max = 5;
};
if ((i==2)||(i==3)) {
o_max = p_max = 10;
r_max = 8;
y_max = 6;
t_max = 5;
}; 
if ((i>=4)&&(i<=6)){
o_max = p_max =12;
r_max =10;
y_max = 8;
t_max = 5;
};
if ((i>=7)&&(i<=14)){
o_max = p_max =12;
r_max =10;
t_max = 8;
y_max = 8;
};



for (Int_t o=1; o<=o_max; o++) {
for (Int_t p=1; p<=p_max; p++) {
for (Int_t r=1; r<=r_max; r++) {
for (Int_t t=1; t<=t_max; t++) {
for (Int_t y=1; y<=y_max; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;


Long64_t tmp_bin = h_rec_1[i]->GetBin(bins);
tmp_eff = h_rec_1[i]->GetBinContent(tmp_bin);

if ((tmp_eff > 0.)&&(tmp_eff < 0.3)) {
//cout << tmp_eff << endl;;
sum_eff = sum_eff + tmp_eff;
num_of_cells = num_of_cells + 1.;
};


};
};
};
}; 
};





//cout << Q2_bin << " " << W_bin[i] << " " <<  sum_eff/num_of_cells <<"\n";
};
};
cout << sum_eff/num_of_cells <<"\n";

//cout << "sum_pim = " << sum_pim << endl;


}; //end of main program


