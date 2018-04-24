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







 
 
 
 
void statistics() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
 


Int_t k,i,j;





Double_t w_err;
Float_t Int[21];


Int_t sum_pim = 0;
Int_t sum_pip = 0;
Int_t sum_pr = 0;
Int_t sum_0 = 0;






//Define input files


TFile *file_data = new TFile("data_new_time_corr.root","READ");


 for (Int_t qq2=0; qq2<12;qq2++) {
 Q2_bin = 0.425 + 0.05*qq2;




for (i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {

 W_bin[i] = 1.3125+0.025*i; 


//cout << i << "\n";

TH1D *h_prj_crs;

//Reading input files

//read_data(file_data,i);


file_data->cd();
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_1[i]);
qqq.str("");

cout << Q2_bin << " " << W_bin[i] << " " << " " << h_data_1[i]->GetEntries() <<"\n";

sum_pim = sum_pim + h_data_1[i]->GetEntries();


qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_2[i]);
qqq.str("");

cout << Q2_bin << " " << W_bin[i] << " " << " " << h_data_2[i]->GetEntries() <<"\n";

sum_pip = sum_pip + h_data_2[i]->GetEntries();


qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_3[i]);
qqq.str("");

cout << Q2_bin << " " << W_bin[i] << " " << " " << h_data_3[i]->GetEntries() <<"\n";

sum_pr = sum_pr + h_data_3[i]->GetEntries();

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_4[i]);
qqq.str("");

cout << Q2_bin << " " << W_bin[i] << " " << " " << h_data_4[i]->GetEntries() <<"\n";

sum_0 = sum_0 + h_data_4[i]->GetEntries();


};
};


cout << "sum_pim = " << sum_pim << endl;
cout << "sum_pip = " << sum_pip << endl;
cout << "sum_pr = " << sum_pr << endl;
cout << "sum_0 = " << sum_0 << endl;
cout << "total = " << sum_0 + sum_pr + sum_pip + sum_pim<< endl;

}; //end of main program


