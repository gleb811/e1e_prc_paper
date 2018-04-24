void theta_var_top() {

gStyle->SetOptStat(0);

THnSparseD *h_pim[3];
THnSparseD *h_pip[3];
THnSparseD *h_pr[3];
THnSparseD *h_0[3];

TH1D *h_th_pim[4];
TH1D *h_th_pr[4];
TH1D *h_th_pip[4];

TCanvas *c = new TCanvas("c","c",900,300);
c->Divide(3,1);

ostringstream qqq;
TFile *file_data = new TFile("data_new_time_corr.root","READ");

//The first set of variables

file_data->cd();
qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pim_1_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pim[0]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pip_1_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pip[0]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pr_1_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pr[0]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_0_1_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_0[0]);
qqq.str("");


c->cd(1);

h_th_pr[0]=h_pim[0]->Projection(2,"");
h_th_pr[1]=h_pip[0]->Projection(2,"");
h_th_pr[2]=h_pr[0]->Projection(2,"");
h_th_pr[3]=h_0[0]->Projection(2,"");

c->cd(1)->SetLeftMargin(0.17);
c->cd(1)->SetTopMargin(0.05);
c->cd(1)->SetRightMargin(0.05);
c->cd(1)->SetBottomMargin(0.17);

h_th_pr[0]->SetTitle("");
h_th_pr[0]->GetXaxis()->SetTitle("#theta_{p}^{cms} (deg)");
h_th_pr[0]->GetXaxis()->SetTitleSize(0.07);
h_th_pr[0]->GetXaxis()->SetTitleOffset(1.1);
h_th_pr[0]->GetXaxis()->SetNdivisions(-506);
h_th_pr[0]->GetXaxis()->SetLabelSize(0.06);
h_th_pr[0]->GetYaxis()->SetTitle("number of events");
h_th_pr[0]->GetYaxis()->SetTitleSize(0.07);
h_th_pr[0]->GetYaxis()->SetTitleOffset(1.3);
h_th_pr[0]->GetYaxis()->SetLabelSize(0.06);
h_th_pr[0]->SetAxisRange(0.,1200.,"Y");

TLegend *leg_pr = new TLegend(0.2,0.45,0.45,0.99);


h_th_pr[0]->SetLineColor(1);
h_th_pr[0]->SetLineWidth(2);
h_th_pr[0]->Draw();
leg_pr->AddEntry(h_th_pr[0],"#pi^{-} miss","l");

h_th_pr[1]->SetLineColor(2);
h_th_pr[1]->SetLineWidth(2);
h_th_pr[1]->Draw("same");
leg_pr->AddEntry(h_th_pr[1],"#pi^{+} miss","l");

h_th_pr[2]->SetLineColor(3);
h_th_pr[2]->SetLineWidth(2);
h_th_pr[2]->Draw("same");
leg_pr->AddEntry(h_th_pr[2],"p miss","l");

h_th_pr[3]->SetLineColor(4);
h_th_pr[3]->SetLineWidth(2);
h_th_pr[3]->Draw("same");
leg_pr->AddEntry(h_th_pr[3],"exclusive","l");


leg_pr->Draw();






//The second set of variables

file_data->cd();
qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pim_2_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pim[1]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pip_2_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pip[1]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pr_2_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pr[1]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_0_2_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_0[1]);
qqq.str("");


c->cd(2);

h_th_pim[0]=h_pim[1]->Projection(2,"");
h_th_pim[1]=h_pip[1]->Projection(2,"");
h_th_pim[2]=h_pr[1]->Projection(2,"");
h_th_pim[3]=h_0[1]->Projection(2,"");

TLegend *leg_pim = new TLegend(0.67,0.45,0.92,0.99);

c->cd(2)->SetLeftMargin(0.17);
c->cd(2)->SetTopMargin(0.05);
c->cd(2)->SetRightMargin(0.05);
c->cd(2)->SetBottomMargin(0.17);

h_th_pim[0]->SetTitle("");
h_th_pim[0]->GetXaxis()->SetTitle("#theta_{#pi^{-}}^{cms} (deg)");
h_th_pim[0]->GetXaxis()->SetTitleSize(0.07);
h_th_pim[0]->GetXaxis()->SetTitleOffset(1.1);
h_th_pim[0]->GetXaxis()->SetNdivisions(-506);
h_th_pim[0]->GetXaxis()->SetLabelSize(0.06);
h_th_pim[0]->GetYaxis()->SetTitle("number of events");
h_th_pim[0]->GetYaxis()->SetTitleSize(0.07);
h_th_pim[0]->GetYaxis()->SetTitleOffset(1.3);
h_th_pim[0]->GetYaxis()->SetLabelSize(0.06);
h_th_pim[0]->SetAxisRange(0.,1200.,"Y");

h_th_pim[0]->SetLineColor(1);
h_th_pim[0]->SetLineWidth(2);
h_th_pim[0]->Draw();
leg_pim->AddEntry(h_th_pim[0],"#pi^{-} miss","l");

h_th_pim[1]->SetLineColor(2);
h_th_pim[1]->SetLineWidth(2);
h_th_pim[1]->Draw("same");
leg_pim->AddEntry(h_th_pim[1],"#pi^{+} miss","l");

h_th_pim[2]->SetLineColor(3);
h_th_pim[2]->SetLineWidth(2);
h_th_pim[2]->Draw("same");
leg_pim->AddEntry(h_th_pim[2],"p miss","l");

h_th_pim[3]->SetLineColor(4);
h_th_pim[3]->SetLineWidth(2);
h_th_pim[3]->Draw("same");
leg_pim->AddEntry(h_th_pim[3],"exclusive","l");


leg_pim->Draw();


//The third set of variables

file_data->cd();
qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pim_3_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pim[2]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pip_3_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pip[2]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_pr_3_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_pr[2]);
qqq.str("");

qqq << "q2_0.475" << "/w_1.6125"  << "/h_5dim_0_3_q2_475" << "_w_1612.5";
gDirectory->GetObject(qqq.str().c_str(),h_0[2]);
qqq.str("");


c->cd(3);

h_th_pip[0]=h_pim[2]->Projection(2,"");
h_th_pip[1]=h_pip[2]->Projection(2,"");
h_th_pip[2]=h_pr[2]->Projection(2,"");
h_th_pip[3]=h_0[2]->Projection(2,"");

TLegend *leg_pip = new TLegend(0.76,0.45,0.99,0.99);


c->cd(3)->SetLeftMargin(0.17);
c->cd(3)->SetTopMargin(0.05);
c->cd(3)->SetRightMargin(0.05);
c->cd(3)->SetBottomMargin(0.17);

h_th_pip[0]->SetTitle("");
h_th_pip[0]->GetXaxis()->SetTitle("#theta_{#pi^{+}}^{cms} (deg)");
h_th_pip[0]->GetXaxis()->SetTitleSize(0.07);
h_th_pip[0]->GetXaxis()->SetTitleOffset(1.1);
h_th_pip[0]->GetXaxis()->SetNdivisions(-506);
h_th_pip[0]->GetXaxis()->SetLabelSize(0.06);
h_th_pip[0]->GetYaxis()->SetTitle("number of events");
h_th_pip[0]->GetYaxis()->SetTitleSize(0.07);
h_th_pip[0]->GetYaxis()->SetTitleOffset(1.3);
h_th_pip[0]->GetYaxis()->SetLabelSize(0.06);
h_th_pip[0]->SetAxisRange(0.,1200.,"Y");

h_th_pip[0]->SetLineColor(1);
h_th_pip[0]->SetLineWidth(2);
h_th_pip[0]->Draw();
leg_pip->AddEntry(h_th_pip[0],"#pi^{-} miss","l");

h_th_pip[1]->SetLineColor(2);
h_th_pip[1]->SetLineWidth(2);
h_th_pip[1]->Draw("same");
leg_pip->AddEntry(h_th_pip[1],"#pi^{+} miss","l");

h_th_pip[2]->SetLineColor(3);
h_th_pip[2]->SetLineWidth(2);
h_th_pip[2]->Draw("same");
leg_pip->AddEntry(h_th_pip[2],"p miss","l");

h_th_pip[3]->SetLineColor(4);
h_th_pip[3]->SetLineWidth(2);
h_th_pip[3]->Draw("same");
leg_pip->AddEntry(h_th_pip[3],"exclusive","l");

leg_pip->Draw();
}
