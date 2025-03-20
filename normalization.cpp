#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH2.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TBenchmark.h>
#include <TSystemDirectory.h>
#include <stdlib.h>
#include <TLorentzVector.h>
#include <clas12reader.h>
#include <HipoChain.h>
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include <TStyle.h>

using namespace clas12;


// 假设 hist 是你已经填充完毕的 TH1D 对象
void ScaleHistogramY(TH1* hist, double scale_factor) {
    // 遍历所有的 bin（包括 underflow 和 overflow）
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double content = hist->GetBinContent(i);        // 获取 bin 的内容
        double error = hist->GetBinError(i);            // 获取 bin 的误差

        // 乘以系数
        hist->SetBinContent(i, content * scale_factor);
        hist->SetBinError(i, error * scale_factor);     // 更新误差值（假设误差也需要被缩放）
    }
}

void normalization(){
    TFile * f = TFile::Open("/work/clas12/yijie/02DVCS_Analysis/auto_bin_analysis/hist_E6535_in_data_corrected_W_cut_apply_thetaHR.root","read");
    //TFile * f = TFile::Open("dst6535.root","read");

	TTree * t = (TTree*) f->Get("trackTree");
    TFile *f_ratio = TFile::Open("/work/clas12/yijie/02DVCS_Analysis/auto_bin_analysis/Acceptance_Rate_RECON.root","read");
    TH1* h_ratio = dynamic_cast<TH1*>(f_ratio->Get("h_ratio"));

    //TFile * histFile = TFile::Open("hist_E7546.root","recreate");
    TFile * histFile = TFile::Open("hist_E6535_in_cross_section_s1_correctedHR.root","recreate");


    TH1* h_unnormalized = dynamic_cast<TH1*>(f->Get("Hist1D_Electron_theta_selected_s1theta5.000000to20.000000"));
    h_unnormalized->SetTitle("h_unnormalized");

    TH1* h_normalized = (TH1*)h_unnormalized->Clone("h_normalized");
    h_normalized->SetTitle("h_normalized");

    


    //Beam Charge
    //float charge_accumlated=133.0*138857.0;//7.5GeV
    //float charge_accumlated=104*259504/(6.6);//104*259504*(255248280/407714413)*(12524/16437);//6.5GeV
    float charge_accumlated=18.480766598083676*pow(10,6)/6.66;//actually should be 17.645871262811298mc
    //float charge_sample=0.879*104*259504*0.477009758418*2379277/255248280;
    cout<<charge_accumlated<<" nC for RGK Fall2018 Pass1 6.5GeV"<<endl;

    float NA=6.022*pow(10,23);//1
    float e_charge=1.602*pow(10,-19);//C
    float target_length=5;//cm
    float target_density=0.07;//g/cm^3
    cout<<"N_A = "<<NA<<endl;
    cout<<"e_charge/C = "<<e_charge<<endl;
    cout<<"target_length/cm = "<<target_length<<endl;
    cout<<"target_density/g/cm^3 = "<<target_density<<endl;

    float luminosity_accumlated=(NA*target_length*target_density*charge_accumlated*pow(10,-9))*pow(10,-40)/e_charge;
    cout<<luminosity_accumlated<<"*10^(40)  cm^(-2) luminosity for RGK Fall2018 Pass1 6.5GeV"<<endl;

    //float luminosity_sample=(NA*target_length*target_density*charge_sample*pow(10,-9))*pow(10,-40)/e_charge;
    //cout<<luminosity_sample<<"*10^(40)  cm^(-2) luminosity for sample"<<endl;


    //float luminosity=luminosity_accumlated*0.00109721342766*1.3*870/94;
    float luminosity=luminosity_accumlated;
    cout<<"Use "<<luminosity<<"*10^(40)  cm^(-2) for this program..."<<endl;


    float PI=3.14159265359;
    float RadtoDeg=180.0/PI;

    float theta_begin=5.0*PI/180;
    float theta_end=20.0*PI/180;
    float theta_bin=75;
    float deltatheta=(theta_end-theta_begin)/theta_bin;

    float begin_phi=-23;
    float end_phi=337;
    //float phi_bins=6;
    //float deltaphibins=(end_phi-begin_phi)/phi_bins;
    float deltaphibins=22;

    float eff=0.8;//useless
    float selected_phi_ratio=deltaphibins/360;
    cout<<"Phi selected Ratio: "<<selected_phi_ratio*100<<"%"<<endl;
    float normalize_factor=pow(10,33)*1.0/(luminosity*pow(10,40)*deltatheta*eff*selected_phi_ratio);
    cout<<"Overall Reconstruction acceptance = "<<eff*100<<"%"<<endl;
    cout<<"Normaliaze factor/nb = "<<normalize_factor<<endl;
    //float event_calculated=0;

    
    //h_normalized->Divide(h_ratio);
    ScaleHistogramY(h_normalized,normalize_factor);
    h_normalized->Write();
    h_unnormalized->Write();
    h_ratio->Write();

    gApplication->Terminate(0);
}