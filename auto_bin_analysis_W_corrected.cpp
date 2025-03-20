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

double bifurGauss(double* x, double* par) {
    // x[0]: 自变量
    // par[0]: mean (均值)
    // par[1]: sigma_L (左侧标准差)
    // par[2]: sigma_R (右侧标准差)
    double mean = par[0];
    double sigma_L = par[1];
    double sigma_R = par[2];
    double mag = par[3];

    if (x[0] < mean) {
        return mag *
               TMath::Exp(-0.5 * TMath::Power((x[0] - mean) / sigma_L, 2));
    } else {
        return mag *
               TMath::Exp(-0.5 * TMath::Power((x[0] - mean) / sigma_R, 2));
    }
}

// 常量定义
const double E1 = 6.535;
const double M = 0.938;



bool isInRange(float num) {
    // 检查 num 是否在指定的范围内
    return (num >= -7.0f && num <= 15.0f) ||
           (num >= 53.0f && num <= 75.0f) ||
           (num >= 113.0f && num <= 135.0f) ||
           (num >= 173.0f && num <= 195.0f) ||
           (num >= 233.0f && num <= 255.0f) ||
           (num >= 293.0f && num <= 315.0f);
}

double theta_deltaP(int sector, double x) {
    switch (sector) {
        case 1:
            return 0.0000038941 * pow(x, 4) - 0.0002038678 * pow(x, 3) + 0.0035549624 * pow(x, 2) - 0.0214692377 * x + 0.0303117187;
        case 2:
            return 0.0000033809 * pow(x, 4) - 0.0001221358 * pow(x, 3) + 0.0009052173 * pow(x, 2) + 0.0068831943 * x - 0.0714679562;
        case 3:
            return 0.0000054717 * pow(x, 4) - 0.0002912588 * pow(x, 3) + 0.0052684411 * pow(x, 2) - 0.0371702272 * x + 0.0758673367;
        case 4:
            return 0.0000052404 * pow(x, 4) - 0.0003044749 * pow(x, 3) + 0.0059914364 * pow(x, 2) - 0.0463647967 * x + 0.1209073661;
        case 5:
            return 0.0000078980 * pow(x, 4) - 0.0004030582 * pow(x, 3) + 0.0070267467 * pow(x, 2) - 0.0483415613 * x + 0.0997526245;
        case 6:
            return 0.0000037816 * pow(x, 4) - 0.0001674624 * pow(x, 3) + 0.0019722212 * pow(x, 2) - 0.0003911925 * x - 0.0511021970;
        default:
            std::cerr << "Invalid sector: " << sector << ". Sector must be between 1 and 6.\n";
            return NAN;  // 返回 NaN 表示无效输入
    }
}

int get_sector(double x) {
    if (x >= -23 && x < 37) return 1;
    if (x >= 37 && x < 97) return 2;
    if (x >= 97 && x < 157) return 3;
    if (x >= 157 && x < 217) return 4;
    if (x >= 217 && x < 277) return 5;
    if (x >= 277 && x <= 337) return 6;
    
    return -1;  // 返回 -1 表示无效的输入
}



void auto_bin_analysis_W_corrected(){
    //TFile * f = TFile::Open("/work/clas12/yijie/02DVCS_Analysis/auto_bin_analysis/dst6535_in_MC_rate.root","read");
    TFile * f = TFile::Open("/work/clas12/yijie/02DVCS_Analysis/Fiducal_cut/job_out/dst6535in.root","read");
	TTree * t = (TTree*) f->Get("trackTree");

    TFile * histFile = TFile::Open("hist_E6535_in_data_highR_all_theta.root","recreate");

    float beam_energy=6.535; //GeV
	float Proton_mass=0.93827208816; //GeV
	float Electron_mass=0.00051099895000; //Gev
    float PI=3.14159265359;
    float RadtoDeg=180.0/PI;

    float theta_begin=5.0*PI/180;
    float theta_end=20.0*PI/180;
    float theta_bin=1;
    float deltatheta=(theta_end-theta_begin)/theta_bin;

    float begin_phi=-23;
    float end_phi=337;
    float phi_bins=6;


    TH1 * Hist1D_Electron_theta = new TH1D("Electron_theta","Electron_theta",theta_bin,theta_begin*RadtoDeg,theta_end*RadtoDeg);
    TH1 * Hist1D_Electron_phi = new TH1D("Electron_phi","Electron_phi",500,-30,360);
    TH1 * Hist1D_Electron_E = new TH1D("Electron_E","Electron_E",500,0,8);
    TH1 * Hist1D_Electron_sector = new TH1D("Electron_sector","Electron_sector",6,1,7);
    TH1 * Hist1D_Electron_vx = new TH1D("Electron_vx","Electron_vx",500,-20,20);
    TH1 * Hist1D_Electron_vy = new TH1D("Electron_vy","Electron_vy",500,-20,20);
    TH1 * Hist1D_Electron_vz = new TH1D("Electron_vz","Electron_vz",500,-20,20);

    

    TH1 * Hist1D_Q2 = new TH1D("Q2","Q2",500,0,6);
    TH1 * Hist1D_W = new TH1D("W","W",500,0.0,2.5);

    //TH1 * Hist1D_Electron_E_corrected_selected = new TH1D("Electron_E_corrected_selected","Electron_E_corrected_selected",500,0,8);
    //TH1 * Hist1D_Q2_corrected_selected = new TH1D("Q2_corrected_selected","Q2_corrected_selected",500,0,6);
    //TH1 * Hist1D_W_corrected_selected = new TH1D("W_corrected_selected","W_corrected_selected",500,0,2.5);
    TH2 * Hist2D_W_theta_S1_selected = new TH2D("W_theta_S1_selected","W_theta_S1_selected",500,5,20,500,-0.3,0.3);
    TH2 * Hist2D_W_theta_S1_corrected_selected = new TH2D("W_theta_S1_corrected_selected","W_theta_S1_corrected_selected",500,5,20,500,-0.3,0.3);

    TH2 * Hist2D_W_theta_S2_selected = new TH2D("W_theta_S2_selected","W_theta_S2_selected",500,5,20,500,-0.3,0.3);
    TH2 * Hist2D_W_theta_S2_corrected_selected = new TH2D("W_theta_S2_corrected_selected","W_theta_S2_corrected_selected",500,5,20,500,-0.3,0.3);

    TH2 * Hist2D_W_theta_S3_selected = new TH2D("W_theta_S3_selected","W_theta_S3_selected",500,5,20,500,-0.3,0.3);
    TH2 * Hist2D_W_theta_S3_corrected_selected = new TH2D("W_theta_S3_corrected_selected","W_theta_S3_corrected_selected",500,5,20,500,-0.3,0.3);

    TH2 * Hist2D_W_theta_S4_selected = new TH2D("W_theta_S4_selected","W_theta_S4_selected",500,5,20,500,-0.3,0.3);
    TH2 * Hist2D_W_theta_S4_corrected_selected = new TH2D("W_theta_S4_corrected_selected","W_theta_S4_corrected_selected",500,5,20,500,-0.3,0.3);

    TH2 * Hist2D_W_theta_S5_selected = new TH2D("W_theta_S5_selected","W_theta_S5_selected",500,5,20,500,-0.3,0.3);
    TH2 * Hist2D_W_theta_S5_corrected_selected = new TH2D("W_theta_S5_corrected_selected","W_theta_S5_corrected_selected",500,5,20,500,-0.3,0.3);

    TH2 * Hist2D_W_theta_S6_selected = new TH2D("W_theta_S6_selected","W_theta_S6_selected",500,5,20,500,-0.3,0.3);
    TH2 * Hist2D_W_theta_S6_corrected_selected = new TH2D("W_theta_S6_corrected_selected","W_theta_S6_corrected_selected",500,5,20,500,-0.3,0.3);

    TLorentzVector p4_electron_corrected;
	TLorentzVector beam;
	TLorentzVector target;
    beam.SetXYZM(0,0,beam_energy,Electron_mass);
	target.SetXYZM(0,0,0,Proton_mass);

    float Q2_readout=0;
    float W_readout=0;
    int Electron_sector_readout=0;
    std::vector<float> Q2_store;
    std::vector<float> W_store;
    std::vector<int> Electron_sector_store;
    t->SetBranchAddress("Q2",&Q2_readout);
    t->SetBranchAddress("W",&W_readout);
    t->SetBranchAddress("Electron_sector",&Electron_sector_readout);

    std::vector<float> *Electron_theta_readout=0;
    std::vector<float> Electron_theta_store;
    std::vector<float> *Electron_phi_readout=0;
    std::vector<float> Electron_phi_store;
    std::vector<float> *Electron_E_readout=0;
    std::vector<float> Electron_E_store;
    std::vector<float> *Electron_vx_readout=0;
    std::vector<float> Electron_vx_store;
    std::vector<float> *Electron_vy_readout=0;
    std::vector<float> Electron_vy_store;
    std::vector<float> *Electron_vz_readout=0;
    std::vector<float> Electron_vz_store;

    std::vector<float> Electron_E_corrected_store;
    std::vector<float> Q2_corrected_store;
    std::vector<float> W_corrected_store;
    
    t->SetBranchAddress("Electron_theta",&Electron_theta_readout); 
    t->SetBranchAddress("Electron_phi",&Electron_phi_readout);
    t->SetBranchAddress("Electron_E",&Electron_E_readout);
    t->SetBranchAddress("Electron_vx",&Electron_vx_readout);
    t->SetBranchAddress("Electron_vy",&Electron_vy_readout);
    t->SetBranchAddress("Electron_vz",&Electron_vz_readout);

    for (int i = 0; i < t->GetEntries() ; i++){
    //for (int i = 0; i < 100000 ; i++){
        t->GetEntry(i);
        Hist1D_Q2->Fill(Q2_readout);
        Q2_store.push_back(Q2_readout);
        Hist1D_W->Fill(W_readout);
        W_store.push_back(W_readout);
        Electron_sector_store.push_back(Electron_sector_readout);
        Hist1D_Electron_sector->Fill(Electron_sector_readout);
        for (int j = 0; j < Electron_theta_readout->size(); j++){
            Electron_theta_store.push_back(Electron_theta_readout->at(j));
            Hist1D_Electron_theta->Fill(Electron_theta_readout->at(j)*RadtoDeg);
            Electron_E_store.push_back(Electron_E_readout->at(j));
            Hist1D_Electron_E->Fill(Electron_E_readout->at(j));
            float temp_phi=0;
            if(Electron_phi_readout->at(j)>=-23*PI/180){
                Electron_phi_store.push_back(Electron_phi_readout->at(j));
                temp_phi=Electron_phi_readout->at(j);
                Hist1D_Electron_phi->Fill(Electron_phi_readout->at(j)*RadtoDeg);
            }
            else{
                Electron_phi_store.push_back(Electron_phi_readout->at(j)+2*PI);
                temp_phi=Electron_phi_readout->at(j)+2*PI;
                Hist1D_Electron_phi->Fill(Electron_phi_readout->at(j)*RadtoDeg+360);
            }
            
            Electron_vx_store.push_back(Electron_vx_readout->at(j));
            Hist1D_Electron_vx->Fill(Electron_vx_readout->at(j));
            Electron_vy_store.push_back(Electron_vy_readout->at(j));
            Hist1D_Electron_vy->Fill(Electron_vy_readout->at(j));
            Electron_vz_store.push_back(Electron_vz_readout->at(j));
            Hist1D_Electron_vz->Fill(Electron_vz_readout->at(j));
        
            float temp_corrected_E=theta_deltaP(get_sector(temp_phi*RadtoDeg),Electron_theta_readout->at(j)*RadtoDeg)+Electron_E_readout->at(j);
            //cout<<Electron_E_readout->at(j)<<","<<temp_corrected_E<<endl;
            float temp_Pz=cos(Electron_theta_readout->at(j))*temp_corrected_E;
            float temp_Px=sin(Electron_theta_readout->at(j))*temp_corrected_E;
            p4_electron_corrected.SetXYZM(temp_Px,0,temp_Pz,Electron_mass);
            float temp_Q2=-(beam-p4_electron_corrected).M2();
		    float temp_W=TMath::Sqrt((beam+target-p4_electron_corrected).M2());
            //cout<<W_readout<<","<<temp_W<<endl;
            Q2_corrected_store.push_back(temp_Q2);
            W_corrected_store.push_back(temp_W);
            Electron_E_corrected_store.push_back(temp_corrected_E);
        }
    }
    Hist1D_Q2->Write();
    Hist1D_W->Write();

    Hist1D_Electron_theta->Write();
    Hist1D_Electron_phi->Write();
    Hist1D_Electron_E->Write();
    Hist1D_Electron_sector->Write();
    Hist1D_Electron_vx->Write();
    Hist1D_Electron_vy->Write();
    Hist1D_Electron_vz->Write();



    std::map<std::string, TH1*> histMap;
    std::string histName;
    std::string histTitle;

    std::map<std::string, TH1*> histMap_phi;
    std::string histName_phi;
    std::string histTitle_phi;

    std::map<std::string, TH1*> histMap_W;
    std::string histName_W;
    std::string histTitle_W; 

    TF1 *f0 = new TF1("f0","gaus",0.8,1.05);
    

    FILE *fp;
    fp = fopen("theta_W_in_data_highR_all_theta.dat","a");

    for (int l=0; l<theta_bin;l++){
        for (int k = 0; k < 6; k++){
            


            histName = "Hist1D_Electron_theta_selected_s" + std::to_string(k+1)+"theta"+std::to_string((theta_begin+l*deltatheta)*RadtoDeg)+"to"+std::to_string((theta_begin+(l+1)*deltatheta)*RadtoDeg);
            histTitle = "E_theta Sector " + std::to_string(k+1)+"theta"+std::to_string((theta_begin+l*deltatheta)*RadtoDeg)+"to"+std::to_string((theta_begin+(l+1)*deltatheta)*RadtoDeg);
            histMap[histName] = new TH1D(histName.c_str(), histTitle.c_str(),100,theta_begin*RadtoDeg,theta_end*RadtoDeg);
            histMap[histName]->SetStats(0);

            histName_W = "Hist1D_W_selected_s" + std::to_string(k+1)+"theta"+std::to_string((theta_begin+l*deltatheta)*RadtoDeg)+"to"+std::to_string((theta_begin+(l+1)*deltatheta)*RadtoDeg);
            histTitle_W = "Normalized W Sector " + std::to_string(k+1)+"theta"+std::to_string((theta_begin+l*deltatheta)*RadtoDeg)+"to"+std::to_string((theta_begin+(l+1)*deltatheta)*RadtoDeg);
            histMap_W[histName_W] = new TH1D(histName_W.c_str(), histTitle_W.c_str(),100,0.5,1.2);
            histMap_W[histName_W]->SetStats(0);

            histName_phi = "Hist1D_Electron_phi_selected_s" + std::to_string(k+1)+"theta"+std::to_string((theta_begin+l*deltatheta)*RadtoDeg)+"to"+std::to_string((theta_begin+(l+1)*deltatheta)*RadtoDeg);
            histTitle_phi = "E_phi Sector " + std::to_string(k+1)+"theta"+std::to_string((theta_begin+l*deltatheta)*RadtoDeg)+"to"+std::to_string((theta_begin+(l+1)*deltatheta)*RadtoDeg);
            histMap_phi[histName_phi] = new TH1D(histName_phi.c_str(), histTitle_phi.c_str(),1000,-30,370);
            histMap_phi[histName_phi]->SetStats(0);


            for (int i = 0; i < Electron_E_store.size(); i++){
                //if ((Electron_phi_store[i]*RadtoDeg>=begin_phi+k*deltaphibins)&&(Electron_phi_store[i]*RadtoDeg<begin_phi+(k+1)*deltaphibins)){
                if (isInRange(Electron_phi_store[i]*RadtoDeg)&&(Electron_phi_store[i]*RadtoDeg>=begin_phi+k*60)&&(Electron_phi_store[i]*RadtoDeg<begin_phi+(k+1)*60)){
                    if ((Electron_theta_store[i]>=theta_begin+l*deltatheta)&&(Electron_theta_store[i]<=theta_begin+(l+1)*deltatheta)){
                        //histMap[histName]->Fill(Electron_theta_store[i]*RadtoDeg,1/norm_S(Electron_theta_store[i]*RadtoDeg));
                        histMap[histName]->Fill(Electron_theta_store[i]*RadtoDeg);
                        //histMap_W[histName_W]->Fill(W_corrected_store[i]);
                        histMap_W[histName_W]->Fill(W_store[i]);
                        histMap_phi[histName_phi]->Fill(Electron_phi_store[i]*RadtoDeg);
                        //event_calculated++;
                        if(k==0){
                            Hist2D_W_theta_S1_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_store[i]-Proton_mass);
                            Hist2D_W_theta_S1_corrected_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_corrected_store[i]-Proton_mass);
                        }
                        if(k==1){
                            Hist2D_W_theta_S2_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_store[i]-Proton_mass);
                            Hist2D_W_theta_S2_corrected_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_corrected_store[i]-Proton_mass);
                        }
                        if(k==2){
                            Hist2D_W_theta_S3_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_store[i]-Proton_mass);
                            Hist2D_W_theta_S3_corrected_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_corrected_store[i]-Proton_mass);
                        }
                        if(k==3){
                            Hist2D_W_theta_S4_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_store[i]-Proton_mass);
                            Hist2D_W_theta_S4_corrected_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_corrected_store[i]-Proton_mass);
                        }
                        if(k==4){
                            Hist2D_W_theta_S5_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_store[i]-Proton_mass);
                            Hist2D_W_theta_S5_corrected_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_corrected_store[i]-Proton_mass);
                        }
                        if(k==5){
                            Hist2D_W_theta_S6_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_store[i]-Proton_mass);
                            Hist2D_W_theta_S6_corrected_selected->Fill(Electron_theta_store[i]*RadtoDeg,W_corrected_store[i]-Proton_mass);
                        }
                        
                        //cout<<W_store[i]<<","<<W_corrected_store[i]<<endl;
                    }
                }
            
            }

            histMap[histName]->SetXTitle("Electron #theta [deg]");  // 电子的θ角度（单位：度）
            histMap[histName]->SetYTitle("count");
            //histMap[histName]->SetYTitle("#frac{d#sigma}{d#theta} [nb]"); // LaTeX 样式的 dσ/dθ
            //histMap[histName]->SetYTitle("Ratio"); // LaTeX 样式的 dσ/dθ
            //histMap[histName]->SetMinimum(1e-3);  // 设置 y 轴最小值，例如 0.01
            //histMap[histName]->SetMaximum(2000);   // 设置 y 轴最大值，例如 1000

            TCanvas* canvas = new TCanvas((histName + "_canvas").c_str(), histTitle.c_str(), 800, 600);
            canvas->SetLogy();
            histMap[histName]->Draw();
            //TLegend *legend = new TLegend(0.4, 0.75, 0.85, 0.85); // 设置图例的位置，稍微缩小和居中
            //legend->AddEntry(histMap[histName], "RGK Fall2018 Pass2 (E=7.546GeV)", "leq");
            //legend->AddEntry(histMap[histName], "MC::Bank (E=7.546GeV)", "leq");
            //legend->AddEntry(func, "Rosenbluth Formula (E=7.546GeV)", "l");
            // 调整图例的外观
            //legend->SetTextSize(0.03);  // 设置字体大小
            //legend->Draw();
            std::string savePath = "hist_output/" + histName + ".png";
            canvas->SaveAs(savePath.c_str());
            histMap[histName]->Write();

            histMap[histName]->SetXTitle("Electron #theta [deg]");  // 电子的θ角度（单位：度）
            histMap[histName]->SetYTitle("count");
            //histMap[histName]->SetYTitle("#frac{d#sigma}{d#theta} [nb]"); // LaTeX 样式的 dσ/dθ
            //histMap[histName]->SetYTitle("Ratio"); // LaTeX 样式的 dσ/dθ
            //histMap[histName]->SetMinimum(1e-3);  // 设置 y 轴最小值，例如 0.01
            //histMap[histName]->SetMaximum(2000);   // 设置 y 轴最大值，例如 1000

            TCanvas* canvas_phi = new TCanvas((histName_phi + "_canvas").c_str(), histTitle_phi.c_str(), 800, 600);
            //canvas_phi->SetLogy();
            histMap_phi[histName_phi]->Draw();
            std::string savePath_phi = "hist_output/" + histName_phi + ".png";
            canvas_phi->SaveAs(savePath_phi.c_str());
            histMap_phi[histName_phi]->Write();


            histMap_W[histName_W]->Fit(f0,"R");
            //Hist1D_BE_sector0->GetListOfFunctions()->Add(f0);
            float temp_mag=f0->GetParameter(0);
            float temp_peak=f0->GetParameter(1);
            float temp_sigma=f0->GetParameter(2);
            TF1 *f1 = new TF1("f0","gaus",temp_peak - 1.0*temp_sigma, temp_peak + 0.8*temp_sigma);
            histMap_W[histName_W]->Fit(f1,"R");
            float temp_mag2=f1->GetParameter(0);
            float temp_peak2=f1->GetParameter(1);
            float temp_sigma2=f1->GetParameter(2);


            TF1* f0_a = new TF1("f0_a", bifurGauss, temp_peak2 - 1.2*temp_sigma2, temp_peak2 + 1.0*temp_sigma2, 4);
            f0_a->SetParameters(temp_peak2, temp_sigma2, temp_sigma2,temp_mag2); // Mean, Sigma_L, Sigma_R
            f0_a->SetParNames("Mean", "Sigma_L", "Sigma_R","Mag");
            histMap_W[histName_W]->Fit(f0_a,"R");
            histMap_W[histName_W]->GetListOfFunctions()->Add(f0_a);
            //float temp_peak2=f0_a->GetParameter(1);
            //float temp_sigma2=f0_a->GetParameter(2);

            


            //histMap_W[histName_W]->Scale(1.0 / histMap_W[histName_W]->Integral());
            TCanvas* canvas_W = new TCanvas((histName_W + "_canvas").c_str(), histTitle_W.c_str(), 800, 600);
            //canvas_W->SetLogy();
            histMap_W[histName_W]->Draw();
            histMap_W[histName_W]->Write();
            std::string savePath_W = "hist_output/" + histName_W + ".png";
            canvas_W->SaveAs(savePath_W.c_str());
            fprintf(fp,"%d\t%f\t%f\t%f\t%f\t%f\t",k+1,(theta_begin+l*deltatheta)*RadtoDeg,(theta_begin+(l+1)*deltatheta)*RadtoDeg,f0_a->GetParameter(0),f0_a->GetParameter(1),f0_a->GetParameter(2));
        }
        fprintf(fp,"\n");
    }
    Hist2D_W_theta_S1_selected->Write();
    Hist2D_W_theta_S1_corrected_selected->Write();

    Hist2D_W_theta_S2_selected->Write();
    Hist2D_W_theta_S2_corrected_selected->Write();

    Hist2D_W_theta_S3_selected->Write();
    Hist2D_W_theta_S3_corrected_selected->Write();

    Hist2D_W_theta_S4_selected->Write();
    Hist2D_W_theta_S4_corrected_selected->Write();

    Hist2D_W_theta_S5_selected->Write();
    Hist2D_W_theta_S5_corrected_selected->Write();

    Hist2D_W_theta_S6_selected->Write();
    Hist2D_W_theta_S6_corrected_selected->Write();

    fclose(fp);


    Electron_theta_store.clear();
    Electron_phi_store.clear();
    Electron_E_store.clear();
    Electron_sector_store.clear();
    Electron_vx_store.clear();
    Electron_vy_store.clear();
    Electron_vz_store.clear();

    gApplication->Terminate(0);
}