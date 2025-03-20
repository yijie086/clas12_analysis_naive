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

using namespace clas12;
bool isInRange(float num) {
    // 检查 num 是否在指定的范围内
    return (num >= -7.0f && num <= 15.0f) ||
           (num >= 53.0f && num <= 75.0f) ||
           (num >= 113.0f && num <= 135.0f) ||
           (num >= 173.0f && num <= 195.0f) ||
           (num >= 233.0f && num <= 255.0f) ||
           (num >= 293.0f && num <= 315.0f);
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

void Draw_Electron_Kinematics(){
    //TFile * f = TFile::Open("dst7546_MC.root","read");
    TFile * f = TFile::Open("/work/clas12/yijie/02DVCS_Analysis/auto_bin_analysis/dst6535_all_e_test.root","read");
    //TFile * f = TFile::Open("dst6535.root","read");
	TTree * t = (TTree*) f->Get("trackTree");

    //TFile * histFile = TFile::Open("hist_E7546.root","recreate");
    TFile * histFile = TFile::Open("hist_E6535_all_e_test_no_cut.root","recreate");

    float PI=3.14159265359;
    float RadtoDeg=180.0/PI;
    TH1 * Hist1D_Electron_theta = new TH1D("Electron_theta","Electron_theta",500,0,0.6*RadtoDeg);
    TH1 * Hist1D_Electron_phi = new TH1D("Electron_phi","Electron_phi",500,-30,340);
    TH1 * Hist1D_Electron_E = new TH1D("Electron_E","Electron_E",500,0,8);
    TH1 * Hist1D_Electron_vx = new TH1D("Electron_vx","Electron_vx",500,-20,20);
    TH1 * Hist1D_Electron_vy = new TH1D("Electron_vy","Electron_vy",500,-20,20);
    TH1 * Hist1D_Electron_vz = new TH1D("Electron_vz","Electron_vz",500,-20,20);

    TH2 * Hist2D_Electron_E_theta = new TH2D("Electron_E_theta","Electron_E_theta",500,5,30,500,0,7);
    TH2 * Hist2D_Electron_E_phi = new TH2D("Electron_E_phi","Electron_E_phi",500,-30,340,500,0,8);
    TH2 * Hist2D_Electron_theta_phi = new TH2D("Electron_theta_phi","Electron_theta_phi",500,-30,340,500,0,0.6*RadtoDeg);

    TH1 * Hist1D_Q2 = new TH1D("Q2","Q2",500,0,6);
    TH1 * Hist1D_W = new TH1D("W","W",500,0.0,1.5);

    TH2 * Hist2D_Q2_W = new TH2D("Q2_W","Q2_W",500,0.5,1.5,500,0,6);
    TH2 * Hist2D_theta_W = new TH2D("theta_W","theta_W",500,0.5,4.5,500,0,40);

    TH2 * Hist2D_Electron_E_theta_s1 = new TH2D("Electron_E_theta_s1","Electron_E_theta_s1",500,5,10,500,4,7);
    TH2 * Hist2D_Electron_E_theta_s2 = new TH2D("Electron_E_theta_s2","Electron_E_theta_s2",500,5,10,500,4,7);
    TH2 * Hist2D_Electron_E_theta_s3 = new TH2D("Electron_E_theta_s3","Electron_E_theta_s3",500,5,10,500,4,7);
    TH2 * Hist2D_Electron_E_theta_s4 = new TH2D("Electron_E_theta_s4","Electron_E_theta_s4",500,5,10,500,4,7);
    TH2 * Hist2D_Electron_E_theta_s5 = new TH2D("Electron_E_theta_s5","Electron_E_theta_s5",500,5,10,500,4,7);
    TH2 * Hist2D_Electron_E_theta_s6 = new TH2D("Electron_E_theta_s6","Electron_E_theta_s6",500,5,10,500,4,7);

    TH1 * Hist1D_Electron_theta_s1 = new TH1D("Electron_theta_s1","Electron_theta_s1",500,5,15);
    TH1 * Hist1D_Electron_theta_s2 = new TH1D("Electron_theta_s2","Electron_theta_s2",500,5,15);
    TH1 * Hist1D_Electron_theta_s3 = new TH1D("Electron_theta_s3","Electron_theta_s3",500,5,15);
    TH1 * Hist1D_Electron_theta_s4 = new TH1D("Electron_theta_s4","Electron_theta_s4",500,5,15);
    TH1 * Hist1D_Electron_theta_s5 = new TH1D("Electron_theta_s5","Electron_theta_s5",500,5,15);
    TH1 * Hist1D_Electron_theta_s6 = new TH1D("Electron_theta_s6","Electron_theta_s6",500,5,15);

    TH1 * Hist1D_Electron_DC_R1_X = new TH1D("Electron_DC_R1_X","Electron_DC_R1_X",500,-100,100);
    TH1 * Hist1D_Electron_DC_R1_Y = new TH1D("Electron_DC_R1_Y","Electron_DC_R1_Y",500,-100,100);

    TH1 * Hist1D_Electron_DC_R1_X_s1 = new TH1D("Electron_DC_R1_X_s1","Electron_DC_R1_X_s1",500,20,50);
    TH1 * Hist1D_Electron_DC_R1_Y_s1 = new TH1D("Electron_DC_R1_Y_s1","Electron_DC_R1_Y_s1",500,-15,15);
    TH1 * Hist1D_Electron_DC_R1_X_s2 = new TH1D("Electron_DC_R1_X_s2","Electron_DC_R1_X_s2",500,20,50);
    TH1 * Hist1D_Electron_DC_R1_Y_s2 = new TH1D("Electron_DC_R1_Y_s2","Electron_DC_R1_Y_s2",500,-15,15);
    TH1 * Hist1D_Electron_DC_R1_X_s3 = new TH1D("Electron_DC_R1_X_s3","Electron_DC_R1_X_s3",500,20,50);
    TH1 * Hist1D_Electron_DC_R1_Y_s3 = new TH1D("Electron_DC_R1_Y_s3","Electron_DC_R1_Y_s3",500,-15,15);
    TH1 * Hist1D_Electron_DC_R1_X_s4 = new TH1D("Electron_DC_R1_X_s4","Electron_DC_R1_X_s4",500,20,50);
    TH1 * Hist1D_Electron_DC_R1_Y_s4 = new TH1D("Electron_DC_R1_Y_s4","Electron_DC_R1_Y_s4",500,-15,15);
    TH1 * Hist1D_Electron_DC_R1_X_s5 = new TH1D("Electron_DC_R1_X_s5","Electron_DC_R1_X_s5",500,20,50);
    TH1 * Hist1D_Electron_DC_R1_Y_s5 = new TH1D("Electron_DC_R1_Y_s5","Electron_DC_R1_Y_s5",500,-15,15);
    TH1 * Hist1D_Electron_DC_R1_X_s6 = new TH1D("Electron_DC_R1_X_s6","Electron_DC_R1_X_s6",500,20,50);
    TH1 * Hist1D_Electron_DC_R1_Y_s6 = new TH1D("Electron_DC_R1_Y_s6","Electron_DC_R1_Y_s6",500,-15,15);

    TH2 * Hist2D_Electron_DC_R1_X_theta_s1 = new TH2D("Electron_DC_R1_X_theta_s1","Electron_DC_R1_X_theta_s1",500,5,15,500,20,50);
    TH2 * Hist2D_Electron_DC_R1_Y_theta_s1 = new TH2D("Electron_DC_R1_Y_theta_s1","Electron_DC_R1_Y_theta_s1",500,5,15,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_theta_s2 = new TH2D("Electron_DC_R1_X_theta_s2","Electron_DC_R1_X_theta_s2",500,5,15,500,20,50);
    TH2 * Hist2D_Electron_DC_R1_Y_theta_s2 = new TH2D("Electron_DC_R1_Y_theta_s2","Electron_DC_R1_Y_theta_s2",500,5,15,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_theta_s3 = new TH2D("Electron_DC_R1_X_theta_s3","Electron_DC_R1_X_theta_s3",500,5,15,500,20,50);
    TH2 * Hist2D_Electron_DC_R1_Y_theta_s3 = new TH2D("Electron_DC_R1_Y_theta_s3","Electron_DC_R1_Y_theta_s3",500,5,15,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_theta_s4 = new TH2D("Electron_DC_R1_X_theta_s4","Electron_DC_R1_X_theta_s4",500,5,15,500,20,50);
    TH2 * Hist2D_Electron_DC_R1_Y_theta_s4 = new TH2D("Electron_DC_R1_Y_theta_s4","Electron_DC_R1_Y_theta_s4",500,5,15,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_theta_s5 = new TH2D("Electron_DC_R1_X_theta_s5","Electron_DC_R1_X_theta_s5",500,5,15,500,20,50);
    TH2 * Hist2D_Electron_DC_R1_Y_theta_s5 = new TH2D("Electron_DC_R1_Y_theta_s5","Electron_DC_R1_Y_theta_s5",500,5,15,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_theta_s6 = new TH2D("Electron_DC_R1_X_theta_s6","Electron_DC_R1_X_theta_s6",500,5,15,500,20,50);
    TH2 * Hist2D_Electron_DC_R1_Y_theta_s6 = new TH2D("Electron_DC_R1_Y_theta_s6","Electron_DC_R1_Y_theta_s6",500,5,15,500,-15,15);

    TH2 * Hist2D_Electron_DC_R1_X_Y_s1 = new TH2D("Electron_DC_R1_X_Y_s1","Electron_DC_R1_X_Y_s1",500,20,50,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_Y_s2 = new TH2D("Electron_DC_R1_X_Y_s2","Electron_DC_R1_X_Y_s2",500,20,50,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_Y_s3 = new TH2D("Electron_DC_R1_X_Y_s3","Electron_DC_R1_X_Y_s3",500,20,50,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_Y_s4 = new TH2D("Electron_DC_R1_X_Y_s4","Electron_DC_R1_X_Y_s4",500,20,50,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_Y_s5 = new TH2D("Electron_DC_R1_X_Y_s5","Electron_DC_R1_X_Y_s5",500,20,50,500,-15,15);
    TH2 * Hist2D_Electron_DC_R1_X_Y_s6 = new TH2D("Electron_DC_R1_X_Y_s6","Electron_DC_R1_X_Y_s6",500,20,50,500,-15,15);


    float Q2_readout=0;
    float W_readout=0;
    std::vector<float> Q2_store;
    std::vector<float> W_store;
    t->SetBranchAddress("Q2",&Q2_readout);
    t->SetBranchAddress("W",&W_readout);
    for (int i = 0; i < t->GetEntries(); i++){
    //for (int i = 0; i < 1000000; i++){
        //t->SetBranchAddress("Q2",&Q2_readout);
        t->GetEntry(i);
        Hist1D_Q2->Fill(Q2_readout);
        Q2_store.push_back(Q2_readout);
        Hist1D_W->Fill(W_readout);
        W_store.push_back(W_readout);
    }  
    Hist1D_Q2->Write();
    Hist1D_W->Write();

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
    t->SetBranchAddress("Electron_theta",&Electron_theta_readout);
    t->SetBranchAddress("Electron_phi",&Electron_phi_readout);
    t->SetBranchAddress("Electron_E",&Electron_E_readout);
    t->SetBranchAddress("Electron_vx",&Electron_vx_readout);
    t->SetBranchAddress("Electron_vy",&Electron_vy_readout);
    t->SetBranchAddress("Electron_vz",&Electron_vz_readout);

    std::vector<float> *Electron_DC_R1_X_readout=0;
    std::vector<float> Electron_DC_R1_X_store;
    std::vector<float> *Electron_DC_R1_Y_readout=0;
    std::vector<float> Electron_DC_R1_Y_store;
    t->SetBranchAddress("Electron_DC_R1_X",&Electron_DC_R1_X_readout);
    t->SetBranchAddress("Electron_DC_R1_Y",&Electron_DC_R1_Y_readout);

    for (int i = 0; i < t->GetEntries() ; i++){
    //for (int i = 0; i < 1000000; i++){
        t->GetEntry(i);
        for (int j = 0; j < Electron_theta_readout->size(); j++){
            Electron_theta_store.push_back(Electron_theta_readout->at(j));
            Hist1D_Electron_theta->Fill(Electron_theta_readout->at(j)*RadtoDeg);
            Electron_E_store.push_back(Electron_E_readout->at(j));
            Hist1D_Electron_E->Fill(Electron_E_readout->at(j));
            //Electron_phi_store.push_back(Electron_phi_readout->at(j));
            //Hist1D_Electron_phi->Fill(Electron_phi_readout->at(j)*RadtoDeg);
            Electron_DC_R1_X_store.push_back(Electron_DC_R1_X_readout->at(j));
            Electron_DC_R1_Y_store.push_back(Electron_DC_R1_Y_readout->at(j));
            Hist1D_Electron_DC_R1_X->Fill(Electron_DC_R1_X_readout->at(j));
            Hist1D_Electron_DC_R1_Y->Fill(Electron_DC_R1_Y_readout->at(j));

            Electron_vx_store.push_back(Electron_vx_readout->at(j));
            Hist1D_Electron_vx->Fill(Electron_vx_readout->at(j));
            Electron_vy_store.push_back(Electron_vy_readout->at(j));
            Hist1D_Electron_vy->Fill(Electron_vy_readout->at(j));
            Electron_vz_store.push_back(Electron_vz_readout->at(j));
            Hist1D_Electron_vz->Fill(Electron_vz_readout->at(j));
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
        }    
    }
    Hist1D_Electron_theta->Write();
    Hist1D_Electron_phi->Write();
    Hist1D_Electron_E->Write();
    Hist1D_Electron_vx->Write();
    Hist1D_Electron_vy->Write();
    Hist1D_Electron_vz->Write();
    Hist1D_Electron_DC_R1_X->Write();
    Hist1D_Electron_DC_R1_Y->Write();

    for (int i = 0; i < Electron_E_store.size(); i++){
        //if(isInRange(Electron_phi_store[i]*RadtoDeg)){
        //Hist2D_BE_Electron_theta->Fill(Electron_theta_store[i]*RadtoDeg,BE_store[i]);
            Hist2D_Electron_E_theta->Fill(Electron_theta_store[i]*RadtoDeg,Electron_E_store[i]);
            Hist2D_Electron_E_phi->Fill(Electron_phi_store[i]*RadtoDeg,Electron_E_store[i]);
            Hist2D_Electron_theta_phi->Fill(Electron_phi_store[i]*RadtoDeg,Electron_theta_store[i]*RadtoDeg);
            Hist2D_Q2_W->Fill(W_store[i],Q2_store[i]);
            Hist2D_theta_W->Fill(W_store[i],Electron_theta_store[i]*RadtoDeg);
            if(get_sector(Electron_phi_store[i]*RadtoDeg)==1){
                Hist2D_Electron_E_theta_s1->Fill(Electron_theta_store[i]*RadtoDeg,Electron_E_store[i]);
                Hist1D_Electron_theta_s1->Fill(Electron_theta_store[i]*RadtoDeg);
                float theta=(0)*3.14159265359/3;
	            float tempX=cos(theta)*Electron_DC_R1_X_store[i]+sin(theta)*Electron_DC_R1_Y_store[i];
	            float tempY=-sin(theta)*Electron_DC_R1_X_store[i]+cos(theta)*Electron_DC_R1_Y_store[i];
                Hist2D_Electron_DC_R1_X_theta_s1->Fill(Electron_theta_store[i]*RadtoDeg,tempX);
                Hist2D_Electron_DC_R1_Y_theta_s1->Fill(Electron_theta_store[i]*RadtoDeg,tempY);
                Hist1D_Electron_DC_R1_X_s1->Fill(tempX);
                Hist1D_Electron_DC_R1_Y_s1->Fill(tempY);
                Hist2D_Electron_DC_R1_X_Y_s1->Fill(tempX,tempY);
            }
            if(get_sector(Electron_phi_store[i]*RadtoDeg)==2){
                Hist2D_Electron_E_theta_s2->Fill(Electron_theta_store[i]*RadtoDeg,Electron_E_store[i]);
                Hist1D_Electron_theta_s2->Fill(Electron_theta_store[i]*RadtoDeg);
                float theta=(1)*3.14159265359/3;
	            float tempX=cos(theta)*Electron_DC_R1_X_store[i]+sin(theta)*Electron_DC_R1_Y_store[i];
	            float tempY=-sin(theta)*Electron_DC_R1_X_store[i]+cos(theta)*Electron_DC_R1_Y_store[i];
                Hist2D_Electron_DC_R1_X_theta_s2->Fill(Electron_theta_store[i]*RadtoDeg,tempX);
                Hist2D_Electron_DC_R1_Y_theta_s2->Fill(Electron_theta_store[i]*RadtoDeg,tempY);
                Hist1D_Electron_DC_R1_X_s2->Fill(tempX);
                Hist1D_Electron_DC_R1_Y_s2->Fill(tempY);
                Hist2D_Electron_DC_R1_X_Y_s2->Fill(tempX,tempY);
            }
            if(get_sector(Electron_phi_store[i]*RadtoDeg)==3){
                Hist2D_Electron_E_theta_s3->Fill(Electron_theta_store[i]*RadtoDeg,Electron_E_store[i]);
                Hist1D_Electron_theta_s3->Fill(Electron_theta_store[i]*RadtoDeg);
                float theta=(2)*3.14159265359/3;
	            float tempX=cos(theta)*Electron_DC_R1_X_store[i]+sin(theta)*Electron_DC_R1_Y_store[i];
	            float tempY=-sin(theta)*Electron_DC_R1_X_store[i]+cos(theta)*Electron_DC_R1_Y_store[i];
                Hist2D_Electron_DC_R1_X_theta_s3->Fill(Electron_theta_store[i]*RadtoDeg,tempX);
                Hist2D_Electron_DC_R1_Y_theta_s3->Fill(Electron_theta_store[i]*RadtoDeg,tempY);
                Hist1D_Electron_DC_R1_X_s3->Fill(tempX);
                Hist1D_Electron_DC_R1_Y_s3->Fill(tempY);
                Hist2D_Electron_DC_R1_X_Y_s3->Fill(tempX,tempY);
            }
            if(get_sector(Electron_phi_store[i]*RadtoDeg)==4){
                Hist2D_Electron_E_theta_s4->Fill(Electron_theta_store[i]*RadtoDeg,Electron_E_store[i]);
                Hist1D_Electron_theta_s4->Fill(Electron_theta_store[i]*RadtoDeg);
                float theta=(3)*3.14159265359/3;
	            float tempX=cos(theta)*Electron_DC_R1_X_store[i]+sin(theta)*Electron_DC_R1_Y_store[i];
	            float tempY=-sin(theta)*Electron_DC_R1_X_store[i]+cos(theta)*Electron_DC_R1_Y_store[i];
                Hist2D_Electron_DC_R1_X_theta_s4->Fill(Electron_theta_store[i]*RadtoDeg,tempX);
                Hist2D_Electron_DC_R1_Y_theta_s4->Fill(Electron_theta_store[i]*RadtoDeg,tempY);
                Hist1D_Electron_DC_R1_X_s4->Fill(tempX);
                Hist1D_Electron_DC_R1_Y_s4->Fill(tempY);
                Hist2D_Electron_DC_R1_X_Y_s4->Fill(tempX,tempY);
            }
            if(get_sector(Electron_phi_store[i]*RadtoDeg)==5){
                Hist2D_Electron_E_theta_s5->Fill(Electron_theta_store[i]*RadtoDeg,Electron_E_store[i]);
                Hist1D_Electron_theta_s5->Fill(Electron_theta_store[i]*RadtoDeg);
                float theta=(4)*3.14159265359/3;
	            float tempX=cos(theta)*Electron_DC_R1_X_store[i]+sin(theta)*Electron_DC_R1_Y_store[i];
	            float tempY=-sin(theta)*Electron_DC_R1_X_store[i]+cos(theta)*Electron_DC_R1_Y_store[i];
                Hist2D_Electron_DC_R1_X_theta_s5->Fill(Electron_theta_store[i]*RadtoDeg,tempX);
                Hist2D_Electron_DC_R1_Y_theta_s5->Fill(Electron_theta_store[i]*RadtoDeg,tempY);
                Hist1D_Electron_DC_R1_X_s5->Fill(tempX);
                Hist1D_Electron_DC_R1_Y_s5->Fill(tempY);
                Hist2D_Electron_DC_R1_X_Y_s5->Fill(tempX,tempY);
            }
            if(get_sector(Electron_phi_store[i]*RadtoDeg)==6){
                Hist2D_Electron_E_theta_s6->Fill(Electron_theta_store[i]*RadtoDeg,Electron_E_store[i]);
                Hist1D_Electron_theta_s6->Fill(Electron_theta_store[i]*RadtoDeg);
                float theta=(5)*3.14159265359/3;
	            float tempX=cos(theta)*Electron_DC_R1_X_store[i]+sin(theta)*Electron_DC_R1_Y_store[i];
	            float tempY=-sin(theta)*Electron_DC_R1_X_store[i]+cos(theta)*Electron_DC_R1_Y_store[i];
                Hist2D_Electron_DC_R1_X_theta_s6->Fill(Electron_theta_store[i]*RadtoDeg,tempX);
                Hist2D_Electron_DC_R1_Y_theta_s6->Fill(Electron_theta_store[i]*RadtoDeg,tempY);
                Hist1D_Electron_DC_R1_X_s6->Fill(tempX);
                Hist1D_Electron_DC_R1_Y_s6->Fill(tempY);
                Hist2D_Electron_DC_R1_X_Y_s6->Fill(tempX,tempY);
            }
        //}
    }

    Hist2D_Electron_E_theta->Write();
    Hist2D_Electron_E_phi->Write();
    Hist2D_Electron_theta_phi->Write();
    Hist2D_Q2_W->Write();
    Hist2D_theta_W->Write();

    Hist2D_Electron_E_theta_s1->Write();
    Hist2D_Electron_E_theta_s2->Write();
    Hist2D_Electron_E_theta_s3->Write();
    Hist2D_Electron_E_theta_s4->Write();
    Hist2D_Electron_E_theta_s5->Write();
    Hist2D_Electron_E_theta_s6->Write();

    Hist1D_Electron_theta_s1->Write();
    Hist1D_Electron_theta_s2->Write();
    Hist1D_Electron_theta_s3->Write();
    Hist1D_Electron_theta_s4->Write();
    Hist1D_Electron_theta_s5->Write();
    Hist1D_Electron_theta_s6->Write();

    Hist2D_Electron_DC_R1_X_theta_s1->Write();
    Hist2D_Electron_DC_R1_Y_theta_s1->Write();
    Hist2D_Electron_DC_R1_X_theta_s2->Write();
    Hist2D_Electron_DC_R1_Y_theta_s2->Write();
    Hist2D_Electron_DC_R1_X_theta_s3->Write();
    Hist2D_Electron_DC_R1_Y_theta_s3->Write();
    Hist2D_Electron_DC_R1_X_theta_s4->Write();
    Hist2D_Electron_DC_R1_Y_theta_s4->Write();
    Hist2D_Electron_DC_R1_X_theta_s5->Write();
    Hist2D_Electron_DC_R1_Y_theta_s5->Write();
    Hist2D_Electron_DC_R1_X_theta_s6->Write();
    Hist2D_Electron_DC_R1_Y_theta_s6->Write();

    Hist1D_Electron_DC_R1_X_s1->Write();
    Hist1D_Electron_DC_R1_Y_s1->Write();
    Hist1D_Electron_DC_R1_X_s2->Write();
    Hist1D_Electron_DC_R1_Y_s2->Write();
    Hist1D_Electron_DC_R1_X_s3->Write();
    Hist1D_Electron_DC_R1_Y_s3->Write();
    Hist1D_Electron_DC_R1_X_s4->Write();
    Hist1D_Electron_DC_R1_Y_s4->Write();
    Hist1D_Electron_DC_R1_X_s5->Write();
    Hist1D_Electron_DC_R1_Y_s5->Write();
    Hist1D_Electron_DC_R1_X_s6->Write();
    Hist1D_Electron_DC_R1_Y_s6->Write();

    Hist2D_Electron_DC_R1_X_Y_s1->Write();
    Hist2D_Electron_DC_R1_X_Y_s2->Write();
    Hist2D_Electron_DC_R1_X_Y_s3->Write();
    Hist2D_Electron_DC_R1_X_Y_s4->Write();
    Hist2D_Electron_DC_R1_X_Y_s5->Write();
    Hist2D_Electron_DC_R1_X_Y_s6->Write();

    Electron_theta_store.clear();
    Electron_phi_store.clear();
    Electron_E_store.clear();
    Electron_vx_store.clear();
    Electron_vy_store.clear();
    Electron_vz_store.clear();

    gApplication->Terminate(0);
}