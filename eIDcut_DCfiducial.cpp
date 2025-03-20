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
#include <filesystem>
#include <TVector3.h>

using namespace clas12;


namespace fs = std::filesystem;

int addFilesToChain(clas12root::HipoChain& chain, const std::string& directory) {
	int num_file=0;
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
        if (entry.path().extension() == ".hipo") {
            chain.Add(entry.path().c_str());
            std::cout << "Added file: " << entry.path() << std::endl;
			num_file++;
        }
    }
	cout<<"**************** "<<num_file<<" Files Loaded ****************"<<endl;
	return num_file;
}

bool Pass_DC_R1_fiducial_cut(float X,float Y,float Z, int S){
	float a1=0.556;
	float b1=-6.878;
	float a2=-0.560;
	float b2=7.482;
	float c=24.052;
	float theta=(S-1)*3.14159265359/3;
	float tempX=cos(theta)*X+sin(theta)*Y;
	float tempY=-sin(theta)*X+cos(theta)*Y;
	//float theta2=-25*3.14159265359/180;
	//float tempX2=cos(theta2)*tempX+sin(theta2)*Z;
	//float tempZ=-sin(theta2)*tempX+cos(theta2)*Z;
	//cout<<tempX2<<","<<tempY<<","<<Z<<":"<<atan(tempY/(tempX2+80))*180/3.14159265359<<":::::"<<S<<endl;
	//cout<<X<<","<<Y<<","<<Z<<":"<<atan(Y/X)*180/3.14159265359<<":"<<atan(Y/X)*180/3.14159265359-atan(tempY/tempX)*180/3.14159265359<<endl;
	if((tempY<a1*tempX+b1)&&(tempY>a2*tempX+b2)&&(tempX>c)){
		//cout<<1;
		return true;
	}
	else{
		//cout<<0;
		return false;
	}
}

bool Pass_DC_R2_fiducial_cut(float X,float Y,float Z, int S){
	float a1=0.578;
	float b1=-13.898;
	float a2=-0.577;
	float b2=14.851;
	float c=39.705;
	float theta=(S-1)*3.14159265359/3;
	float tempX=cos(theta)*X+sin(theta)*Y;
	float tempY=-sin(theta)*X+cos(theta)*Y;
	//float theta2=-25*3.14159265359/180;
	//float tempX2=cos(theta2)*tempX+sin(theta2)*Z;
	//float tempZ=-sin(theta2)*tempX+cos(theta2)*Z;
	//cout<<tempX2<<","<<tempY<<","<<Z<<":"<<atan(tempY/(tempX2+80))*180/3.14159265359<<":::::"<<S<<endl;
	//cout<<X<<","<<Y<<","<<Z<<":"<<atan(Y/X)*180/3.14159265359<<":"<<atan(Y/X)*180/3.14159265359-atan(tempY/tempX)*180/3.14159265359<<endl;
	if((tempY<a1*tempX+b1)&&(tempY>a2*tempX+b2)&&(tempX>c)){
		//cout<<1;
		return true;
	}
	else{
		//cout<<0;
		return false;
	}
}

bool Pass_DC_R3_fiducial_cut(float X,float Y,float Z, int S){
	float a1=0.591;
	float b1=-27.459;
	float a2=-0.588;
	float b2=26.912;
	float c=77.755;
	float theta=(S-1)*3.14159265359/3;
	float tempX=cos(theta)*X+sin(theta)*Y;
	float tempY=-sin(theta)*X+cos(theta)*Y;
	//float theta2=-25*3.14159265359/180;
	//float tempX2=cos(theta2)*tempX+sin(theta2)*Z;
	//float tempZ=-sin(theta2)*tempX+cos(theta2)*Z;
	//cout<<tempX2<<","<<tempY<<","<<Z<<":"<<atan(tempY/(tempX2+80))*180/3.14159265359<<":::::"<<S<<endl;
	//cout<<X<<","<<Y<<","<<Z<<":"<<atan(Y/X)*180/3.14159265359<<":"<<atan(Y/X)*180/3.14159265359-atan(tempY/tempX)*180/3.14159265359<<endl;
	if((tempY<a1*tempX+b1)&&(tempY>a2*tempX+b2)&&(tempX>c)){
		//cout<<1;
		return true;
	}
	else{
		//cout<<0;
		return false;
	}
}


int main(){
	clas12root::HipoChain chain;

	//chain.Add("/volatile/clas12/osg/yijie/job_8292/output/*.hipo");

	//chain.Add("/cache/clas12/rg-k/production/recon/fall2018/torus+1/6535MeV/pass2/v12/train/elastic/elastic_005969.hipo");

	chain.Add("/cache/clas12/rg-k/production/recon/fall2018/torus+1/6535MeV/pass2/v0/dst/recon/005951/rec_clas_005951.evio.00*.hipo");
	//chain.Add("/cache/clas12/rg-k/production/recon/fall2018/torus+1/6535MeV/pass2/v0/dst/recon/005995/rec_clas_005995.evio.*.hipo");//empty target

	//chain.Add("/w/hallb-scshelf2102/clas12/yijie/02DVCS_Analysis/elastic_genRad/elastgen/job_8530_tail/output/*.hipo");
	//chain.Add("/w/hallb-scshelf2102/clas12/yijie/02DVCS_Analysis/elastic_genRad/elastgen/job_8484_elastic/output/*.hipo");
	//chain.Add("/w/hallb-scshelf2102/clas12/yijie/job_8547/output/*.hipo");
	//chain.Add("/w/hallb-scshelf2102/clas12/yijie/job_8548/output/*.hipo");

	//int num_file_total=addFilesToChain(chain, "/cache/clas12/rg-k/production/recon/fall2018/torus+1/6535MeV/pass2/v0/dst/recon");

	//chain.Add("/volatile/clas12/osg/yijie/job_8719/output/inclusivedis6535_online*.hipo");

	//chain.Add("/work/clas12/yijie/02DVCS_Analysis/job_8597/output/elas6535_1000*.hipo");
	//chain.Add("/work/clas12/yijie/02DVCS_Analysis/job_8596/output/elastail6535_1000*.hipo");

	//int num_file_total=addFilesToChain(chain, "/cache/clas12/rg-k/production/recon/fall2018/torus+1/7546MeV/pass2/v0/dst/recon");
	//int num_file_total=addFilesToChain(chain, "/mss/clas12/rg-k/production/recon/fall2018/torus+1/6535MeV/pass2/v0/dst/recon");

	//chain.Add("/cache/clas12/rg-k/production/recon/fall2018/torus+1/7546MeV/pass2/v0/dst/recon/005864/rec_clas_005864.evio.00*.hipo");
	//chain.Add("/cache/clas12/rg-k/production/recon/fall2018/torus+1/7546MeV/pass2/v0/dst/recon/005867/*.hipo");

	//TFile * rootFile = TFile::Open("dst7546.root","recreate");
	TFile * rootFile = TFile::Open("dst6535_all_e_test.root","recreate");
	TTree * trackTree = new TTree("trackTree","elastic_scattering");
	
	TLorentzVector p4_electron;
	TLorentzVector p4_MC_electron;
	TLorentzVector p4_proton;
	TLorentzVector beam;
	TLorentzVector target;
	TVector3 p3_uncharged;

	//float beam_energy=10.604; //GeV
	//float beam_energy=7.546; //GeV
	float beam_energy=6.535; //GeV

	float Proton_mass=0.93827208816; //GeV
	float Electron_mass=0.00051099895000; //Gev
	
	float Q2;
	float Nu;
	float xB;
	float W;
	float MC_W;
	float Beam_energy;
    float PI=3.14159265359;

	float Electron_phi_begin_angle=-23*PI/180;
    float RadtoDeg=180.0/PI;
    float DegtoRad=PI/180.0;

    float deltaphi;

	float beam_charge;
	int run_number;
	int event_number;

	beam.SetXYZM(0,0,beam_energy,Electron_mass);
	target.SetXYZM(0,0,0,Proton_mass);

	std::vector<float> Electron_Px;
	std::vector<float> Electron_Py;
	std::vector<float> Electron_Pz;
	std::vector<float> Electron_vx;
	std::vector<float> Electron_vy;
	std::vector<float> Electron_vz;
	std::vector<float> Electron_theta;
	std::vector<float> Electron_phi;
	std::vector<float> Electron_E;
/*
	std::vector<float> MC_Electron_Px;
	std::vector<float> MC_Electron_Py;
	std::vector<float> MC_Electron_Pz;
	std::vector<float> MC_Electron_vx;
	std::vector<float> MC_Electron_vy;
	std::vector<float> MC_Electron_vz;
	std::vector<float> MC_Electron_theta;
	std::vector<float> MC_Electron_phi;
	std::vector<float> MC_Electron_E;
*/
	std::vector<int> Electron_DC_R1_Pid;
	std::vector<int> Electron_DC_R1_Did;
	std::vector<int> Electron_DC_R1_Layer;
	std::vector<int> Electron_DC_R1_Sector;
	std::vector<float> Electron_DC_R1_X;
	std::vector<float> Electron_DC_R1_Y;
	std::vector<float> Electron_DC_R1_Z;
	std::vector<float> Electron_DC_R1_Edge;
	std::vector<float> Electron_DC_R1_X2;
	std::vector<float> Electron_DC_R1_Y2;
	std::vector<float> Electron_DC_R1_Z2;
	int Electron_DC_R1_hit_size=0;

	std::vector<int> Electron_DC_R2_Pid;
	std::vector<int> Electron_DC_R2_Did;
	std::vector<int> Electron_DC_R2_Layer;
	std::vector<int> Electron_DC_R2_Sector;
	std::vector<float> Electron_DC_R2_X;
	std::vector<float> Electron_DC_R2_Y;
	std::vector<float> Electron_DC_R2_Z;
	std::vector<float> Electron_DC_R2_Edge;
	std::vector<float> Electron_DC_R2_X2;
	std::vector<float> Electron_DC_R2_Y2;
	std::vector<float> Electron_DC_R2_Z2;
	int Electron_DC_R2_hit_size=0;

	std::vector<int> Electron_DC_R3_Pid;
	std::vector<int> Electron_DC_R3_Did;
	std::vector<int> Electron_DC_R3_Layer;
	std::vector<int> Electron_DC_R3_Sector;
	std::vector<float> Electron_DC_R3_X;
	std::vector<float> Electron_DC_R3_Y;
	std::vector<float> Electron_DC_R3_Z;
	std::vector<float> Electron_DC_R3_Edge;
	std::vector<float> Electron_DC_R3_X2;
	std::vector<float> Electron_DC_R3_Y2;
	std::vector<float> Electron_DC_R3_Z2;
	int Electron_DC_R3_hit_size=0;

	std::vector<float> HTCC_Response_X;
	std::vector<float> HTCC_Response_Y;
	std::vector<float> HTCC_Response_Z;
	std::vector<float> HTCC_Response_hX;
	std::vector<float> HTCC_Response_hY;
	std::vector<float> HTCC_Response_hZ;

	std::vector<float> FTOF_Response_X;
	std::vector<float> FTOF_Response_Y;
	std::vector<float> FTOF_Response_Z;

	int Electron_sector;

	float Electron_Nphe=0;
	float ECAL_total_depE=0;
	float ECAL_SF=0;
	
	float Electron_ECAL_ECin_E;
	float Electron_ECAL_ECin_X;
	float Electron_ECAL_ECin_Y;
	float Electron_ECAL_ECin_Z;
	
	float Electron_ECAL_ECout_E;
	float Electron_ECAL_ECout_X;
	float Electron_ECAL_ECout_Y;
	float Electron_ECAL_ECout_Z;
	
	float Electron_ECAL_PCAL_E;
	float Electron_ECAL_PCAL_X;
	float Electron_ECAL_PCAL_Y;
	float Electron_ECAL_PCAL_Z;
	float Electron_ECAL_PCAL_U;
	float Electron_ECAL_PCAL_V;
	float Electron_ECAL_PCAL_W;
	//std::vector<float>
	//std::vector<float>

	std::vector<float> Proton_Px;
	std::vector<float> Proton_Py;
	std::vector<float> Proton_Pz;
	std::vector<float> Proton_vx;
	std::vector<float> Proton_vy;
	std::vector<float> Proton_vz;
	std::vector<float> Proton_theta;
	std::vector<float> Proton_phi;
	std::vector<float> Proton_E; 

	std::vector<int> Proton_DC_R1_Pid;
	std::vector<int> Proton_DC_R1_Did;
	std::vector<int> Proton_DC_R1_Layer;
	std::vector<int> Proton_DC_R1_Sector;
	std::vector<float> Proton_DC_R1_X;
	std::vector<float> Proton_DC_R1_Y;
	std::vector<float> Proton_DC_R1_Z;
	int Proton_DC_R1_hit_size=0;

	int Sum_DC_hit_size=0;

	trackTree->Branch("Q2",&Q2);
	trackTree->Branch("Nu",&Nu);
	trackTree->Branch("xB",&xB);
	trackTree->Branch("W",&W);
	trackTree->Branch("Beam_energy",&Beam_energy);
    trackTree->Branch("deltaphi",&deltaphi);

	trackTree->Branch("beam_charge",&beam_charge);
/*
	trackTree->Branch("MC_W",&MC_W);
*/
	trackTree->Branch("Electron_Px",&Electron_Px);
	trackTree->Branch("Electron_Py",&Electron_Py);
	trackTree->Branch("Electron_Pz",&Electron_Pz);
	trackTree->Branch("Electron_vx",&Electron_vx);
	trackTree->Branch("Electron_vy",&Electron_vy);
	trackTree->Branch("Electron_vz",&Electron_vz);
	trackTree->Branch("Electron_theta",&Electron_theta);
	trackTree->Branch("Electron_phi",&Electron_phi);
	trackTree->Branch("Electron_E",&Electron_E);
/*
	trackTree->Branch("MC_Electron_Px",&MC_Electron_Px);
	trackTree->Branch("MC_Electron_Py",&MC_Electron_Py);
	trackTree->Branch("MC_Electron_Pz",&MC_Electron_Pz);
	trackTree->Branch("MC_Electron_vx",&MC_Electron_vx);
	trackTree->Branch("MC_Electron_vy",&MC_Electron_vy);
	trackTree->Branch("MC_Electron_vz",&MC_Electron_vz);
	trackTree->Branch("MC_Electron_theta",&MC_Electron_theta);
	trackTree->Branch("MC_Electron_phi",&MC_Electron_phi);
	trackTree->Branch("MC_Electron_E",&MC_Electron_E);
*/
	trackTree->Branch("Electron_DC_R1_Pid",&Electron_DC_R1_Pid);
	trackTree->Branch("Electron_DC_R1_Did",&Electron_DC_R1_Did);
	trackTree->Branch("Electron_DC_R1_Layer",&Electron_DC_R1_Layer);
	trackTree->Branch("Electron_DC_R1_Sector",&Electron_DC_R1_Sector);
	trackTree->Branch("Electron_DC_R1_X",&Electron_DC_R1_X);
	trackTree->Branch("Electron_DC_R1_Y",&Electron_DC_R1_Y);
	trackTree->Branch("Electron_DC_R1_Z",&Electron_DC_R1_Z);
	trackTree->Branch("Electron_DC_R1_Edge",&Electron_DC_R1_Edge);
	trackTree->Branch("Electron_DC_R1_X2",&Electron_DC_R1_X2);
	trackTree->Branch("Electron_DC_R1_Y2",&Electron_DC_R1_Y2);
	trackTree->Branch("Electron_DC_R1_Z2",&Electron_DC_R1_Z2);
	trackTree->Branch("Electron_DC_R1_hit_size",&Electron_DC_R1_hit_size);

	trackTree->Branch("Electron_DC_R2_Pid",&Electron_DC_R2_Pid);
	trackTree->Branch("Electron_DC_R2_Did",&Electron_DC_R2_Did);
	trackTree->Branch("Electron_DC_R2_Layer",&Electron_DC_R2_Layer);
	trackTree->Branch("Electron_DC_R2_Sector",&Electron_DC_R2_Sector);
	trackTree->Branch("Electron_DC_R2_X",&Electron_DC_R2_X);
	trackTree->Branch("Electron_DC_R2_Y",&Electron_DC_R2_Y);
	trackTree->Branch("Electron_DC_R2_Z",&Electron_DC_R2_Z);
	trackTree->Branch("Electron_DC_R2_Edge",&Electron_DC_R2_Edge);
	trackTree->Branch("Electron_DC_R2_X2",&Electron_DC_R2_X2);
	trackTree->Branch("Electron_DC_R2_Y2",&Electron_DC_R2_Y2);
	trackTree->Branch("Electron_DC_R2_Z2",&Electron_DC_R2_Z2);
	trackTree->Branch("Electron_DC_R2_hit_size",&Electron_DC_R2_hit_size);

	trackTree->Branch("Electron_DC_R3_Pid",&Electron_DC_R3_Pid);
	trackTree->Branch("Electron_DC_R3_Did",&Electron_DC_R3_Did);
	trackTree->Branch("Electron_DC_R3_Layer",&Electron_DC_R3_Layer);
	trackTree->Branch("Electron_DC_R3_Sector",&Electron_DC_R3_Sector);
	trackTree->Branch("Electron_DC_R3_X",&Electron_DC_R3_X);
	trackTree->Branch("Electron_DC_R3_Y",&Electron_DC_R3_Y);
	trackTree->Branch("Electron_DC_R3_Z",&Electron_DC_R3_Z);
	trackTree->Branch("Electron_DC_R3_Edge",&Electron_DC_R3_Edge);
	trackTree->Branch("Electron_DC_R3_X2",&Electron_DC_R3_X2);
	trackTree->Branch("Electron_DC_R3_Y2",&Electron_DC_R3_Y2);
	trackTree->Branch("Electron_DC_R3_Z2",&Electron_DC_R3_Z2);
	trackTree->Branch("Electron_DC_R3_hit_size",&Electron_DC_R3_hit_size);

	trackTree->Branch("HTCC_Response_X",&HTCC_Response_X);
	trackTree->Branch("HTCC_Response_Y",&HTCC_Response_Y);
	trackTree->Branch("HTCC_Response_Z",&HTCC_Response_Z);
	trackTree->Branch("HTCC_Response_hX",&HTCC_Response_hX);
	trackTree->Branch("HTCC_Response_hY",&HTCC_Response_hY);
	trackTree->Branch("HTCC_Response_hZ",&HTCC_Response_hZ);

	trackTree->Branch("FTOF_Response_X",&FTOF_Response_X);
	trackTree->Branch("FTOF_Response_Y",&FTOF_Response_Y);
	trackTree->Branch("FTOF_Response_Z",&FTOF_Response_Z);

	trackTree->Branch("Electron_Nphe",&Electron_Nphe);
	trackTree->Branch("ECAL_total_depE",&ECAL_total_depE);
	trackTree->Branch("ECAL_SF",&ECAL_SF);

	trackTree->Branch("Electron_sector",&Electron_sector);
	
	trackTree->Branch("Electron_ECAL_ECin_E",&Electron_ECAL_ECin_E);
	trackTree->Branch("Electron_ECAL_ECin_X",&Electron_ECAL_ECin_X);
	trackTree->Branch("Electron_ECAL_ECin_Y",&Electron_ECAL_ECin_Y);
	trackTree->Branch("Electron_ECAL_ECin_Z",&Electron_ECAL_ECin_Z);

	trackTree->Branch("Electron_ECAL_ECout_E",&Electron_ECAL_ECout_E);
	trackTree->Branch("Electron_ECAL_ECout_X",&Electron_ECAL_ECout_X);
	trackTree->Branch("Electron_ECAL_ECout_Y",&Electron_ECAL_ECout_Y);
	trackTree->Branch("Electron_ECAL_ECout_Z",&Electron_ECAL_ECout_Z);

	trackTree->Branch("Electron_ECAL_PCAL_E",&Electron_ECAL_PCAL_E);
	trackTree->Branch("Electron_ECAL_PCAL_X",&Electron_ECAL_PCAL_X);
	trackTree->Branch("Electron_ECAL_PCAL_Y",&Electron_ECAL_PCAL_Y);
	trackTree->Branch("Electron_ECAL_PCAL_Z",&Electron_ECAL_PCAL_Z);
	trackTree->Branch("Electron_ECAL_PCAL_U",&Electron_ECAL_PCAL_U);
	trackTree->Branch("Electron_ECAL_PCAL_V",&Electron_ECAL_PCAL_V);
	trackTree->Branch("Electron_ECAL_PCAL_W",&Electron_ECAL_PCAL_W);

	trackTree->Branch("Proton_Px",&Proton_Px);
	trackTree->Branch("Proton_Py",&Proton_Py);
	trackTree->Branch("Proton_Pz",&Proton_Pz);
	trackTree->Branch("Proton_vx",&Proton_vx);
	trackTree->Branch("Proton_vy",&Proton_vy);
	trackTree->Branch("Proton_vz",&Proton_vz);
	trackTree->Branch("Proton_theta",&Proton_theta);
	trackTree->Branch("Proton_phi",&Proton_phi);
	trackTree->Branch("Proton_E",&Proton_E);

	trackTree->Branch("Proton_DC_R1_Pid",&Proton_DC_R1_Pid);
	trackTree->Branch("Proton_DC_R1_Did",&Proton_DC_R1_Did);
	trackTree->Branch("Proton_DC_R1_Layer",&Proton_DC_R1_Layer);
	trackTree->Branch("Proton_DC_R1_Sector",&Proton_DC_R1_Sector);
	trackTree->Branch("Proton_DC_R1_X",&Proton_DC_R1_X);
	trackTree->Branch("Proton_DC_R1_Y",&Proton_DC_R1_Y);
	trackTree->Branch("Proton_DC_R1_Z",&Proton_DC_R1_Z);
	trackTree->Branch("Proton_DC_R1_hit_size",&Proton_DC_R1_hit_size);

	trackTree->Branch("Sum_DC_hit_size",&Sum_DC_hit_size);

	int number_proton=0;
	int number_electron=0;
	int number_charged=0;

	trackTree->Branch("number_Proton",&number_proton);
	trackTree->Branch("number_Electron",&number_electron);
	trackTree->Branch("number_charged",&number_charged);

	int number_uncharged=0;
	std::vector<float> uncharged_sector;
	std::vector<float> uncharged_px;
	std::vector<float> uncharged_py;
	std::vector<float> uncharged_pz;
	std::vector<float> uncharged_P;

	std::vector<float> uncharged_P_different_sector;
	std::vector<float> uncharged_P_same_sector;

	std::vector<float> diff_sector;

	trackTree->Branch("number_uncharged",&number_uncharged);
	trackTree->Branch("uncharged_sector",&uncharged_sector);
	trackTree->Branch("uncharged_px",&uncharged_px);
	trackTree->Branch("uncharged_py",&uncharged_py);
	trackTree->Branch("uncharged_pz",&uncharged_pz);
	trackTree->Branch("uncharged_P",&uncharged_P);
	trackTree->Branch("uncharged_P_same_sector",&uncharged_P_same_sector);
	trackTree->Branch("uncharged_P_different_sector",&uncharged_P_different_sector);

	trackTree->Branch("diff_sector",&diff_sector);

	bool no_unacceptable_uncharged=true;
	

	int number_events=0;
	int number_events_selected=0;
	
	std::ofstream Pfile;
	//Pfile.open("selected_events.dat",std::ios::app);
	FILE *fp;
    fp = fopen("selected_events.dat","w");

	for (int ifile=0; ifile<chain.GetNFiles();ifile++) {
		cout<<"6535_incluesive_version"<<endl;
		//cout<<"Reading: "<<ifile<<"/"<<num_file_total<<" of all files..."<<endl;
		cout<<"Reading: "<<ifile<<" files..."<<endl;
		cout<<"Reading: "<<number_events_selected<<" events..."<<endl;
		clas12::clas12reader c12{chain.GetFileName(ifile).Data()};

		// main particle bank ========
        auto idx_RecPart=c12.addBank("REC::Particle");
        auto iPid=c12.getBankOrder(idx_RecPart,"pid");
		auto iPx=c12.getBankOrder(idx_RecPart,"px");
		auto iPy=c12.getBankOrder(idx_RecPart,"py");
		auto iPz=c12.getBankOrder(idx_RecPart,"pz");
		auto ivx=c12.getBankOrder(idx_RecPart,"vx");
		auto ivy=c12.getBankOrder(idx_RecPart,"vy");
		auto ivz=c12.getBankOrder(idx_RecPart,"vz");
		auto icharge=c12.getBankOrder(idx_RecPart,"charge");
/*
		auto idx_MCPart=c12.addBank("MC::Particle");
		auto iMCPid=c12.getBankOrder(idx_MCPart,"pid");
		auto iMCPx=c12.getBankOrder(idx_MCPart,"px");
		auto iMCPy=c12.getBankOrder(idx_MCPart,"py");
		auto iMCPz=c12.getBankOrder(idx_MCPart,"pz");
*/
		auto idx_RunConfig=c12.addBank("RUN::config");
		auto irun_num=c12.getBankOrder(idx_RunConfig,"run");
		auto ievent_num=c12.getBankOrder(idx_RunConfig,"event");


		auto idx_RecEvent=c12.addBank("REC::Event");
		auto ibeam_charge=c12.getBankOrder(idx_RecEvent,"beamCharge");

		// Read banks: with DC, CVT, FTOF, LTCC, HTCC, ECAL, CTOF, CND 
		auto idx_Traj = c12.addBank("REC::Traj");
        auto iPindex = c12.getBankOrder(idx_Traj,"pindex");
        auto iDetector = c12.getBankOrder(idx_Traj,"detector");
		auto iLayer = c12.getBankOrder(idx_Traj,"layer");
		//auto iSector = c12.getBankOrder(idx_Traj,"sector");
        auto iX = c12.getBankOrder(idx_Traj,"x");
        auto iY = c12.getBankOrder(idx_Traj,"y");
        auto iZ = c12.getBankOrder(idx_Traj,"z");
		auto iEdge = c12.getBankOrder(idx_Traj,"edge");

		//Read banks REC::Track
		auto idx_Track = c12.addBank("REC::Track");
		auto iSector = c12.getBankOrder(idx_Track,"sector");

		auto idx_Cal = c12.addBank("REC::Calorimeter"); //ECAL
        auto iPindex_Cal = c12.getBankOrder(idx_Cal,"pindex");
        auto iDetector_Cal = c12.getBankOrder(idx_Cal,"detector");
		auto iLayer_Cal = c12.getBankOrder(idx_Cal,"layer");
		auto iSector_Cal = c12.getBankOrder(idx_Cal,"sector");
        auto iX_Cal = c12.getBankOrder(idx_Cal,"x");
        auto iY_Cal = c12.getBankOrder(idx_Cal,"y");
        auto iZ_Cal = c12.getBankOrder(idx_Cal,"z");
		auto iEnergy_Cal = c12.getBankOrder(idx_Cal,"energy");
		auto ilu_Cal = c12.getBankOrder(idx_Cal,"lu");
		auto ilv_Cal = c12.getBankOrder(idx_Cal,"lv");
		auto ilw_Cal = c12.getBankOrder(idx_Cal,"lw");

		auto idx_Cherenkov = c12.addBank("REC::Cherenkov"); //HTCC
        auto iPindex_HTCC = c12.getBankOrder(idx_Cherenkov,"pindex");
        auto iDetector_HTCC = c12.getBankOrder(idx_Cherenkov,"detector");
        auto inphe_HTCC = c12.getBankOrder(idx_Cherenkov,"nphe");
		auto iX_HTCC = c12.getBankOrder(idx_Cherenkov,"x");
		auto iY_HTCC = c12.getBankOrder(idx_Cherenkov,"y");
		auto iZ_HTCC = c12.getBankOrder(idx_Cherenkov,"z");
		//auto ihX_HTCC = c12.getBankOrder(idx_Cherenkov,"hx");
		//auto ihY_HTCC = c12.getBankOrder(idx_Cherenkov,"hy");
		//auto ihZ_HTCC = c12.getBankOrder(idx_Cherenkov,"hz");

		auto idx_Scintillator = c12.addBank("REC::Scintillator"); //FTOF
		auto iPindex_Scintillator = c12.getBankOrder(idx_Scintillator,"pindex");
		auto iDetector_Scintillator = c12.getBankOrder(idx_Scintillator,"detector");
		auto iX_Scintillator = c12.getBankOrder(idx_Scintillator,"x");
		auto iY_Scintillator = c12.getBankOrder(idx_Scintillator,"y");
		auto iZ_Scintillator = c12.getBankOrder(idx_Scintillator,"z");


		while (c12.next()==true) {
			number_events++;

			for (int i = 0; i < c12.getBank(idx_RecEvent)->getRows(); i++){
				beam_charge=c12.getBank(idx_RecEvent)->getFloat(ibeam_charge,i);
			}

			for (int i = 0; i < c12.getBank(idx_RunConfig)->getRows(); i++){
				run_number=c12.getBank(idx_RunConfig)->getInt(irun_num,i);
				event_number=c12.getBank(idx_RunConfig)->getInt(ievent_num,i);
			}
/*
			for (int i = 0; i < c12.getBank(idx_MCPart)->getRows(); i++){
				auto temp_Px=c12.getBank(idx_MCPart)->getFloat(iPx,i);
				auto temp_Py=c12.getBank(idx_MCPart)->getFloat(iPy,i);
				auto temp_Pz=c12.getBank(idx_MCPart)->getFloat(iPz,i);
				auto temp_Pid=c12.getBank(idx_MCPart)->getInt(iPid,i);
				if (temp_Pid==11){		//electron
					MC_Electron_Px.push_back(temp_Px);
					MC_Electron_Py.push_back(temp_Py);
					MC_Electron_Pz.push_back(temp_Pz);
					MC_Electron_E;
					p4_MC_electron.SetXYZM(temp_Px,temp_Py,temp_Pz,Electron_mass);
					MC_Electron_theta.push_back(acos(p4_MC_electron.CosTheta()));
					MC_Electron_phi.push_back(p4_MC_electron.Phi());
					MC_Electron_E.push_back(p4_MC_electron.E());
				}
		
			}
*/
			float temp_electron_E=0;
			for (int i = 0; i < c12.getBank(idx_RecPart)->getRows(); i++){
				auto temp_Px=c12.getBank(idx_RecPart)->getFloat(iPx,i);
				auto temp_Py=c12.getBank(idx_RecPart)->getFloat(iPy,i);
				auto temp_Pz=c12.getBank(idx_RecPart)->getFloat(iPz,i);
				auto temp_vx=c12.getBank(idx_RecPart)->getFloat(ivx,i);
				auto temp_vy=c12.getBank(idx_RecPart)->getFloat(ivy,i);
				auto temp_vz=c12.getBank(idx_RecPart)->getFloat(ivz,i);
				auto temp_Pid=c12.getBank(idx_RecPart)->getInt(iPid,i);
				auto temp_charge=c12.getBank(idx_RecPart)->getInt(icharge,i);
				if(temp_charge!=0){
					number_charged++;
				}
				

				if (temp_Pid==11){		//electron
					number_electron++;
					Electron_Px.push_back(temp_Px);
					Electron_Py.push_back(temp_Py);
					Electron_Pz.push_back(temp_Pz);
					Electron_vx.push_back(temp_vx);
					Electron_vy.push_back(temp_vy);
					Electron_vz.push_back(temp_vz);
					p4_electron.SetXYZM(temp_Px,temp_Py,temp_Pz,Electron_mass);
					
					Electron_theta.push_back(acos(p4_electron.CosTheta()));
					Electron_phi.push_back(p4_electron.Phi());

					Electron_E.push_back(p4_electron.E());
					temp_electron_E=p4_electron.E();

					if (p4_electron.Phi()>Electron_phi_begin_angle){
                		Electron_sector=floor((p4_electron.Phi()*RadtoDeg+23)/60)+1;
            		}
            		else{
                		Electron_sector=floor((p4_electron.Phi()*RadtoDeg+23+360)/60)+1;
            		}
				}

				if(temp_charge==0){		//uncharged particle
					number_uncharged++;
					p3_uncharged.SetXYZ(temp_Px,temp_Py,temp_Pz);
					if(p3_uncharged.Mag()>=0.1){
						uncharged_px.push_back(temp_Px);
						uncharged_py.push_back(temp_Py);
						uncharged_pz.push_back(temp_Pz);
						float temp_theta=p3_uncharged.Theta();
						uncharged_P.push_back(p3_uncharged.Mag());
						float temp_uncharged_sector=0;

						if (p3_uncharged.Phi()>Electron_phi_begin_angle){
                			uncharged_sector.push_back(floor((p3_uncharged.Phi()*RadtoDeg+23)/60)+1);
							temp_uncharged_sector=floor((p3_uncharged.Phi()*RadtoDeg+23)/60)+1;
            			}
            			else{
                			uncharged_sector.push_back(floor((p3_uncharged.Phi()*RadtoDeg+23+360)/60)+1);
							temp_uncharged_sector=floor((p3_uncharged.Phi()*RadtoDeg+23+360)/60)+1;
            			}
						diff_sector.push_back(Electron_sector-temp_uncharged_sector);
						if(Electron_sector!=temp_uncharged_sector){
							no_unacceptable_uncharged=false;
							uncharged_P_different_sector.push_back(p3_uncharged.Mag());
						}
						if(Electron_sector==temp_uncharged_sector){
							uncharged_P_same_sector.push_back(p3_uncharged.Mag());
						}
					}
				}

				if (temp_Pid==2212){		//proton
					number_proton++;
					Proton_Px.push_back(temp_Px);
					Proton_Py.push_back(temp_Py);
					Proton_Pz.push_back(temp_Pz);
					Proton_vx.push_back(temp_vx);
					Proton_vy.push_back(temp_vy);
					Proton_vz.push_back(temp_vz);
					p4_proton.SetXYZM(temp_Px,temp_Py,temp_Pz,Proton_mass);
					
					Proton_theta.push_back(acos(p4_proton.CosTheta()));
					Proton_phi.push_back(p4_proton.Phi());

					Proton_E.push_back(p4_proton.E());
				}
			}

			for(int j=0;j<c12.getBank(idx_Traj)->getRows();j++){
				auto temp_traj_Pid = c12.getBank(idx_Traj)->getInt(iPindex,j);//the row index (starting from 0) of this response's particle in the REC::Particle bank
				auto temp_traj_Did = c12.getBank(idx_Traj)->getInt(iDetector,j);
				auto temp_traj_Layer = c12.getBank(idx_Traj)->getInt(iLayer,j);
				//auto temp_traj_Sector = c12.getBank(idx_Traj)->getInt(iSector,j);
				auto temp_traj_X = c12.getBank(idx_Traj)->getFloat(iX,j);
				auto temp_traj_Y = c12.getBank(idx_Traj)->getFloat(iY,j);
				auto temp_traj_Z = c12.getBank(idx_Traj)->getFloat(iZ,j);
				auto temp_traj_Edge = c12.getBank(idx_Traj)->getFloat(iEdge,j);
				if(temp_traj_Did == 6){ //DC
					if(c12.getBank(idx_RecPart)->getInt(iPid,temp_traj_Pid)==11){
						//Electron_DC_R1_Sector.push_back(c12.getBank(idx_Track)->getInt(iSector,temp_traj_Pid));
						if(temp_traj_Layer == 6){//R1
							Electron_DC_R1_Pid.push_back(c12.getBank(idx_RecPart)->getInt(iPid,temp_traj_Pid));
							Electron_DC_R1_Did.push_back(temp_traj_Did);
							Electron_DC_R1_Layer.push_back(temp_traj_Layer);
							Electron_DC_R1_Sector.push_back(c12.getBank(idx_Track)->getInt(iSector,temp_traj_Pid));
							float theta=(c12.getBank(idx_Track)->getInt(iSector,temp_traj_Pid)-1)*3.14159265359/3;
							float tempX=cos(theta)*temp_traj_X+sin(theta)*temp_traj_Y;
							float tempY=-sin(theta)*temp_traj_X+cos(theta)*temp_traj_Y;
							//float theta2=-25*3.14159265359/180;
							//float tempX2=cos(theta2)*tempX+sin(theta2)*temp_traj_Z;
							//float tempZ=-sin(theta2)*tempX+cos(theta2)*temp_traj_Z;
							Electron_DC_R1_X.push_back(temp_traj_X);
							Electron_DC_R1_Y.push_back(temp_traj_Y);
							Electron_DC_R1_Z.push_back(temp_traj_Z);
							Electron_DC_R1_Edge.push_back(temp_traj_Edge);
							Electron_DC_R1_X2.push_back(tempX);
							Electron_DC_R1_Y2.push_back(tempY);
							Electron_DC_R1_Z2.push_back(temp_traj_Z);
							Electron_DC_R1_hit_size++;
						}
						if(temp_traj_Layer == 18){//R2
							Electron_DC_R2_Pid.push_back(c12.getBank(idx_RecPart)->getInt(iPid,temp_traj_Pid));
							Electron_DC_R2_Did.push_back(temp_traj_Did);
							Electron_DC_R2_Layer.push_back(temp_traj_Layer);
							Electron_DC_R2_Sector.push_back(c12.getBank(idx_Track)->getInt(iSector,temp_traj_Pid));
							float theta=(c12.getBank(idx_Track)->getInt(iSector,temp_traj_Pid)-1)*3.14159265359/3;
							float tempX=cos(theta)*temp_traj_X+sin(theta)*temp_traj_Y;
							float tempY=-sin(theta)*temp_traj_X+cos(theta)*temp_traj_Y;
							//float theta2=-25*3.14159265359/180;
							//float tempX2=cos(theta2)*tempX+sin(theta2)*temp_traj_Z;
							//float tempZ=-sin(theta2)*tempX+cos(theta2)*temp_traj_Z;
							Electron_DC_R2_X.push_back(temp_traj_X);
							Electron_DC_R2_Y.push_back(temp_traj_Y);
							Electron_DC_R2_Z.push_back(temp_traj_Z);
							Electron_DC_R2_Edge.push_back(temp_traj_Edge);
							Electron_DC_R2_X2.push_back(tempX);
							Electron_DC_R2_Y2.push_back(tempY);
							Electron_DC_R2_Z2.push_back(temp_traj_Z);
							Electron_DC_R2_hit_size++;
						}
						if(temp_traj_Layer == 36){//R2
							Electron_DC_R3_Pid.push_back(c12.getBank(idx_RecPart)->getInt(iPid,temp_traj_Pid));
							Electron_DC_R3_Did.push_back(temp_traj_Did);
							Electron_DC_R3_Layer.push_back(temp_traj_Layer);
							Electron_DC_R3_Sector.push_back(c12.getBank(idx_Track)->getInt(iSector,temp_traj_Pid));
							float theta=(c12.getBank(idx_Track)->getInt(iSector,temp_traj_Pid)-1)*3.14159265359/3;
							float tempX=cos(theta)*temp_traj_X+sin(theta)*temp_traj_Y;
							float tempY=-sin(theta)*temp_traj_X+cos(theta)*temp_traj_Y;
							//float theta2=-25*3.14159265359/180;
							//float tempX2=cos(theta2)*tempX+sin(theta2)*temp_traj_Z;
							//float tempZ=-sin(theta2)*tempX+cos(theta2)*temp_traj_Z;
							Electron_DC_R3_X.push_back(temp_traj_X);
							Electron_DC_R3_Y.push_back(temp_traj_Y);
							Electron_DC_R3_Z.push_back(temp_traj_Z);
							Electron_DC_R3_Edge.push_back(temp_traj_Edge);
							Electron_DC_R3_X2.push_back(tempX);
							Electron_DC_R3_Y2.push_back(tempY);
							Electron_DC_R3_Z2.push_back(temp_traj_Z);
							Electron_DC_R3_hit_size++;
						}
						
					}
					if(c12.getBank(idx_RecPart)->getInt(iPid,temp_traj_Pid)==2212){
						Proton_DC_R1_Pid.push_back(c12.getBank(idx_RecPart)->getInt(iPid,temp_traj_Pid));
						Proton_DC_R1_Did.push_back(temp_traj_Did);
						Proton_DC_R1_Layer.push_back(temp_traj_Layer);
						//Proton_DC_R1_Sector.push_back(temp_traj_Sector);
						if(temp_traj_Layer == 6){
							Proton_DC_R1_Sector.push_back(c12.getBank(idx_Track)->getInt(iSector,temp_traj_Pid));
							Proton_DC_R1_X.push_back(temp_traj_X);
							Proton_DC_R1_Y.push_back(temp_traj_Y);
							Proton_DC_R1_Z.push_back(temp_traj_Z);
						}
						Proton_DC_R1_hit_size++;
					}
				}
			}

			ECAL_SF=0;
			ECAL_total_depE=0;
			float temp_ECAL_depE=0;
			float temp_ECin_depE=0;
			float temp_ECout_depE=0;
			float temp_PCAL_depE=0;
			for(int j=0;j<c12.getBank(idx_Cal)->getRows();j++){
				auto temp_cal_Pid = c12.getBank(idx_Cal)->getInt(iPindex_Cal,j);//the row index (starting from 0) of this response's particle in the REC::Particle bank
				auto temp_cal_Did = c12.getBank(idx_Cal)->getInt(iDetector_Cal,j);
				auto temp_cal_Layer = c12.getBank(idx_Cal)->getInt(iLayer_Cal,j);
				auto temp_cal_X = c12.getBank(idx_Cal)->getFloat(iX_Cal,j);
				auto temp_cal_Y = c12.getBank(idx_Cal)->getFloat(iY_Cal,j);
				auto temp_cal_Z = c12.getBank(idx_Cal)->getFloat(iZ_Cal,j);
				auto temp_cal_U = c12.getBank(idx_Cal)->getFloat(ilu_Cal,j);
				auto temp_cal_V = c12.getBank(idx_Cal)->getFloat(ilv_Cal,j);
				auto temp_cal_W = c12.getBank(idx_Cal)->getFloat(ilw_Cal,j);
				auto temp_cal_E = c12.getBank(idx_Cal)->getFloat(iEnergy_Cal,j);

				if(temp_cal_Did == 7){ //ECAL
					if(c12.getBank(idx_RecPart)->getInt(iPid,temp_cal_Pid)==11){
						//cout<<number_electron<<":"<<c12.getBank(idx_RecPart)->getInt(iPid,temp_cal_Pid)<<":"<<temp_cal_X<<":"<<temp_cal_Layer<<endl;
						if(temp_cal_Layer==1){//PCAL
							Electron_ECAL_PCAL_E=temp_cal_E;
							Electron_ECAL_PCAL_X=temp_cal_X;
							Electron_ECAL_PCAL_Y=temp_cal_Y;
							Electron_ECAL_PCAL_Z=temp_cal_Z;
							Electron_ECAL_PCAL_U=temp_cal_U;
							Electron_ECAL_PCAL_V=temp_cal_V;
							Electron_ECAL_PCAL_W=temp_cal_W;
							temp_PCAL_depE=temp_cal_E;
						}
						if(temp_cal_Layer==4){//EC-inner
							Electron_ECAL_ECin_E=temp_cal_E;
							Electron_ECAL_ECin_X=temp_cal_X;
							Electron_ECAL_ECin_Y=temp_cal_Y;
							Electron_ECAL_ECin_Z=temp_cal_Z;
							temp_ECin_depE=temp_cal_E;
						}
						if(temp_cal_Layer==7){//EC-outer
							Electron_ECAL_ECout_E=temp_cal_E;
							Electron_ECAL_ECout_X=temp_cal_X;
							Electron_ECAL_ECout_Y=temp_cal_Y;
							Electron_ECAL_ECout_Z=temp_cal_Z;
							temp_ECout_depE=temp_cal_E;
						}
					}
				}
			}
			temp_ECAL_depE=temp_PCAL_depE+temp_ECin_depE+temp_ECout_depE;
			ECAL_total_depE=temp_ECAL_depE;
			ECAL_SF=ECAL_total_depE/temp_electron_E;

			for(int j=0;j<c12.getBank(idx_Cherenkov)->getRows();j++){
				auto temp_HTCC_Pid = c12.getBank(idx_Cherenkov)->getInt(iPindex_HTCC,j);//the row index (starting from 0) of this response's particle in the REC::Particle bank
				auto temp_HTCC_Did = c12.getBank(idx_Cherenkov)->getInt(iDetector_HTCC,j);
				auto temp_HTCC_nphe = c12.getBank(idx_Cherenkov)->getFloat(inphe_HTCC,j);
				auto temp_HTCC_X = c12.getBank(idx_Cherenkov)->getFloat(iX_HTCC,j);
				auto temp_HTCC_Y = c12.getBank(idx_Cherenkov)->getFloat(iY_HTCC,j);
				auto temp_HTCC_Z = c12.getBank(idx_Cherenkov)->getFloat(iZ_HTCC,j);
				auto temp_HTCC_hX = c12.getBank(idx_Traj)->getFloat(iX,temp_HTCC_Pid);
				auto temp_HTCC_hY = c12.getBank(idx_Traj)->getFloat(iY,temp_HTCC_Pid);
				auto temp_HTCC_hZ = c12.getBank(idx_Traj)->getFloat(iZ,temp_HTCC_Pid);
				//cout<<temp_HTCC_hX<<endl;
				//auto temp_HTCC_hX = c12.getBank(idx_Cherenkov)->getFloat(ihX_HTCC,j);
				//auto temp_HTCC_hY = c12.getBank(idx_Cherenkov)->getFloat(ihY_HTCC,j);
				//auto temp_HTCC_hZ = c12.getBank(idx_Cherenkov)->getFloat(ihZ_HTCC,j);
				
				if(temp_HTCC_Did == 15){ //HTCC
					if(c12.getBank(idx_RecPart)->getInt(iPid,temp_HTCC_Pid)==11){
						//cout<<number_electron<<":"<<c12.getBank(idx_RecPart)->getInt(iPid,temp_HTCC_Pid)<<":"<<temp_HTCC_nphe<<endl;
						//cout<<c12.getBank(idx_RecPart)->getInt(iPid,temp_HTCC_Pid)<<":\t"<<temp_HTCC_X<<":\t"<<temp_HTCC_Y<<":\t"<<temp_HTCC_Z<<endl;
						Electron_Nphe=temp_HTCC_nphe;
						HTCC_Response_X.push_back(temp_HTCC_X);
						HTCC_Response_Y.push_back(temp_HTCC_Y);
						HTCC_Response_Z.push_back(temp_HTCC_Z);
						HTCC_Response_hX.push_back(temp_HTCC_hX);
						HTCC_Response_hY.push_back(temp_HTCC_hY);
						HTCC_Response_hZ.push_back(temp_HTCC_hZ);
					}
				}
			}

			for(int j=0;j<c12.getBank(idx_Scintillator)->getRows();j++){
				auto temp_FTOF_Pid = c12.getBank(idx_Scintillator)->getInt(iPindex_Scintillator,j);
				auto temp_FTOF_Did = c12.getBank(idx_Scintillator)->getInt(iDetector_Scintillator,j);
				auto temp_FTOF_X = c12.getBank(idx_Scintillator)->getFloat(iX_Scintillator,j);
				auto temp_FTOF_Y = c12.getBank(idx_Scintillator)->getFloat(iY_Scintillator,j);
				auto temp_FTOF_Z = c12.getBank(idx_Scintillator)->getFloat(iZ_Scintillator,j);

				if(temp_FTOF_Did == 12){ //FTOF
					if(c12.getBank(idx_RecPart)->getInt(iPid,temp_FTOF_Pid)==11){
						FTOF_Response_X.push_back(temp_FTOF_X);
						FTOF_Response_Y.push_back(temp_FTOF_Y);
						FTOF_Response_Z.push_back(temp_FTOF_Z);
					}
				}
			}

			Q2=-(beam-p4_electron).M2();
			Nu=beam.E()-p4_electron.E();
			xB=Q2/(2*Proton_mass*Nu);
			W=TMath::Sqrt((beam+target-p4_electron).M2());
			/*MC_W=TMath::Sqrt((beam+target-p4_MC_electron).M2());*/
            deltaphi=TMath::Abs(p4_proton.Phi()-p4_electron.Phi());
			Sum_DC_hit_size=Proton_DC_R1_hit_size+Electron_DC_R1_hit_size;
			
			Beam_energy=p4_electron.CosTheta()+p4_proton.CosTheta()*sin(acos(p4_electron.CosTheta()))/sin(acos(p4_proton.CosTheta()))-1;
			Beam_energy=Beam_energy*Proton_mass/(2*sin(acos(p4_electron.CosTheta())/2)*sin(acos(p4_electron.CosTheta())/2));

			bool DC_R1_fid=false;
			bool DC_R2_fid=false;
			bool DC_R3_fid=false;
			if(Electron_DC_R1_hit_size==1){
				DC_R1_fid=Pass_DC_R1_fiducial_cut(Electron_DC_R1_X[0],Electron_DC_R1_Y[0],Electron_DC_R1_Z[0],Electron_DC_R1_Sector[0]);
			}
			if(Electron_DC_R2_hit_size==1){
				DC_R2_fid=Pass_DC_R2_fiducial_cut(Electron_DC_R2_X[0],Electron_DC_R2_Y[0],Electron_DC_R2_Z[0],Electron_DC_R2_Sector[0]);
			}
			if(Electron_DC_R3_hit_size==1){
				DC_R3_fid=Pass_DC_R3_fiducial_cut(Electron_DC_R3_X[0],Electron_DC_R3_Y[0],Electron_DC_R3_Z[0],Electron_DC_R3_Sector[0]);
			}
			
			//if ((number_electron==1)&&(number_proton==1)&&(TMath::Abs(deltaphi-PI)<2*PI/180)&&(W>0.8)&&(W<1.05)){
			//if ((number_electron==1)&&(number_proton==1)&&(number_charged==2)&&(Electron_Nphe>=2)&&DC_R1_fid&&DC_R2_fid&&DC_R3_fid&&(W>0.0)&&(W<10)&&(TMath::Abs(deltaphi-PI)<2*PI/180)&&no_unacceptable_uncharged){
			//if ((number_electron==1)&&(number_proton==1)&&(number_charged==2)&&(Electron_Nphe>=2)&&DC_R1_fid&&DC_R2_fid&&DC_R3_fid&&(W>0.0)&&(W<10)&&(TMath::Abs(deltaphi-PI)<2*PI/180)){
			//if ((number_electron==1)&&(number_proton==1)&&(Electron_Nphe>=2)&&DC_R1_fid&&DC_R2_fid&&DC_R3_fid&&(W>0.8)&&(W<4)&&(TMath::Abs(deltaphi-PI)<2*PI/180)){
			//if ((number_electron==1)&&(Electron_Nphe>=2)&&DC_R1_fid&&DC_R2_fid&&DC_R3_fid&&(number_proton==1)&&(TMath::Abs(deltaphi-PI)<2*PI/180)&&(W>0.8)&&(W<1.05)){
			//if ((number_electron==1)&&(Electron_Nphe>=2)&&DC_R1_fid&&DC_R2_fid&&DC_R3_fid&&(number_proton==1)&&(number_charged==2)){
			//if ((number_electron==1)&&(Electron_Nphe>=2)&&DC_R1_fid&&DC_R2_fid&&DC_R3_fid&&(W>0.0)&&(W<1.5)){
			//if ((number_electron==1)&&(Electron_Nphe>=2)&&DC_R1_fid&&DC_R2_fid&&DC_R3_fid&&(W>0.0)&&(W<5.0)){
			if ((number_electron==1)&&(Electron_Nphe>=2)&&DC_R1_fid&&DC_R2_fid&&DC_R3_fid){
			//if ((number_electron==1)&&(Electron_Nphe>=2)){//&&(number_charged==2)
			//if ((number_electron==1)){
			//if (1){
				//cout<<DC_R1_fid<<endl;
				trackTree->Fill();
				number_events_selected++;
				//Pfile<<std::fixed<<run_number<<"\t";
				//Pfile<<std::fixed<<event_number<<"\t";
				//Pfile<<"\n";
				fprintf(fp,"%d\t%d\n",run_number,event_number);
			}
			Electron_Px.clear();
			Electron_Py.clear();
			Electron_Pz.clear();
			Electron_vx.clear();
			Electron_vy.clear();
			Electron_vz.clear();
			Electron_theta.clear();
			Electron_phi.clear();
			Electron_E.clear();
/*
			MC_Electron_Px.clear();
			MC_Electron_Py.clear();
			MC_Electron_Pz.clear();
			MC_Electron_vx.clear();
			MC_Electron_vy.clear();
			MC_Electron_vz.clear();
			MC_Electron_theta.clear();
			MC_Electron_phi.clear();
			MC_Electron_E.clear();
*/
			Electron_DC_R1_Pid.clear();
			Electron_DC_R1_Did.clear();
			Electron_DC_R1_Layer.clear();
			Electron_DC_R1_Sector.clear();
			Electron_DC_R1_X.clear();
			Electron_DC_R1_Y.clear();
			Electron_DC_R1_Z.clear();
			Electron_DC_R1_Edge.clear();
			Electron_DC_R1_X2.clear();
			Electron_DC_R1_Y2.clear();
			Electron_DC_R1_Z2.clear();

			Electron_DC_R2_Pid.clear();
			Electron_DC_R2_Did.clear();
			Electron_DC_R2_Layer.clear();
			Electron_DC_R2_Sector.clear();
			Electron_DC_R2_X.clear();
			Electron_DC_R2_Y.clear();
			Electron_DC_R2_Z.clear();
			Electron_DC_R2_Edge.clear();
			Electron_DC_R2_X2.clear();
			Electron_DC_R2_Y2.clear();
			Electron_DC_R2_Z2.clear();

			Electron_DC_R3_Pid.clear();
			Electron_DC_R3_Did.clear();
			Electron_DC_R3_Layer.clear();
			Electron_DC_R3_Sector.clear();
			Electron_DC_R3_X.clear();
			Electron_DC_R3_Y.clear();
			Electron_DC_R3_Z.clear();
			Electron_DC_R3_Edge.clear();
			Electron_DC_R3_X2.clear();
			Electron_DC_R3_Y2.clear();
			Electron_DC_R3_Z2.clear();

			HTCC_Response_X.clear();
			HTCC_Response_Y.clear();
			HTCC_Response_Z.clear();
			HTCC_Response_hX.clear();
			HTCC_Response_hY.clear();
			HTCC_Response_hZ.clear();

			FTOF_Response_X.clear();
			FTOF_Response_Y.clear();
			FTOF_Response_Z.clear();
			
			Electron_ECAL_ECin_E=0;
			Electron_ECAL_ECin_X=0;
			Electron_ECAL_ECin_Y=0;
			Electron_ECAL_ECin_Z=0;
			Electron_ECAL_ECout_E=0;
			Electron_ECAL_ECout_X=0;
			Electron_ECAL_ECout_Y=0;
			Electron_ECAL_ECout_Z=0;
			Electron_ECAL_PCAL_E=0;
			Electron_ECAL_PCAL_X=0;
			Electron_ECAL_PCAL_Y=0;
			Electron_ECAL_PCAL_Z=0;
			Electron_ECAL_PCAL_U=0;
			Electron_ECAL_PCAL_V=0;
			Electron_ECAL_PCAL_W=0;

			Proton_Px.clear();
			Proton_Py.clear();
			Proton_Pz.clear();
			Proton_vx.clear();
			Proton_vy.clear();
			Proton_vz.clear();
			Proton_theta.clear();
			Proton_phi.clear();
			Proton_E.clear();

			Proton_DC_R1_Pid.clear();
			Proton_DC_R1_Did.clear();
			Proton_DC_R1_Layer.clear();
			Proton_DC_R1_Sector.clear();
			Proton_DC_R1_X.clear();
			Proton_DC_R1_Y.clear();
			Proton_DC_R1_Z.clear();

			number_electron=0;
			number_proton=0;
			number_charged=0;

			number_uncharged=0;
			uncharged_sector.clear();
			uncharged_px.clear();
			uncharged_py.clear();
			uncharged_pz.clear();
			uncharged_P.clear();

			uncharged_P_same_sector.clear();
			uncharged_P_different_sector.clear();

			diff_sector.clear();

			no_unacceptable_uncharged=true;

			Electron_DC_R1_hit_size=0;
			Electron_DC_R2_hit_size=0;
			Electron_DC_R3_hit_size=0;

			Proton_DC_R1_hit_size=0;
			Sum_DC_hit_size=0;
			Electron_Nphe=0;
		}
	}
	cout<<number_events<<" events in total"<<endl;
	cout<<number_events_selected<<" events selected in total"<<endl;
	rootFile->Write();
	rootFile->Close();
	//Pfile.close();
	fclose(fp);
}