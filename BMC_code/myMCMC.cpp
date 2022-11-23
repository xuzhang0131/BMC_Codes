/*
 * main.cpp
 *
 *  Created on: Oct 15, 2018
 *      Author: zhangxu
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <bits/stdc++.h>
//#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <random>
//#include <boost/random.hpp>
//#include <boost/random/normal_distribution.hpp>

#include "Constants.h"
#include "genDRAM_GibbsWithMH.h"

using namespace std;
using namespace boost::numeric::odeint;

//starting values of parameters;
const double Optimize_orig [57]= {Kf_GLC_Transport,
    Kr_GLC_Transport,
    Kf_GLC_SerPool1,
    Kf_GLC_PYR,
    Kf_GLC_PRPP ,
    Kf_SerPool3_Synthesis                   ,
    Kf_SerPool3_Degradation                 ,
    Kf_GlyPool3_Synthesis                   ,
    Kf_GlyPool3_Degradation                 ,
    Kf_SerPool1_GlyPool1                      ,
    Kr_SerPool1_GlyPool1                      ,
    Kf_Ser_Transport                        ,
    Kr_Ser_Transport                        ,
    Kf_Gly_Transport                        ,
    Kr_Gly_Transport                        ,
    Kf_SerPool2_SerPoolMitochon             ,
    Kr_SerPool2_SerPoolMitochon             ,
    Kf_GlyPool2_GlyPoolMitochon             ,
    Kr_GlyPool2_GlyPoolMitochon             ,
    Kf_SerPoolMitochon_GlyPoolMitochon      ,
    Kr_SerPoolMitochon_GlyPoolMitochon      ,
    Kf_GlyPool1_CO2                         ,
    Kr_GlyPool1_CO2                          ,
    Kf_MethyleneTHF_FormylTHF_Cytoplasm     ,
    Kr_MethyleneTHF_FormylTHF_Cytoplasm     ,
    Kf_FormylTHF_Formate_Cytoplasm          ,
    Kr_FormylTHF_Formate_Cytoplasm          ,
    Kf_GlyPoolMitochon_CO2                  ,
    Kr_GlyPoolMitochon_CO2                  ,
    Kf_MethyleneTHF_FormylTHF_Mitochon      ,
    Kr_MethyleneTHF_FormylTHF_Mitochon      ,
    Kf_FormylTHF_Formate_Mitochon           ,
    Kr_FormylTHF_Formate_Mitochon           ,
    Kf_PRPP1_GAR                            ,
    Kf_GAR_FGAR                             ,
    Kf_FGAR_AMP                             ,
    Kf_IMP_AMP                              ,
    Kf_AMP_Degradation                      ,
    Kr_MethyleneTHF_Transport               ,
    Kf_MethyleneTHF_Transport               ,
    Kr_FormylTHF_Transport                  ,
    Kf_FormylTHF_Transport                  ,
    Kr_THF_Transport                        ,
    Kf_THF_Transport                        ,
    Kr_Formate_Transport                    ,
    Kf_Formate_Transport                    ,
    Kf_PRPP2_GAR                            ,
    Kf_PRPP3_GAR                            ,
    Kf_SerPool2_GlyPool2                    ,
    Kr_SerPool2_GlyPool2                    ,
    Kf_GlyPool2_CO2                         ,
    Kr_GlyPool2_CO2                         ,
    CellVolume                                ,
    MediumVolume                            ,
    TissueVolume                            ,
    CytoplasmVolume                         ,
    MitochonVolume                          };

double odeParameters[57] = {1.148784e-04, 1.851001e-01, 1.982816e-05, 3.470530e-01, 1.925921e-03, 5.223675e-06,
    2.586479e-02, 3.144930e-05, 5.218906e-04, 8.974777e+01, 8.051776e+00, 5.776139e-05,
    3.086536e-02, 3.121908e-05, 7.905642e-03, 2.412366e-03, 9.247545e-03, 3.138743e-03,
    1.435322e-02, 3.541535e+02, 8.662046e+00, 4.026545e+01, 1.768642e-02, 1.005021e+00,
    5.562776e-01, 3.315666e-01, 5.352891e+05, 2.505246e+02, 2.301281e-01, 2.607113e+00,
    1.974221e-01, 2.061819e+00, 2.411940e+04, 4.191163e+01, 4.281330e+04, 4.259399e+04,
    1.305757e-02, 5.887159e-04, 3.844901e-02, 5.783046e-03, 3.168669e+00, 2.850576e-01,
    7.568579e-01, 4.662015e-02, 3.153770e+00, 2.282994e-01, 2.504351e+01, 1.487911e+01,
    2.700213e+02, 2.265308e+00, 6.708198e+01, 1.047755e-02,
    CellVolume  ,  MediumVolume  ,  TissueVolume ,   CytoplasmVolume ,   MitochonVolume};

extern state_type Xd_flag;
state_type Xd_flag;

const double t0 (0*60), tf (24*60), df (60), MoleculeNumberInOneNanoMole (6.02214129e14);

double mmYobs [3][10] ={
    {log(16.72720275),log(0.008553777),log(0.008553777),log(0.008553777),log(133.3070782),log(0.008553777),log(0.316733723),log(2895.962684),log(8352.975818),log(0.008553777)},
    {log(14.613105),log(0.008553777),log(0.008553777),log(0.008553777),log(86.80758989),log(0.008553777),log(0.316733723),log(2330.656011),log(7553.664549),log(0.008553777)},
    {log(15.75248504),log(0.008553777),log(0.008553777),log(0.030020138),log(140.9276324),log(0.008553777),log(0.316733723),log(3425.228351),log(12198.03665),log(0.008553777)}
};
double Yode [10], startp [62], muprior [62], sigmaprior [62];
double beta_tmp[10], alpha_tmp[10], beta_err[10], alpha_err[10];
int observed [10];
vector<int> ifFit;
vector<double> sd_tune;

const int nruns (20000);
const  int nburnin (6000);
const int nthinning (1);

const int N_tmp(62);
const int n_Y(10);
const int n_err(10);
//extern const int fitsize(5);

//void write_forsave( const state_type &x , const double t )
//{
//    cout << t << '\t' << x[2] << '\t' << x[3] << endl;
//}

int main(int argc, char **argv)
{
    char *argu = argv[1];
    int seq_from = atoi(argu);
    char *argu2 = argv[2];
    int seq_to = atoi(argu2);
    
    int i,j; //for loop
    //initial data values
    state_type Xd;
    for(i=0; i<512; ++i) Xd[i] = 0;
    Xd[2] = 0.3 * 6.02214129e17;      // 300 n moles unlabeled cGLC, cGLC_13C0
    Xd[84] = 0.01 * 6.02214129e17; // cSerPool1_13C000D000 = 0.01 * 6.02214129e17; 0.01 u moles unlabeled cSerPool1
    Xd[438] = 0.0001 * 6.02214129e17; //cPRPP_13C0 = 0.0001 * 6.02214129e17; 0.1 n moles unlabeled cPRPP
    Xd[148] = 0.1 * 6.02214129e17; //cGlyPool1_13C00D00 = 0.1 * 6.02214129e17; 0.1 u moles unlabeled cGlyPool1
    Xd[420] = 0.00005 * 6.02214129e17; //cTHF = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled cTHF
    Xd[404] = 0.00005 * 6.02214129e17; //cMethyleneTHF_13C0D00 = 0.00005 * 6.02214129e17; 0.05 n moles unlabeled cMethyleneTHF (no carbon label, no deterium label)
    Xd[164] = 0.01 * 6.02214129e17; //cSerPool2_13C000D000 = 0.01 * 6.02214129e17;      % 0.01 u moles unlabeled cSerPool2 (no carbon label, no deterium label)
    Xd[228] = 0.1 * 6.02214129e17; //cGlyPool2_13C00D00 = 0.1 * 6.02214129e17;      % 0.1 u moles unlabeled cGlyPool2 (no carbon label, no deterium label)
    Xd[244] = 0.01 * 6.02214129e17; //cSerPool3_13C000D000 = 0.01 * 6.02214129e17;      % 0.01 u moles unlabeled cSerPool3 (no carbon label, no deterium label)
    Xd[308] = 0.1 * 6.02214129e17; //cGlyPool3_13C00D00 = 0.1 * 6.02214129e17;      % 0.1 u moles unlabeled cGlyPool3 (no carbon label, no deterium label)
    Xd[324] = 0.002 * 6.02214129e17; //mSerPool_13C000D000 = 0.002 * 6.02214129e17;      % 0.002 u moles unlabeled mSerPool (no carbon label, no deterium label)
    Xd[388] = 0.02 * 6.02214129e17; //mGlyPool_13C00D00 = 0.02 * 6.02214129e17;      % 0.02 u moles unlabeled mGlyPool (no carbon label, no deterium label)
    Xd[437] = 0.00005 * 6.02214129e17; //mTHF = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled cTHF
    Xd[421] = 0.00005 * 6.02214129e17; //mMethyleneTHF_13C0D00 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled mMethyleneTHF (no carbon label, no deterium label)
    Xd[412] = 0.00005 * 6.02214129e17; //cFormylTHF_13C0D0 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled cFormylTHF (no carbon label, no deterium label)
    Xd[416] = 0.00005 * 6.02214129e17; //cFormate_13C0D0 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled cFormate (no carbon label, no deterium label)
    Xd[433] = 0.00005 * 6.02214129e17; //mFormate_13C0D0 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled mFormate (no carbon label, no deterium label)
    Xd[429] = 0.00005 * 6.02214129e17; //mFormylTHF_13C0D0 = 0.00005 * 6.02214129e17;      % 0.05 n moles unlabeled mFormylTHF (no carbon label, no deterium label)
    Xd[440] = 0.0001 * 6.02214129e17; //cGAR_13C000 = 0.0001 * 6.02214129e17;      % 0.1 n moles unlabeled cGAR (no carbon label, no deterium label)
    Xd[448] = 0.0001 * 6.02214129e17; //cFGAR_13C0000D0 = 0.0001 * 6.02214129e17;      % 0.1 n moles unlabeled cFGAR (no carbon label, no deterium label)
    Xd[480] = 0.02 * 6.02214129e17; //cAMP_13C00000D00 = 0.02 * 6.02214129e17;      % 20 n moles unlabeled cAMP (no carbon label, no deterium label)
    Xd[1] = 236.310 * 6.02214129e17; //eGLC_13C1 = 236.310 * 6.02214129e17;      % 236310 n moles labeled eGLC
    Xd[4] = 2.987748945 * 6.02214129e17; //eSer_13C000D000 = 2.987748945 * 6.02214129e17;      % 2987.748945 n moles unlabeled eSer (no carbon label, no deterium label)
    Xd[68] = 4.496768784 * 6.02214129e17; //eGly_13C00D00 = 4.496768784 * 6.02214129 * 10^17;      % 4496.768784 n moles unlabeled eGly (no carbon label, no deterium label)
    
    
    for(i=0; i<NumberOfIsotopomers; ++i) Xd_flag[i] = Xd[i];//cout<<Xd_flag[i]<<'\t';}
    
    //    integrate(PurineSynthesis, Xd, t0, tf, df);
    //    cout << Xd[0] << '\t' << Xd[1] << '\t' << Xd[2] << '\t' << Xd[607] << endl;
    
    
    
    vector<vector<double> > content;
    vector<double> row;
    string line, word;
    
    fstream csvfile ("cppinput.csv", ios::in);
    if(csvfile.is_open())
    {
        while(getline(csvfile, line))
        {
            row.clear();
            
            stringstream str(line);
            
            while(getline(str, word, ','))row.push_back(stod(word));
            
            content.push_back(row);
        }
    }
    //first row is the data start.
    
    
    double out_for_save[seq_to - seq_from + 1][10]; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    //for the first sample result, should also change rownum here to start a new loop $$$$$$$$$$$$$$$$$$$$$$$$$$;
    for(int colnum=0;colnum<52;colnum++)
    {
        odeParameters[colnum] = content[seq_from - 1][colnum];
    }
    
    integrate(PurineSynthesis, Xd, t0, tf, df);//, write_forsave);
    
//    cout << Xd[2] << '\t' << Xd[3] << endl;
    
    state_type Yd; // the observations at 24h;
    for(i = 0; i < Yd.size() ; ++i ){
        Yd[i] = Xd[i]/MoleculeNumberInOneNanoMole; // unit is nano mole
    }
    
    double data_arr[34];
    double *data_obs = data_process(data_arr,Yd);
    
    j = 0;
    for (i = 0;i <34; ++i){
        if(data_obs[i] != 0) {observed[j] = i; ++j;}
    }
    
    for(i=0; i<10; ++i)
    {
        Yode[i] = log(data_obs[observed[i]]); // number of actually observed Y is 34.
        out_for_save[0][i] = Yode[i];
    }
    
    
    
    for(int rownum=seq_from;rownum<seq_to;rownum++)//$$$$$$$$$$$$$$$$$ adjust the rownum $$$$$$$$$$$$$$$$$$
    {
        for(int colnum=0;colnum<52;colnum++)
        {
            odeParameters[colnum] = content[rownum][colnum];
        }
        integrate(PurineSynthesis, Xd, t0, tf, df);//, write_forsave);

        for(i = 0; i < Yd.size() ; ++i ){
            Yd[i] = Xd[i]/MoleculeNumberInOneNanoMole; // unit is nano mole
        }

        data_obs = data_process(data_arr,Yd);

        j = 0;
        for (i = 0;i <34; ++i){
            if(data_obs[i] != 0) {observed[j] = i; ++j;}
        }

        for(i=0; i<10; ++i)
        {
            Yode[i] = log(data_obs[observed[i]]); // number of actually observed Y is 34.
            
            out_for_save[rownum - seq_from + 1][i] = Yode[i];
        }

//        cout << rownum + 1 << '\t';
    }


    string file_name0 = "ode_est_data10_2022_"+ std::to_string(seq_from) + "_" + std::to_string(seq_to)+ ".txt";;   //change txt name;$$$$$$$$$$$$$$$$$$$
    ofstream myfile0;
    myfile0.open (file_name0);//************************************************************

    for(i = 0; i<seq_to - seq_from + 1; ++i){

    for(j = 0; j<10; ++j){
        myfile0 << out_for_save[i][j] <<'\t';// exp(ps[file_save*(save_q-1)+i][ifFit[j]]);
    }
    myfile0 << '\n';
    }

    myfile0.close();//*********************************************************************
    
//    cout << endl;
    
    
    
    
    
    
    double Y_var[10] = {0.8994403,  3.8721452,  1.2737354, 49.1109880,  2.4902220,
                        1.4848506,  3.2646328,  9.8501356, 0.8180014,   1.5405173};
    
    //vector<double> startp;
    for(i=0; i<52; ++i)
    {
        startp[i] = log(OptimizeParameters[i]);
        muprior[i] = log(Optimize_orig[i]);
    }
    for(i=52; i<N_tmp; ++i)
    {
        startp[i] = Y_var[i-52];//startp.insert(startp.end(),Y_var.begin(),Y_var.end());
        //        muprior[i] = Y_var[i-52];
    }
    
    double sigma_par [52] = {log(0.01)-log(0.005),0-log(0.5),log(0.1)-log(0.05),0-log(0.5),log(0.01)-log(0.005),
        log(5e-4)-log(2.5e-4),log(1)-log(0.5),log(2e-4)-log(1e-4),log(0.02)-log(0.01),
        log(500)-log(225),log(10)-log(5),log(0.01)-log(0.005),log(0.1)-log(0.05),log(0.001)-log(0.0005),
        log(0.02)-log(0.01),log(0.02)-log(0.01),log(0.03)-log(0.015),log(0.01)-log(0.005),
        log(0.05)-log(0.025),log(1000)-log(450),log(20)-log(10),log(100)-log(45),log(2)-log(1),
        log(3)-log(1.5),log(10)-log(5),log(3)-log(1.5),log(1e7)-log(5e6),log(1e3)-log(500),
        log(4)-log(2),log(10)-log(5),0-log(0.5),log(10)-log(5),log(1e5)-log(5e4),log(100)-log(50),
        log(2e5)-log(0.9e5),log(2e5)-log(0.9e5),log(0.25)-log(0.15),log(0.01)-log(0.005),
        log(0.5)-log(0.25),log(0.3)-log(0.15),log(50)-log(25),log(5)-log(2.5),log(8)-log(4),
        log(1)-log(0.5),log(10)-log(5),log(10)-log(5),log(100)-log(50),log(50)-log(25),
        log(500)-log(225),log(500)-log(250),log(100)-log(50),log(3)-log(1.5)};
    for(i = 0; i<52; ++i) sigmaprior[i] = pow(sigma_par[i]*1.5,2);
    //    for(i = 52; i<86; ++i) sigmaprior[i] = 0.01*Y_var[i-52]*Y_var[i-52]; // sigma^2
    
    
    std::mt19937 rng;
    rng.seed(5);
    
    //with 3 replicates;
    
    //set up the prior for sigma^2
    double themean;
    Eigen::VectorXd thevar(n_Y);
    for(i=52; i<N_tmp; ++i)
    {
        themean = (mmYobs[0][i-52]+mmYobs[1][i-52]+mmYobs[2][i-52])/3;
        thevar(i-52)=((mmYobs[0][i-52]-themean)*(mmYobs[0][i-52]-themean)
                      +(mmYobs[1][i-52]-themean)*(mmYobs[1][i-52]-themean)
                      +(mmYobs[2][i-52]-themean)*(mmYobs[2][i-52]-themean) )/2; //unbiased sample variance;
    }
    
    double varmean = thevar.mean();
    Eigen::VectorXd varcent = thevar.rowwise()-thevar.colwise().mean();
    Eigen::MatrixXd medvar = (varcent.adjoint()*varcent)/double(n_Y-1);
    for(i=52; i<N_tmp; ++i) {
        muprior[i] = varmean;
        sigmaprior[i] = medvar(0,0);
    }
    //    cout << varmean <<'\t'<<medvar(0,0)<<endl;
    
    //assume the prior of var(epsilon) is inverse Gamma distribution
    for(i=0; i<n_Y; ++i){
        alpha_tmp[i] = pow(muprior[N_tmp-n_Y+i],2)/sigmaprior[N_tmp-n_Y+i] + 2;
        beta_tmp[i] = muprior[N_tmp-n_Y+i] * (alpha_tmp[i] - 1);
        alpha_err[i] = alpha_tmp[i] + 1.5;
        beta_err[i] = beta_tmp[i];
    }
    //the posterior is also inverse Gamma
    
    int indexFit [62];
    for(i=0; i<52; ++i) {indexFit[i] = 1; startp[i] = log(Optimize_orig[i]);}
    for(i=52; i<N_tmp; ++i) indexFit[i] = 0;  // only ode par part, no err var here;
    
    for(i=0; i<NumberOfIsotopomers; ++i) Xd[i] = Xd_flag[i];
    for(i=0; i<52; i++) odeParameters[i] = exp(startp[i]);
    integrate(PurineSynthesis, Xd, t0, tf, df);
    for(i = 0; i < Yd.size() ; ++i ){
        Yd[i] = Xd[i]/MoleculeNumberInOneNanoMole; // unit is nano mole
    }
    data_obs = data_process(data_arr,Yd);
    for(i=0; i<n_Y; ++i) Yode[i]= log(data_obs[observed[i]]);// cout<<Yode_tmp[i]<<'\t';}
    double sum_rand (0);
    for(i=52; i<52+n_err; ++i) {
        sum_rand = 0;
        for(j=0; j<3; ++j) sum_rand += pow(mmYobs[j][i-52]-Yode[i-52],2);
        std::gamma_distribution<double> gamma(alpha_err[i-52], 1/(beta_tmp[i-52]+sum_rand/2));
        beta_err[i-52] = beta_tmp[i-52]+sum_rand/2;
        startp[i] = 1/gamma(rng);
    }
    
    for(i = 0; i <N_tmp; ++i) if(indexFit[i] != 0) ifFit.push_back(i);
    
    for(i=0; i<ifFit.size();++i) sd_tune.push_back(2.38*2.38/ifFit.size());
    
//        genDRAM_GibbsWithMH();
    
    return 0;
} // end of main function;



