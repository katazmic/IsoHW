//
//  HW_iso.cpp
//  
//
//  Created by katy ghantous on 5/6/15.
//  Copyright 2015 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <complex>
#include <stdio.h>
#include<vector> 
#include<math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include "HW_iso.h"
using namespace std;




V::V(){};



void V::setV(int N, double gIn, double k_0i,double qIn,double alphbarIn,double nuZFIn,int ZFi){
    
    int i;
    double n; 
    
    
    q= qIn;
    alphbar= alphbarIn ;nuZF = nuZFIn;
    ZF = ZFi;
    
    
    ipB = 1; ipE = N; 
    
    inB = N+1; inE = 2*N;
    
    Nshls = N;
    g = gIn;
    k_0 = k_0i;
    
    
    
    for(i=0;i<N;i++){
        n = (double) i ;
        kn.push_back(pow(g,n)*k_0);
    }    
    
    indx.push_back(-1);//phi zonal flow
    for(i=0;i<2*N;i++){
        indx.push_back((int) i%N); 
    } 
    indx.push_back(-1); //n zonal flow
    
    
};


void V::set_alph_HW(double alph,double alph_n){
    
    
    ap = alph; bp = -alph/g; cp = alph/g/g;
    an = alph_n; bn = -alph_n/g; cn = alph_n/g/g;
    
    
};

void V::setCoef_HW(){
    int i;
    for(i=0;i<Nshls;i++){
        
        phi_an.push_back(pow(kn[i],2)*(g*g-1)/pow(g,7));
        phi_bn.push_back(pow(kn[i],2)*(g*g*g*g-1)/pow(g,2));
        phi_cn.push_back(pow(kn[i],2)*(g*g-1)*pow(g,5));
        
        n_an.push_back(pow(kn[i],2)/pow(g,3));
	    n_bn.push_back(pow(kn[i],2));
        n_cn.push_back(pow(kn[i],2)*pow(g,3));
        
    }
    
    
};


V::~V(){};


int func_hmdsi(double t, const double y[], double dydt[], void *params){
    
    
    V *VIn = (V*)params;
    FILE * fInp;
    
    
    complex<double> *dphidt=(complex<double> *)&dydt[0];
    
    
    
    int tid,nthreads;
    int i,n,N_V,fi,N;
    complex<double> C_phi_phi, C_phi_n,Sm,disp,Z_fl,frcng;
    double musm,mubg,Fn,Fp,C2,vn;
    double musmFac,mubgFac;
    
    complex<double> I (0.0,1.0);
    
    complex<double> SpZF (0.0,0.0);
    complex<double> SnZF (0.0,0.0);
    
    //-------------------------------------------------------
    //                  Setting up the equations    
    //-------------------------------------------------------
    
    /////   values from input file 
    
    
    fInp = fopen("INPUT","r");   
    C2 = get_data_NAME(fInp,"C");
    vn = get_data_NAME(fInp,"kappa"); 
    Fp = get_data_NAME(fInp,"forcing");       
    fclose(fInp);
    
    
    Fn= 0.0;
    musm =musmFac*pow(VIn->kn[0],6);  1.e-18;//instead of -13.. 10*pow(VIn->kn[0],6);
    mubg =mubgFac*pow(VIn->kn[VIn->Nshls-1],-4); 1.e-24; // instead of -15;//100*pow(VIn->kn[VIn->Nshls-1],-4);
    
    //C2 = 1; 
    
    N = VIn->Nshls;
    N_V = 2+2*N;
    fi = (int) (log(1/VIn->kn[0])/log(VIn->kn[1]/VIn->kn[0]));
    
    //vn = 1.;
    
    complex<double> *phiZF = (complex<double> *)&y[0];
    
    complex<double> *phi = (complex<double> *)&y[2*VIn->ipB];
    
    complex<double> *nn = (complex<double> *)&y[2*VIn->inB];
    
    complex<double> *nZF = (complex<double> *)&y[2*VIn->inE+2];
    
    for(n=0;n<N-1;n++){
        SpZF = SpZF + pow(VIn->kn[n],3)*((conj(phi[n])*conj(phi[n+1])) + (conj(phi[n])*conj(phi[n+1])));// FIX THIS
        SnZF = SnZF + VIn->kn[n]*((conj(phi[n])*conj(nn[n+1]) - conj(phi[n+1])*conj(nn[n])) + (conj(phi[n])*conj(nn[n+1]) - conj(phi[n+1])*conj(nn[n])));//FIX THIS
    }
    
    
    
    
    
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
#pragma omp parallel shared (phi,nn,dphidt,VIn,Sm,C2,vn,fi,N_V,Fp,Fn,musm,mubg) private(tid,i,n,C_phi_phi,C_phi_n,disp,Z_fl,frcng)  num_threads(2)
    {
#pragma omp for 
        for(i=1;i<N_V-1;i++){ 
            dphidt[i]=0;
            
            n = VIn->indx[i]; // 1+i%N
            
            
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
            //                  Potential
            /////////////////////////////////  PHI /////////////////////
            
            if((i>=(VIn->ipB)) && (i<=(VIn->ipE))){
                
                C_phi_phi = (((n>1) ? VIn->ap*VIn->phi_an[n]*(conj(phi[n-2])*conj(phi[n-1]))  :0 )
                             +( (n>0 &&n<(N-1)) ?   VIn->bp*VIn->phi_bn[n]*(conj(phi[n-1])*conj(phi[n+1]))  : 0 )
                             +((n<(N-2)) ?  VIn->cp*VIn->phi_cn[n]*(conj(phi[n+1])*conj(phi[n+2])) : 0)
                             );
                
                
                
                // dissipations 
                disp  =  -(musm*pow(VIn->kn[n],-6)+mubg*pow(VIn->kn[n],4))*phi[n];
                // forcing terms
                frcng = C2*(phi[n]-nn[n])/(-pow(VIn->kn[n],2));
                if(n==fi||n==fi+1)
                    frcng=frcng+Fp; 
                
                //zonal part  :alpbr*q/k*(k^2/g^2-q^2/g)*[phibr*ph0_(n-1)]                                                                                                                                        
                Z_fl = (n<(N-1) && VIn->ZF!=0 ?  VIn->alphbar*VIn->q/VIn->kn[n]*(VIn->kn[n]*VIn->kn[n]*(VIn->g*VIn->g) -VIn->q*VIn->q)*(conj(phiZF[0])*conj(phi[n+1])) : 0); // FIX THIS!!!
                
                dphidt[i] = C_phi_phi + frcng  + disp + Z_fl;
                
            }
            
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
            //                  Denisty 
            ///////////////////////////////////  N /////////////////////
            if((i>=(VIn->inB)) && (i<=(VIn->inE))){
                C_phi_n = (((n>1) ?   VIn->an*VIn->n_an[n]*(conj(nn[n-1])*conj(phi[n-2])-conj(phi[n-1])*conj(nn[n-2])) :0 )
                           +( (n>0 &&n<(N-1)) ?  VIn->bn*VIn->n_bn[n]*(conj(phi[n-1])*conj(nn[n+1])-conj(phi[n+1])*conj(nn[n-1])) : 0 )
                           +((n<(N-2)) ?   VIn->cn*VIn->n_cn[n]*(conj(phi[n+1])*conj(nn[n+2])-conj(phi[n+2])*conj(nn[n+1])) : 0)
                           );
                
                
                // dissipations 
                disp  =  -(musm*pow(VIn->kn[n],-6)+mubg*pow(VIn->kn[n],4))*nn[n];    
                // forcing terms
                frcng = C2*(phi[n]-nn[n]) - I*vn*VIn->kn[n]*phi[n];
                
                if(n==fi||n==fi+1)
                    frcng=frcng+Fn;
                //zonal part
                Z_fl =((n>0 && VIn->ZF!=0) ? VIn->alphbar*VIn->q*VIn->kn[n]*( (conj(phiZF[0])*conj(nn[n+1])-conj(nZF[0])*conj(phi[n+1]))) : 0);
                
                dphidt[i] =   C_phi_n + frcng  + disp + Z_fl;
            }
            //  cout<<t<<"\t"<<i<<"\t"<<n<<"\t"<<abs(dphidt[i])*abs(dphidt[i])<<"\n";      
            
        }
        
    }
    if(VIn->ZF!=0){
        dphidt[0] =  -VIn->alphbar*(VIn->g*VIn->g-1)/VIn->q*SpZF -VIn->nuZF*phiZF[0];//VIn->alphbar*(VIn->g*VIn->q-1)/VIn->q*SpZF  -VIn->nuZF*phiZF[0];                                  
        dphidt[N_V-1] = -VIn->alphbar*VIn->q*SnZF;// no -VIn->nuZF*x[i] for n;                                                                                                                         
    }
    
    
    return GSL_SUCCESS;
};



///////////////////////////  INPUT DATA /////////////////////////

double get_data_NAME(FILE *someFile, const char * nameData)
{
    double vlu;
    int dataExists;
    char dat;
    int dg,i;
    char data[20]; // = malloc(20*sizeof(char));
    
    dataExists = 0;
    dat = fgetc(someFile);
    
    while(dat != EOF && dataExists == 0)
    {
        i=0;     
        dat = fgetc(someFile);
        
        if(dat == nameData[i])
        {
            while(dat == nameData[i])
            {	    
                dat = fgetc(someFile);
                i++;
            }
            
            if(nameData[i] == '\0')
            { 
                dataExists = 1;		  
            }
        }
    }
    
    dg =0;  
    dat = fgetc(someFile);
    
    while((dat == '\t') || (dat == ' '))
    {
        dat = fgetc(someFile);		 
    }
    while((dat != '\t') && (dat != ' ') && (dat != EOF) && (dat != '\0')  )
    {
        data[dg] = dat;
        dat = fgetc(someFile);		 
        dg++;
    }
    
    vlu =(double) atof(data);
    
    if(dataExists ==0){
        printf("\n You goofed! Make sure the correct spelling of %s exists in the INPUT file. \n \n",nameData);
        abort();
    }
    
    return vlu;
    
    
}







