//
//  HW_iso.h
//  
//
//  Created by katy ghantous on 5/6/15.
//  Copyright 2015 __MyCompanyName__. All rights reserved.
//

#ifndef _HW_iso_h
#define _HW_iso_h



#include <iostream>
#include <complex>
#include <stdio.h>
#include<vector> 
#include<math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
using namespace std;




class V{
    
public:
    vector<double> kn;
    vector<double> phi_an,phi_bn,phi_cn, n_an, n_bn,n_cn;
    double g, k_0;
    double ap,bp, cp;
    double an,bn, cn;
    vector<int> indx;
    
    double q;
    double alphbar,nuZF;
    int Nshls;
    int ZF;
    int ipB,ipE;
    int inB,inE;
    
    void setV(int N, double gIn, double k_0i, double qIn,double alphbarIn,double nuZFIn,int ZFi);
    void set_alph_HW(double alph, double alph_n);
    void setCoef_HW();
    
    V();
    ~V();
    
};


int func_hmdsi(double t, const double y[], double dydt[], void *params);

double get_data_NAME(FILE *someFile, const char * nameData);


#endif
