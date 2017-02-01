//#include "/home/haoyu/BHsim/CrossSections/BH_cross_sections/BH_cross_sections.h"
#include "BH_cross_sections.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>  //so thats why cmath included?
#include <time.h>    // for random number generatiion
#include <vector>
#include <string>
#include <sstream>


//#include "Riostream.h"
#include <TCanvas.h>
//#include "TH1D.h"
//#include "TGraphErrors.h"
#include <TGraph.h>
#include <TObject.h>


using namespace std;

// Original idea:
// 1: setup experimental data
// 2: use setup parameters to compute differential cross sections at:
//    theta = 3/5/7... deg, phi = 180*deg, fill histogram TH1D
// 3: use formula Asym = dSigmaA/dSigmaS to compute asymmetry curve
// 4: add standard deviation


// New porcess:
// 1. setup e_p, e_q as a pair of array
// 2. from keyboard input get theta_p/q, and phi
// 3. use those numbers to compute (relative cross sections?) cross sections
// 4. use formula Asym = dSigmaA/dSigmaS to compute asymmetry curve
// 5. plot results on the same TCanvas 

#define PI 3.1415926
#define DEG PI/180.0

// Note:
// 1. Units conventions in BH_cross_sections.h are in MeV and Radians

void TEST (){
	Double_t Theta;
	Double_t Phi = 180*DEG;

	cerr << "Theta for both particle is (in degree): " << endl;
	cin >> Theta;
	Theta *= DEG;
	//cerr << "Phi of positron is: " << endl;\
	//cin >> Phi
	Int_t NPts;
	cerr << "Numnber of points in TGraph is: " << endl;
	cin >> NPts;

	// Cross sections constructor needs 
	Double_t Egamma = 60.0;  // 60MeV beam
	Double_t Z = 92.0; // Uranium

	Double_t E_p[NPts];
	Double_t E_q[NPts];
	Double_t Del[NPts];
	Double_t Xsec[NPts];
	Double_t Xsec_s[NPts];
	Double_t Xsec_a[NPts];
	Double_t Asym[NPts];

	BH_cross_sections* XS = new BH_cross_sections(Z, Egamma);	

	for ( Int_t i = 0.0 ; i < NPts ; ++i){
		E_p[i] = (Double_t)i * Egamma/((Double_t)NPts);
		E_q[i] = Egamma - E_p[i];
		Del[i] = E_p[i] - E_q[i];
		Xsec[i] = XS->xsec_full(E_p[i], E_q[i], Theta, Theta, Phi);
		Xsec_s[i] = XS->xsec_s(E_p[i], E_q[i], Theta, Theta, Phi);
		Xsec_a[i] = XS->xsec_a(E_p[i], E_q[i], Theta, Theta, Phi);
		Asym[i] = Xsec_a[i]/Xsec_s[i];
	}

	TCanvas* TC = new TCanvas("TC"," Cross sections plots", 1500, 800);
	TC->Divide(2,2);
	vector<TGraph*> TG;
	TG.push_back(new TGraph(NPts, Del, Xsec));
	TG.push_back(new TGraph(NPts, Del, Xsec_s));
	TG.push_back(new TGraph(NPts, Del, Xsec_a));
	TG.push_back(new TGraph(NPts, Del, Asym));
	
	for( Int_t i = 0; i < 4; i++){
		TC->cd(i+1);
		TG.at(i)->Draw("AC*");
	}
}
