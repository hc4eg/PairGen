#include "BH_cross_sections.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

double solve_e_positron(double e_gamma, double m_e, double m_r,
	double e_p, double th_p, double th_q, double phi_q);
// Double check correctness
double solve_e_positron2(double e_gamma, double m_e, double m_r,
	double e_p, double th_p, double th_q, double phi_q);
double ke_recoil(double e_gamma, double e_p, double e_q, double m_e, double m_r,
	double th_p, double th_q, double phi_q);
double angle_gen(double th_min, double th_max, double& th, double& phi);

double randfloat(void);
//double RandGaussian(double s);
//vector<string> split(string str, char delimiter);

int main(int argc, char *argv[])
{
	double pi=3.14159265358979323846264338;
	double deg = pi/180.;
	double m_e = 0.511;

	double e_gamma = 60.; // MeV
	double Z = 92; // Z uranium
	double A = 238; // A u  U238
	double m_r = A*931.5; // recoil mass
	
	double th_p, th_q, delta, e_p, e_q, ke_r, phi_q;

	char filename[100];

	BH_cross_sections * xs = new BH_cross_sections(Z, e_gamma);

	srandom((unsigned int) time(NULL));
	long events_to_do = 1000000;
//	long events_to_do = 1000;
	if(argc > 1) events_to_do = atoi(argv[1]);
	long total_count = 0;

	// experiment parameters
	double central_energy = 30.;
	double max_percent = 40.;
	double min_percent = -40.;
	double max_eng = central_energy * (1.+max_percent/100.);
	double min_eng = central_energy * (1.+min_percent/100.);
	double max_t_theta = 15.*deg;
	double max_t_phi = 2.2*deg;
	// Note: t_theta and t_phi ARE NOT SPHERICAL ANGLES, they're just angles atan(x/z) and atan(y/z)
	// are used to calculate x and y from z, which determined by acceptance.
	double min_polar = 4.*deg;
	double th_min = min_polar;
	double th_max = atan( sqrt(tan(max_t_theta)*tan(max_t_theta) + tan(max_t_phi)*tan(max_t_phi)) );

	ifstream runin;
	ofstream runout;
	int runno;
	// if the file "events.runno" exists - use its contents as the run number
	runin.open("events.runno");
	if(runin.good())
		{
		runin>>runno;
		runin.close();
		}
	else
		{
		runno = 0;
		}
	runno++;
	runout.open("events.runno", ios::trunc);
	if(runout.good())
		{
		runout << runno << endl;
		runout.close();
		}
	else
		{
		cout << "Cannot open events.runno" << endl;
		}



	// double max_cos_polar = 1. - cos(max_t_theta);
	double t_theta_p, t_theta_q, t_phi_p, t_phi_q;
	double phi_p_actual, phi_q_actual;

	FILE *fp;
	sprintf(filename, "events.run%.3d.dat",runno);
	fp = fopen(filename, "w");

	FILE *stdfp;
	sprintf(filename, "output.run%.3d.dat",runno);
	stdfp = fopen(filename, "w");

	fprintf(stdfp, "Run %d: Events to do = %d\n", runno, events_to_do);
	fflush(stdfp);

	double maximum_xsec = 5.e-6;
	double max_e_p, max_e_q, max_th_p, max_th_q, max_phi_q;
	int found_in_last = 0;
	long loops = 0;
	long starttime = time(NULL);
	while(total_count < events_to_do)
		{
		loops++;
		// choose random energies and angles for electron and positron
		// choose electron angles within spectrometer acceptance

		/*	
		do	{
			// It looks like here is to generate angle with uniformly random spherical distribution
			// in theta  range (0, max_t_theta) and phi range (0, 360*deg)
			double cos_polar = 1. - randfloat()*max_cos_polar;
			th_p = acos(cos_polar);
			phi_p_actual = randfloat()*360.*deg;
			double xx = sin(th_p)*cos(phi_p_actual);
			double yy = sin(th_p)*sin(phi_p_actual);
			double zz = cos_polar;
			t_theta_p = atan(xx/zz);
			t_phi_p = atan(yy/zz);
		// It only use t_theta_p and t_phi_p to check max_t_theta and max_t_phi condition. THEY'RE NOT SPHERICAL ANGLES
		// Later t_theta_p and t_phi_p are write into dat file.
			} while ( abs(t_theta_p) > max_t_theta
				|| abs(t_phi_p) > max_t_phi
		//Q:  Why choose t_phi_p ( phi_p? ) < 2.2*deg? ie. Why there is no large phi coming out?
		//A: t_phi_p/q ARE NOT SPHERICAL ANGLES
				|| th_p < min_polar);
		//Q:  Here th_p < min_polar disposed, but why in .dat file it can have nagative value? 
		// (A: .dat file prints t_theta_p than th_p)
		// Here th_p/q generate from acos() function, it must in range (0,acos(max_cos_polar)), and it's positive.
		// But not for t_theta_p/q.

		// Q; why needs a min_polar?
		// A: Based on Rob's obervation, small polar angles have bad measurable asymmetryi.
		// choose positron angles within spectrometer acceptance by adding max_t_theta and max_t_phi
		do	{
			double cos_polar = 1. - randfloat()*max_cos_polar;
			th_q = acos(cos_polar);
			phi_q_actual = randfloat()*360.*deg;
			double xx = sin(th_q)*cos(phi_q_actual);
			double yy = sin(th_q)*sin(phi_q_actual);
			double zz = cos_polar;
			t_theta_q = atan(xx/zz);
			t_phi_q = atan(yy/zz);
			} while ( abs(t_theta_q) > max_t_theta
				|| abs(t_phi_q) > max_t_phi
				|| th_q < min_polar);
		*/
	
		// Now randomly generate a pair of directions 
		// (1) Uniform random in phase space
		// (2) Within accpectance
		// Using function Angle_Gen 

		// Q: If phi_q is determined below, then phi_q range is (0.deg , 4.4deg) and (355.6deg, 360deg), is this intented ?
		// A: I didn't see the difference between t_phi_p/q and phi_p/q_actuall

		// Note: I found out the issue:
		// 1. The spherical angles here are th_p/q,  phi_q_actual and phi_q.

		// 2. t_theta_p/q are from atan(xx/zz) is X projection/ Z projection, also t_phi_p/q.
		// They're used as X = Z*tan(t_theta_p/q) and Y = Z*tan(t_phi_p/q). So are not spherical angle

		// 3. From 1.2. only th_p/q and phi_p/q_actual or phi_q should used in function solve_e_positron() and ke_recoil(), which is correct.
		angle_gen(th_min, th_max, th_p, phi_p_actual);
		angle_gen(th_min, th_max, th_q, phi_q_actual);
	
		phi_q = phi_q_actual - phi_p_actual;
		if(phi_q < 0.) phi_q += 360.*deg;

		// choose electron energy within spectrometer acceptance
		// ke_p from 18 to 42 MeV?
		double ke_p = min_eng + (max_eng-min_eng)*randfloat();
		e_p = ke_p + m_e;

		// solve for allowed positron energy
		// If angle information is not concerned, e_p + e_q in .dat file will almost always be 60.0 MeV
		// since the recoil nucleus is large, recoil KE is very small.
	//	double e_q = solve_e_positron(e_gamma, m_e, m_r, e_p, th_p, th_q, phi_q);
		double e_q = solve_e_positron2(e_gamma, m_e, m_r, e_p, th_p, th_q, phi_q);
	//	cerr << "e_q  = " << e_q << ", and e_q2 = " << e_q2 << endl;		
	
		// calculate the cross section
		double xsec = xs->xsec_full(e_p, e_q, th_p, th_q, phi_q);
		//fp << phi_q/deg << endl;
		double ratio = xsec/maximum_xsec;
		if(ratio > 1.)
			{
			fprintf(stdfp, "Cross Section = %.8g > maximum cross section = %.8g\n", xsec, maximum_xsec);
			// (I'm dubious about changing maximum_xsec during simulation, though chances are small)
			// but the former generated pairs data will have higher frequency
			// need a function to find maximum_xsec in the region
			maximum_xsec = xsec;
			fflush(fp);
			fflush(stdfp);
			}
		if( randfloat() > ratio) continue;

		// write the event to a file
		double ke_q = e_q - m_e;
		//fprintf(fp, "%.4f %.4f %.4f %.4f %.4f %.4f %.4g\n",
		//		ke_p, ke_q, t_theta_p/deg, t_phi_p/deg, t_theta_q/deg, t_phi_q/deg, xsec);
		fprintf(fp, "%.4f %.4f %.4f %.4f %.4f %.4f %.4g\n",
				ke_p, ke_q, th_p/deg, phi_p_actual/deg, th_q/deg, phi_q_actual/deg, xsec);
		fflush(fp);
		total_count++;
		cerr << "phi_p_actual = " << phi_p_actual/deg << " ,t_phi_p = " << t_phi_p/deg << endl
		     <<	", phi_q_actual = " << phi_q_actual/deg <<" , t_phi_q = " << t_phi_q/deg << endl
		     <<	", phi_q = " << phi_q/deg << endl;


		// Rate: ~ 3hr 16000 events.
		if(total_count%10 == 0)
			{
			double percent = double(total_count)/double(events_to_do)*100.;
			long now = time(NULL);
			long seconds = now - starttime;
			fprintf(stdfp, "Total Counts = %d (%.0f%%) Loops since last = %d Runtime = %d sec.\n", total_count, percent, loops, seconds);
			fflush(fp);
			fflush(stdfp);
			loops = 0;
			}
		}
	fclose(fp);
	
	fprintf(stdfp, "\nFinished: total_count = %d\n", total_count);
	long now = time(NULL);
	long seconds = now - starttime;
	fprintf(stdfp, "Run time = %d sec\n", seconds);
	fclose(stdfp);
}
double
solve_e_positron(double e_gamma, double m_e, double m_r,
	double e_p, double th_p, double th_q, double phi_q)
{
	// solve for allowed positron energy
	double ke_r = 0.1; // first guess
	double e_q = e_gamma - e_p - ke_r;
	// solve for actual ke_r and hence actual e_q
	double tol = 1.0e-7;
	double diff;
	double step = 0.01;
	// first guess
	double ke_r_cal = ke_recoil(e_gamma, e_p, e_q, m_e, m_r, th_p, th_q, phi_q);
	double old_diff = ke_r - ke_r_cal;
	//cout << "th_p_deg = " << th_p_deg << " delta = " << delta << endl;
	int N = 0;

	while(1)
		{
		ke_r += step;
		e_q = e_gamma - e_p - ke_r;
		ke_r_cal = ke_recoil(e_gamma, e_p, e_q, m_e, m_r, th_p, th_q, phi_q);
		diff = ke_r - ke_r_cal;
		if(fabs(diff) < tol) break;
		if( diff*old_diff > 0.)
			{
			// same sign
			if(fabs(diff) > fabs(old_diff))
				{
				// we are going in the wrong direction
				step = -step;
				}
			}
		else
			{
			// we passed through the solution
			step = -step/2.;
			}
		old_diff = diff;
		N++;
		}
	//cerr << "N = " << N;
	//cerr << "ke_r = " << ke_r << ", " ;
	return(e_q);
}
double
solve_e_positron2(double e_gamma, double m_e, double m_r,
	double e_p, double th_p, double th_q, double phi)
{
	double e_q[2] = { m_e, e_gamma-e_p};
	double minke = 1.0e-7;
	double ke_r[2];
	ke_r[0] = ke_recoil(e_gamma, e_p, e_q[0], m_e, m_r, th_p, th_q, phi);
	ke_r[1] = ke_recoil(e_gamma, e_p, e_q[1], m_e, m_r, th_p, th_q, phi);

	double val[2];
	val[0] = e_gamma-e_p-e_q[0]-ke_r[0];
	val[1] = e_gamma-e_p-e_q[1]-ke_r[1];
	double val_next = val[0];
	double e_q_next, ke_next;
	int N = 0;

	while( (e_q[1] - e_q[0]) > minke && val_next*val[0]*val[1]!=0. ){
		e_q_next = ( e_q[1]+e_q[0] )/2.;
		ke_next = ke_recoil(e_gamma, e_p, e_q_next, m_e, m_r, th_p, th_q, phi);
		val_next = e_gamma - e_p - e_q_next -ke_next;
		if( val_next*val[0] > 0 && val_next*val[1] < 0){
			e_q[0] = e_q_next;
			ke_r[0] = ke_next;
		}
		else if( val_next*val[0] < 0 && val_next*val[1] > 0){
			e_q[1] = e_q_next;
			ke_r[1] = ke_next;
		}
		else{
			cerr << "solve e_q error." << endl;	
			return -1;
			break;
		}
		N++;	
	}
	// cerr << "N2 = " << N;	
	//cerr << " ke_r 2 = " << ke_r[1] << endl;
	return e_q[1];
}

double ke_recoil(double e_gamma, double e_p, double e_q, double m_e, double m_r,
	double th_p, double th_q, double phi_q)
	{
	//if(e_gamma - e_p - e_q < 0.) {
	//	printf("E_gamma = %.2f E_p = %.2f E_q = %.2f\n",e_gamma, e_p, e_q);
	//	return -1.;
	//	}
	//if(e_p < m_e) return -1.;
	//if(e_q < m_e) return -1.;
	double p_p = sqrt(e_p*e_p - m_e*m_e);
	double p_q = sqrt(e_q*e_q - m_e*m_e);
	double p_p_x = p_p * sin(th_p);
	double p_p_z = p_p * cos(th_p);
	double p_p_y = 0.;
	double p_q_x = p_q * sin(th_q) * cos(phi_q);
	double p_q_y = p_q * sin(th_q) * sin(phi_q);
	double p_q_z = p_q * cos(th_q);
	double p_r_z = e_gamma - p_p_z - p_q_z;
	double p_r_x = -(p_p_x + p_q_x);
	double p_r_y = -(p_p_y + p_q_y);
	double p_r_2 = p_r_x*p_r_x + p_r_y*p_r_y + p_r_z*p_r_z;

//      Phase space issue addressed by Blaine here, the theta, phi of recoil nucleus
//      is still undetermined, may affect final value of cross section/ rate observed by experiment.

//	double p_r = sqrt(p_r_2);
//	printf("Electron Momentum = %f\n", p_p);
//	printf("  (x,y,z) = (%f,%f,%f)\n", p_p_x, p_p_y, p_p_z);
//	printf("Positron Momentum = %f\n", p_q);
//	printf("  (x,y,z) = (%f,%f,%f)\n", p_q_x, p_q_y, p_q_z);
//	printf("Recoil Momentum = %f\n", p_r);
//	printf("  (x,y,z) = (%f,%f,%f)\n", p_r_x, p_r_y, p_r_z);
	double ke_r = sqrt(p_r_2 + m_r*m_r) - m_r;	
	return(ke_r);
	}

double angle_gen(double th_min, double th_max, double& theta, double& phi){
	if(th_min >= th_max){
		cerr << "Wrong Theta min and max value";
		return 0;
	}
	double pi=3.14159265358979323846264338;
	double deg = pi/180.;
	double max_t_theta = 15.0*deg;
	double max_t_phi = 2.2*deg;	
	double min_polar = 4.0*deg;
	double t_theta, t_phi;
	double norm = cos(th_min) - cos(th_max);	
	double x,y,z;

	long starttime = time(NULL);
	do{
		double cos_polar = cos(th_min)-norm*randfloat();
		theta = acos(cos_polar);
		phi = randfloat()*360.0*deg;
		x = sin(theta)*cos(phi);
		y = sin(theta)*sin(phi);
		z = cos_polar;
		t_theta = atan(x/z);
		t_phi = atan(y/z);
	}while( abs(t_theta) > max_t_theta
		|| abs(t_phi) > max_t_phi
		|| theta < min_polar);
}

double randfloat(void)
{
	double r = random()/double(RAND_MAX);
	return r;
}
/*
#include <limits>
double RandGaussian(double sigma)
{
	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	if(sigma == 0.) return(0.);
	generate = !generate;

	if (!generate)
	   return z1 * sigma;

	double u1, u2;
	do
	 {
	   u1 = randfloat();
	   u2 = randfloat();
	 }
	while ( u1 <= epsilon );
	double A = sqrt(-2.0 * log(u1));
	z0 = A * cos(two_pi * u2);
	z1 = A * sin(two_pi * u2);
	return z0 * sigma;
}

vector<string> split(string str, char delimiter) {
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  
  while(getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }
  
  return internal;
}
*/
