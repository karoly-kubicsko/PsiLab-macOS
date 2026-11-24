// DOPI3R.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <locale>
#include <sstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <time.h>
#include <random>

using namespace std;

extern "C" void potv_h216o_cvrqd_(double *vp, double *x1, double *x2, double *x3);
extern "C" void potv_d216o_cvrqd_(double *vp, double *x1, double *x2, double *x3);
extern "C" void potv_h217o_cvrqd_(double *vp, double *x1, double *x2, double *x3);
extern "C" void potv_h218o_cvrqd_(double *vp, double *x1, double *x2, double *x3);
//extern "C" void potv_(double *vp, double *x1, double *x2, double *x3, double *g1, double *g2);
/*
typedef vector<int> Column;
typedef vector<Column> Matrix;
#include "Basic.h"

typedef double DP;
typedef NRVec<int> Vec_INT;
typedef NRVec<double> Vec_DP, VEC_IO_DP, Vec_IO_DP;
typedef NRMat<double> Mat_DP, Mat_IO_DP;
typedef NRMat<int> Mat_INT;
*/

void potv(double &vpot, double x1, double x2, double x3);
void jacobi(int n, double **a, double *d, double **v, int &nrot);
//void eigsrt(double *d, double **v, int n);
void legmat(int np, int km, double **h3, double *q, int np_lmat, ofstream &grid);
void eigsrt2(int n, double *d, double **v);
void jacbolbelso2(double tm3,double tm2,double R1j,double R2j,double xcosj,double &R1b,double &R2b,double &xcosb);
void MMUL(double **A, double **B, double **C, int m, int n, int nm);
void legmat2(int np, int km, int km1, double ajm, double **h4, double **ht3, int np_cmat);
void precond(int n, double *x, double *b ,int *nrow, int *ncol, double *val, double *diag);
void norman(double **r, double &beta, int n, int j);
void norm(double *v, int n, double &orma);
void tqli(double *d, double *e, double **z, int n);
void tqli2(double *d, double *e, double **z, int n);
void tred2(int n, double **a, double *d, double *e);
void SVLprecgm(int n, int nmax, int maxiter, double hiba, double *de, int ndont,int *nrow, int *ncol, double *val, double *diag);

void lanczos(int ns1, int np, int jmax, int prodnDIM, int maxiter, int neig, int bsvektor, double *val, int *ncol, int *nrow, double *diag, double *ud, int *diagind, double lambda_min, double lambda_max, string &eivectype);
void precond_nostore(int prodnDIM, double *z, double *v, double *r, double *val, int *ncol, int *nrow, int *diagind, int *IOPT);
void axpy(int n, double da, double *x, double *y, ofstream &out1, int *IOPT);
void daxpy(int n, double da, double *dx, int incx, double *dy, int incy);
void mult(int, double *x, double *y, double *val, int *ncol, int *nrow);
double norm2(int n, double norma2, double *v, int *IOPT);
double dnrm2(int n, double *x, int incx);
void scal(int n, double da, double *x, int *IOPT);
void dscal(int n, double da, double *dx, int incx);
double proda(int n, double *x, double *y, double prod, int *IOPT);
double ddot(int n, double *dx, int incx, double *dy, int incy);
void eigsrt(double *d, double **v, int n);
void expectval(int nmax, int ns1, int np, int nt, double Re, double Rmax, double tm2atom, double tm3atom, string eivectype);
void expectval_rot(int nmax, int ns1, int np, int jmin, int jmax, int nt, double Re, double Rmax, double tm2atom, double tm3atom, string eivectype);
double power(double *val, int *ncol, int *nrow, int *diagind, int n, double eps, int maxiter, int nonzero, string &eivectype);
void powerfordavidson(double **A, double *val, int *ncol, int *nrow, double *y, double &lambda2, int n, double eps, int maxiter);
void comp2mtx(int nt, int nonzero, double **A, double *val, double *diag, double *ud, int *ncol, int *nrow);
void minres(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *ud, double *x2, double *b, ofstream &out2);
void minresforjacobidavidson(int n, int maxiter, double eps, double theta, double mu1, double *y, double *k, double *u, double *val, int *ncol, int *nrow, int *diagind, double *x2, double *b);
void minresfordavidson(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double **A, double *x2, double *b, int param);
//void minresfordavidson(int n, int maxiter, double eps, double **A, double *x2, double *b);
void CGM(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *ud, int *diagind, double *x2, double *b, ofstream &out1, int precondtype, ofstream &out2);
void CGMfordavidson(int n, int maxiter, double eps, double **A, double *x2, double *b, double theta, ofstream &out1);
void mult2(int n, double *x, double *y, double *val, int *ncol, int *nrow);
void chebyshev(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *ud, double *x2, double *b, double lambda_min, double lambda_max);
void BiCGM(int n, int maxiter, double eps, double *val2, int *ncol, int *nrow, double *ud, double *x2, double *b);
double invpower(double *val, int *ncol, int *nrow, double *ud, int *diagind, int n, double eps, int maxiter, string &eivectype);
void MinvA(int n, double *val, double *val2, int *ncol, int *nrow, double *diag);
void transform(int q, double *val, double *val2, int *ncol, int *nrow, double *diag, int *diagind, double theta);
void transform2(int n, double *diag, double *diagtemp, double theta);
double RQiteration(double *val, double *val2, int *ncol, int *nrow, double *ud, int *diagind, int n, double eps, int maxiter, int q, string &eivectype);

void davidson(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *diag, int *diagind, double *ud, string &eivectype, int vector);
void ILUprecond(int n, double *val, int *ncol, int *nrow, double *ud, int *diagind, double *y, double *x);
void ILUprecondforjacobidavidson(int n, double theta, double *val, int *ncol, int *nrow, double *ud, int *diagind, double *y, double *x);
void ILUprecondfordavidson(int n, double **A, double *x, double *y);
void jacobifordavidson(int n, double **a, double *d, double **v, int &nrot);
void tred2fordavidson(int n, double **a, double *d, double *e);
void tqli2fordavidson(double *d, double *e, double **z, int n);
void eigsrtfordavidson(double *d, double **v, int n);
void multfordavidson1(int, double *x, double *y, double *val, int *ncol, int *nrow);
void multfordavidson2(int, int m, double **X, double **Y, double *val, int *ncol, int *nrow);
void multfordavidson3(int, int m, double **X, double **Y, double *val, int *ncol, int *nrow);
void mult_tuti(int n, double *x, double *y, double *val, int *ncol, int *nrow, int *diagind, double theta);

void psikiind(int nt, int nsp1nt, int ns1, int np, double Re, double Rmax, double tm1atom, double tm2atom, double tm3atom, double *val, double *diag);
double gauss1(double sign, double scal, double gscal, double var, double mu, double sigma, double pot);
double gauss2(double sign, double scal, double gscal, double var, double mu, double sigma, double pot);
void jacbolbelso2_r1(double tm3,double tm2,double R1j,double R2j,double xcosj,double &R1b,double &R2b,double &xcosb);
void jacbolbelso2_r2(double tm3,double tm2,double R1j,double R2j,double xcosj,double &R1b,double &R2b,double &xcosb);
void jacbolbelso2_theta(double tm3,double tm2,double R1j,double R2j,double xcosj,double &R1b,double &R2b,double &xcosb);

void trialbasis(int n, int bsize, double **G, double *val, int *ncol, int *nrow);
void DirectIteration(int n, int nonzero, int maxiter, double *val, int *ncol, int *nrow, double *ud, int *diagind);

void jacobidavidson(int n, double m, int maxiter, double eps, double *val, int *ncol, int *nrow, double *diag, int *diagind, string &eivectype, int vector);
void jdcorrectionequation(int n, int nonzero, int maxiter, double eps, double lambda, double* v, double* r, double *y, double* val, double* ud, int* ncol, int* nrow, int* diagind, ofstream &file);
void jdcorrectionequation_tuti(int n, int nonzero, int maxiter, double eps, double lambda, double* v, double* r, double *y, double* val, double* ud, int* ncol, int* nrow, int* diagind, int step, ofstream &file);

void varhatoertektolelteres(int nmax, int ns1, int np, int jmin, int jmax, int nt, double Re, double Rmax, double tm2atom, double tm3atom, string eivectype);

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dsqrarg;
#define SQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

int main()
{
	int np = 2;
	int ns1 = 2;
	int jmax = 0;
	int maxiter = 1000;
	int nd=1000000,nd1=5000;

	double *de = new double[ns1];
	double *q = new double[np];
	double *DEh = new double[maxiter];
	double *TERM = new double[maxiter];
	double *xmass = new double[3];
	double *r1 = new double[ns1];
	double *r2 = new double[ns1];
	double *r3 = new double[ns1];
	double **hr1 = new double *[ns1];
	for(int i = 0; i < ns1; i++)
		hr1[i] = new double[ns1];
	double **h3 = new double *[np];
	for(int i = 0; i < np; i++)
		h3[i] = new double[np];
	double **h4 = new double *[np];
	for(int i = 0; i < np; i++)
		h4[i] = new double[np];
	double **ht4 = new double *[np];
	for(int i = 0; i < np; i++)
		ht4[i] = new double[np];
	double **tv2 = new double *[ns1];
	for(int i = 0; i < ns1; i++)
		tv2[i] = new double[ns1];
	double **hr2 = new double *[ns1];
	for(int i = 0; i < ns1; i++)
		hr2[i] = new double[ns1];
	
	double *val = new double[nd*3*nd1];
	double *val2 = new double[nd*3*nd1];
	double *diag = new double[nd+1];
	int *ncol = new int[nd*3*nd1];
	int *nrow = new int[nd+1];
	int *diagind = new int[nd+1];
/*	Vec_DP de(ns1),q(np),DEh(maxiter),TERM(maxiter),xmass(3);
	Vec_DP r1(ns1),r2(ns1),r3(ns1);
	Mat_DP hr1(ns1,ns1),h3(np,np),h4(np,np),ht4(np,np);
	Mat_DP tv2(ns1,ns1),hr2(ns1,ns1);

	Vec_DP val(nd*3*nd1),val2(nd*3*nd1),diag(nd+1);
	Vec_INT ncol(nd*3*nd1),nrow(nd+1),diagind(nd+1);	*/

	double amtoau = 1.8228883E3;
	double x1 = 1.0, x0 = 0.00;

	int npar = 1; // 0 = odd (páratlan); 1 = even (páros)

	int nmax = 0;
	if(jmax > 0)
		nmax = 1;
	if(jmax == 0) 
		nmax = 8;
	
	int sv = 1; // 0: nem kell sajátvektor; 1: számol sajátvektort

	double hiba = pow(10.0,-10.0);

	// 1=Jacobi; 2=Radau
	int nkoord = 1;

	string par;

	if(npar == 0) par = "odd";
	if(npar == 1) par = "even";

	int kmax, jmin, ntol;

	if(npar == 0)
	{
		kmax=jmax+1;
		jmin=0;
		ntol=0;
	}
	else
	{
		kmax=jmax;
		jmin=1;
		ntol=0;
	}

	if(jmax == 1)
	{
		kmax=2*jmax+1;
		jmin=-1;
		ntol=1;
		npar=0;
	}

	if(jmax == 0)
	{
		kmax=1;
		jmin=0;
		ntol=0;
		npar=0;
	}

	int nx=np*ns1*ns1;
	int nt=nx*kmax;
	int nsp=(ns1+ns1+np+np-2)*nt;
	int nsp1=ns1+ns1+np-2;
	int nsp1nt = nsp1*nt;
	double  ajm=(double)jmax;

	if(jmax != 0)
	{
		cout << "Nonzero Elements of the Hamiltonian Matrix: nsp = "<< nsp << endl;
	}
	else
	{
		cout << "Nonzero Elements of the Hamiltonian Matrix: nsp1nt = "<< nsp1nt << endl;
	}
	
	int h2o_type;
	cout << "H2(16)O: 1 || H2(17)O: 2 || D2(16)O: 3 || H2(18)O: 4" << endl;
	h2o_type = 1;
	if(h2o_type == 1)
	{
		xmass[0] = 1.00727647000000;
		xmass[1] = 1.00727647000000;
//		xmass[1] = 2.01355321400000;
//		xmass[2] = 17.9947714000000;
		xmass[2] = 15.9905260000000;
	}
	else if(h2o_type == 2) 
	{
		xmass[0] = 1.00727647000000;
		xmass[1] = 1.00727647000000;
		xmass[2] = 16.9947425000000;
	}
	else if(h2o_type == 3)
	{
		xmass[0] = 2.01355321400000;
		xmass[1] = 2.01355321400000;
		xmass[2] = 15.9905260000000;
	}
	else if(h2o_type == 4)
	{
		xmass[0] = 1.00727647000000;
		xmass[1] = 1.00727647000000;
		xmass[2] = 17.9991610000000;
	}
	for(int i = 0; i < 3; i++)
		cout << xmass[i] << endl;
//	m_H = 1.00727647000000
//	m_D = 2.01355321400000
//	m_16O = 15.9905260000000
//	m_17O = 16.9947425000000
//	m_18O = 17.9991610;
	
// HD16O H218O HD18O kellenek!!!

	double tm1atom=xmass[0]*amtoau;
	double tm2atom=xmass[1]*amtoau;
	double tm3atom=xmass[2]*amtoau;

	double g1, g2;

	if(nkoord ==1) {
         g1 = xmass[1] / (xmass[1] + xmass[2]);
         g2 = x0;
	}

	if(nkoord == 2) {
	//        radau coordinates
         double a = sqrt(xmass[2] / (xmass[0]+xmass[1]+xmass[2]));
         double b = xmass[1] / (xmass[0]+xmass[1]);
         g1 = x1 - a / (a+b-a*b);
         g2 = x1 - a / (x1-b+a*b);
	}

	double tred1 = amtoau/(g2*g2/xmass[0]+x1/xmass[1]+(x1-g2)*(x1-g2)/xmass[2]);
	double tred2 = amtoau/(x1/xmass[0]+g1*g1/xmass[1]+(x1-g1)*(x1-g1)/xmass[2]);

	//cout<<tred1<<"  "<<tred2<<endl; 

//*****************************************************************
//				   HERMITE
//*****************************************************************


//******A -d2/dR1**2 matrixelemei DVR-ben*******************

	for(int i = 0; i < ns1; i++ ) {
		for(int j = 0; j < ns1; j++ ) {
			hr1[i][j] =0.00;
		}
	}

	for(int i = 0; i < ns1-1; i++ ) {
		hr1[i][i+1] = sqrt(((double)i+1.0)/2.0);
	}
	
	for(int i = 1; i < ns1; i++ ) {
		hr1[i][i-1] = sqrt(((double)i)/2.0);
		//	valamiért anno erre atirtam: hr1[i][i-1] = sqrt(((double)i)/2.0);
	}

	int nrot = 100;
	jacobi(ns1, hr1, de, tv2, nrot);
	eigsrt(de, tv2, ns1);

	double Re=2.000;
	double Rmax=0.800;
//	1.8200
//	0.7
	double allk = de[ns1-1]/Rmax;
	double coef = allk*allk/(2.00*tred1);
	
	ofstream grid;
	grid.open ("gridfile.txt");
	grid.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < ns1; i++)
	{
		grid << setprecision(30) << de[i] << endl;
	}

	for(int i = 0; i < ns1; i++ ) {
		for(int j = 0; j < ns1; j++ ) {
			hr1[i][j] =0.00;
		}
	}

	for(int j = 0; j < ns1; j++)
	{
		for(int k = 0; k < ns1; k++)
		{
			if(j == k)
			{
				hr1[j][j]=coef*(2.00*(ns1-1.0)/3.0+0.500-(de[j]*de[j])/3.00);
			} else {
				hr1[j][k]=(-0.500+2.00/((de[j]-de[k])*(de[j]-de[k])))*coef;
			}
		}
	}

	

//******A -d2/dR2**2 matrixelemei DVR-ben*******************

	if(nkoord == 1) {
		for(int i = 0; i < ns1; i++ ) {
			for(int j = 0; j < ns1; j++ ) {
				hr2[i][j] =0.00;
			}
		}

		coef=allk*allk/(2.00*tred2);

		for(int j = 0; j < ns1; j++)
		{
			for(int k = 0; k < ns1; k++)
			{
				if(j == k)
				{
					hr2[j][j]=coef*(2.00*(ns1-1.0)/3.0+0.500-(de[j]*de[j])/3.00);
				} else {
					hr2[j][k]=(-0.500+2.00/((de[j]-de[k])*(de[j]-de[k])))*coef;
				}
			}
		}
	} else
	{
		for(int i = 0; i < ns1; i++ ) {
			for(int j = 0; j < ns1; j++ ) {
				hr2[i][j] =hr1[i][j];
			}
		}
	}
/*	for(int i = 0; i < ns1; i++)
	{
		for(int j = 0; j < ns1; j++)
		{
			cout << setprecision(10) << hr1[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;
	for(int i = 0; i < ns1; i++)
	{
		for(int j = 0; j < ns1; j++)
		{
			cout << setprecision(10) << hr2[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;	*/

//******Az (1/R1**2) matrixelemei DVR-ben*********************

	for(int i= 0; i < ns1; i++)
	{
		//r1(i)=(1.0d0/(2.0d0*tred1))/(de(i)/allk+RE)**2
		r1[i] = (1.00/(2.00*tred1))/((de[i]/allk+Re)*(de[i]/allk+Re));
//		cout << setprecision(10) << r1[i] << " ";
	}

//******Az (1/R2**2) matrix elemei DVR-ben*********************
    
	for(int i= 0; i < ns1; i++)
	{
		//r2(i)=(1.0d0/(2*tred2))/(de(i)/allk+RE)**2
		r2[i] = (1.00/(2.00*tred2))/((de[i]/allk+Re)*(de[i]/allk+Re));
//		cout << setprecision(10) << r2[i] << " ";
	}
     
	int ki=0;
	int kiv=0;
	int kif1=0;
	int kif2=0;
	int kcik=0;
	int id= 0;
	 
//******Az (1/R3**2) matrix elemei DVR-ben*********************
	
	ofstream potencial;
	potencial.open ("POTENCIAL.txt");
	potencial.setf(ios_base::fixed,ios_base::floatfield);

	double akm, var1, var2, var3, Vpot, potval = 1000.0;
	int kmm, kmm1, nxx, z = 0;
	for(int km = jmin; km <= jmax; km++)
	{
		akm = (double)km;
		kmm= abs(km);
		for(int i = 0; i < ns1; i++)
		{
			r3[i]=(ajm*(ajm+1.00)-2.00*akm*akm)*r1[i];
		}
		legmat(np,kmm,h3,q,0,grid);
		int kc=0;
		for(int km1 = jmin; km1 <= jmax; km1++)
		{
			kmm1 = abs(km1);
			if(km == km1)
			{
				for(int i = 0; i < ns1; i++)
				{
					for(int j = 0; j < ns1; j++)
					{
						for(int l = 0; l < np; l++)
						{
							nxx=(km-npar+ntol)*nx;
							if((km-npar+ntol) != 0) 
								ki=ki+np;
							for(int iv = 0; iv < ns1; iv++)
							{
								if(iv == i)
								{
									for(int jv = 0; jv < ns1; jv++)
									{
										if(jv == j)
										{
											for(int lv = 0; lv < np; lv++)
											{
												//ncol(ki)=(lv+(j-1)*np+(i-1)*ns1*np)+nxx
												ncol[ki]=(lv+(j)*np+(i)*ns1*np)+nxx;
												if(lv == l)
												{
													val[ki]=hr1[i][i]+hr2[j][j]+h3[l][l]*r1[i]+h3[l][l]*r2[j];
													// cout<<hr1[i][i]<<"  "<<hr2[j][j]<<"  "<<h3[l][l]<<"  "<<r1[i]<<"  "<<r2[j]<<endl;
													jacbolbelso2(tm3atom,tm2atom,de[i]/allk+Re,de[j]/allk+Re,q[l],var1,var2,var3);
//													potv(Vpot, var1, var2, var3);
//													Vpot = 0.00;
													if(h2o_type == 1)
													{
														potv_h216o_cvrqd_(&Vpot, &var1, &var2, &var3);
//														if(Vpot < potval)
//															potval = Vpot;
														Vpot = Vpot - 0.00358028360318796825754850132739;
													}
													if(h2o_type == 2)
													{
														potv_h217o_cvrqd_(&Vpot, &var1, &var2, &var3);
//														if(Vpot < potval)
//															potval = Vpot;
														Vpot = Vpot - 0.00343991010054164597375307543814;
													}
													if(h2o_type == 3)
													{
														potv_d216o_cvrqd_(&Vpot, &var1, &var2, &var3);
//														if(Vpot < potval)
//															potval = Vpot;
														Vpot = Vpot - 0.00337543999746865977740843511867;
													}
													if(h2o_type == 4)
													{
//														potv_h218o_cvrqd_(&Vpot, &var1, &var2, &var3);
//														if(Vpot < potval)
//															potval = Vpot;
//														Vpot = Vpot - 0.00337543999746865977740843511867;
													}
//													potv_(&Vpot, &var1, &var2, &var3, &g1, &g2);
//													cout << "Vpot " << Vpot * 219474.630 << endl;
//													if(Vpot * 219474.630 >= 438949.26)
//														Vpot = (438949.26 / 219474.630);
													potencial << setprecision(30) << Vpot << endl;
													if(jmax > 0)
													{
														val[ki]=val[ki]+r3[i]+Vpot; 
													} 
													else
													{
														val[ki]=val[ki]+Vpot;
													}
													diag[kiv]=val[ki];
													diagind[z] = ki;
													z++;
													kiv=kiv+1;
												} 
												else
												{   
													val[ki]= h3[l][lv]*(r1[i]+r2[j]);
												}
												ki=ki+1;
											} // endfor lv
										}
										else
										{  
											ncol[ki]=((jv)*np+l+(i)*ns1*np)+nxx;  
											val[ki]=hr2[j][jv];
											ki=ki+1;
										}
									} // endfor jv
								} 
								else
								{    
									ncol[ki]=((iv)*ns1*np+l+(j)*np)+nxx;  
									val[ki]=hr1[i][iv];
									ki=ki+1;
								}
							}// endfor iv
							if(km!= jmax) 
								ki=ki+np;
						} //endfor l
					}//endfor j
				}// endfor i
			} // if km == km1 

			if(jmax > 0)
			{
				 if(km1 == km+1) 
				 {
					 nxx=(km1-npar+ntol)*nx;
					 legmat2(np,kmm,kmm1,ajm,h4,ht4,0);

					 for(int i =0; i < ns1; i++)
					 {
						 for(int j =0; j < ns1; j++)
						 {
							 for(int l = 0; l < np; l++)
							 {
								 if((km+ntol-npar) == 0) 
									 kif1=kif1+nsp1;

								 if(((km+ntol-npar) > 0) && ((km-npar) < jmax)) 
									 kif1=kif1+nsp1+np;

								 for(int iv = 0; iv < ns1; iv++)
								 {
									 if(iv == i)
									 {
										 for(int jv = 0; jv < ns1; jv++)
										 {
											 if(jv == j)
											 {
												 for(int lv=0; lv<np; lv++)
												 {
													ncol[kif1]=(lv+j*np+i*ns1*np)+nxx;
													
													 if(lv == l)
													{
														val[kif1] = h4[l][l]*r1[i];

													} else 
													{
														val[kif1]=h4[l][lv]*r1[i];
													}// lv == l
													kif1=kif1+1;
												 }//endfor lv
											 }// jv == j
										 } // endfor jv
									 } // iv = i
								 }//endfor iv

							 } //endofr l
						 } // endfor j
					 } // endfor i
				 } //endif km1 == km+1

				 if(km1 == km-1) 
				 {
					  if(kcik == 0)
					  {
						kif2=(nsp1+np)*nx;
						kcik=1;
					  }

					  nxx=(km1-npar+ntol)*nx;
					  legmat2(np,kmm1,kmm,ajm,h4,ht4,0);

					   for(int i =0; i < ns1; i++)
						 {
							 for(int j =0; j < ns1; j++)
							 {
								 for(int l = 0; l < np; l++)
								 {
									 for(int iv = 0; iv < ns1; iv++)
									 {
										 if(iv == i)
										 {
											 for(int jv = 0; jv < ns1; jv++)
											 {
												 if(jv == j)
												 {
													 for(int lv=0; lv<np; lv++)
													 {
														 //ncol(kif2)=(lv+(j-1)*np+(i-1)*ns1*np)+nxx
														ncol[kif2]=(lv+(j)*np+(i)*ns1*np)+nxx;
														
														 if(lv == l)
														{
															val[kif2]=ht4[l][l]*r1[i];

														} else 
														{
															 val[kif2]=ht4[l][lv]*r1[i];
														}// lv == l
														kif2=kif2+1;
													 }//endfor lv
												 }// jv == j
											 } // endfor jv
										 } // iv = i
									 }//endfor iv

									 if((km) < jmax) 
									 {
										kif2=kif2+nsp1+np;
									 }
									 else
									 {
										kif2=kif2+nsp1;
									 }

								 } //endofr l
							 } // endfor j
						 } // endfor i

				 }//endif km1 == km-1

			 } // endif jmax > 0



		} // endfor km1
		
	} // endfor km
	
	potencial.close();

//*********** create nrow vector *************

	if(jmax > 0)
	{
		nrow[0]=0;

		for(int i = 0; i < nx; i++)
		{
			nrow[i+1]= nrow[i]+(2*ns1+2*np-2);
		}

		for(int i = nx; i < nt-nx; i++)
		{
			nrow[i+1]= nrow[i]+(2*ns1+3*np-2);
		}

		for(int i = nt-nx; i < nt; i++)
		{
			nrow[i+1]= nrow[i]+(2*ns1+2*np-2);
		}
		
	} else 
	{
		
		nrow[0]=0;

		for(int i = 0; i < nx; i++)
		{
			nrow[i+1]= nrow[i]+ns1+ns1+np-2;
		}
	}
	
//	SVLprecgm(nt,nmax,maxiter,hiba ,DEh,0, nrow, ncol, val, diag);

//	!!!!!!!!!!!!!!új részek!!!!!!!!!!!!!!
	
	double *ud = new double[nd+1];
//	Vec_DP ud(nd+1);
	for(int i = 0; i < nt; i++)
	{
		if(abs(diag[i]) >= 1.0E-10)
		{
			ud[i] = 1.0 / diag[i];
		}
		else
		{
			ud[i] = 1.0;
//			ud[i] = 1.0E10;
		}
	}
	
	cout << "nt = " << nt << endl;
	
//	MinvA(nt,val,val2,ncol,nrow,diag);
	double **A = new double *[nt];
	for(int i = 0; i < nt; i++)
		A[i] = new double[nt];
//	Mat_DP A(nt,nt);
	int szam = 0;
	for(int i = 0; i < nt*nt; i++)
	{
		if(abs(val[i]) != 0.000000)
			szam++;
//		cout << "val[" << i << "] = " << setprecision(30) << val[i] << endl;
	}
	cout << "szam " << szam << endl;	
	comp2mtx(nt,szam,A,val,diag,ud,ncol,nrow);
	
/*	ofstream szamok;
	szamok.open ("szamok.txt");
	szamok.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < nt; i++)
		szamok << i + 1 << endl;
	szamok.close();	*/
//	cout << setprecision(30) << potval << endl;	
	
	string eivectype = " ";
	double lambda_min, lambda_max;
//	varhatoertektolelteres(nmax,ns1,np,jmin,jmax,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
	psikiind(nt,nsp1nt,ns1,np,Re,Rmax,tm1atom,tm2atom,tm3atom,val,diag);
//	trialbasis(nt,bsize,G,val,ncol,nrow);
//	lambda_min = invpower(val,ncol,nrow,ud,diagind,nt,hiba,maxiter,eivectype);
//	expectval_rot(nmax,ns1,np,jmin,jmax,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
//	RQiteration(val,val2,ncol,nrow,ud,diagind,nt,hiba,maxiter,nsp1nt,eivectype);
//	expectval_rot(nmax,ns1,np,jmin,jmax,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
//	lambda_max = power(val,ncol,nrow,diagind,nt,hiba,10000000,nsp1nt,eivectype);
//	expectval_rot(nmax,ns1,np,jmin,jmax,nt,Re,Rmax,tm2atom,tm3atom,eivectype);

	lanczos(ns1,np,jmax,nt,maxiter,nmax,sv,val,ncol,nrow,diag,ud,diagind,lambda_min,lambda_max,eivectype);
//	expectval(nmax,ns1,np,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
	expectval_rot(nmax,ns1,np,jmin,jmax,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
//	davidson(nt,maxiter,hiba,val,ncol,nrow,diag,diagind,ud,eivectype,sv);
//	expectval(nmax,ns1,np,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
//	expectval_rot(nmax,ns1,np,jmin,jmax,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
	jacobidavidson(nt,nsp1nt,maxiter,hiba,val,ncol,nrow,diag,diagind,eivectype,sv);
//	expectval(nmax,ns1,np,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
	expectval_rot(nmax,ns1,np,jmin,jmax,nt,Re,Rmax,tm2atom,tm3atom,eivectype);
	
	cout << "DONE!" << endl;
	return 0;
}

// új fgv-k

void lanczos(int ns1, int np, int jmax, int prodnDIM, int maxiter, int neig, int bsvektor, double *val, int *ncol, int *nrow, double *diag, double *ud, int *diagind, double lambda_min, double lambda_max, string &eivectype)
{
	clock_t t;
	t = clock();
	double *q2help = new double[prodnDIM];
	double *qtmp = new double[prodnDIM];
	double **lanvec = new double *[prodnDIM];
	for(int i = 0; i < prodnDIM; i++)
		lanvec[i] = new double[maxiter];
	double **tv = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		tv[i] = new double[maxiter];
	double **svmat = new double *[prodnDIM];
	for(int i = 0; i < prodnDIM; i++)
		svmat[i] = new double[neig];
	double *alpha = new double[maxiter];
	double *beta = new double[maxiter];
	double *elozo = new double[maxiter];
	double *de = new double[maxiter];
	double *ee = new double[maxiter];
	double *diff = new double[maxiter];
	double *er = new double[maxiter];
	double *temp = new double[prodnDIM];
	double *temp2 = new double[prodnDIM];
/*	Vec_DP q2help(prodnDIM), qtmp(prodnDIM);
	Mat_DP lanvec(prodnDIM,maxiter), tv(maxiter,maxiter), svmat(prodnDIM,neig);
	Vec_DP alpha(maxiter),beta(maxiter),elozo(maxiter),de(maxiter),ee(maxiter);
	Vec_DP diff(maxiter),er(maxiter);	*/
	
	double alpha1,beta1,gamma1,hlpp,hiba,vhiba;
	int nmax,nstep,nconv,n,nkezd,nstop,i,j,k,l;
	int iter,iterszam,ncheckeigvalues;
	int mmem,ihlpp;
	string svec;
	
	int *IOPT = new int[100];
//	Vec_INT IOPT(100);
	
	ofstream out1;
	out1.open ("lanczos.txt");
	out1.setf(ios_base::fixed,ios_base::floatfield);
	ofstream out2;
//	out2.open ("lanczoslinearsolveriter.txt");
	out2.setf(ios_base::fixed,ios_base::floatfield);
	
	const double ehcm = 219474.630; // hartree -> cm-1
	
//	IOPT[74] = 2;			// reorthogonalization,0:partial,1:local,2:periodic
	IOPT[96] = 10;		// 10E-IOPT[96] !convergence cond for eigenvalues
	IOPT[69] = 10;		// 10E-IOPT[69] !convergence cond for eigenvectors
	IOPT[64] = 12; 		// 10E-IOPT[64] !convergence cond in CGM
	IOPT[66] = 12;		// 10E-IOPT[66] !convergence cond in MINRES
	IOPT[0] = 11;			// 10E-IOPT[0] !convergence cond in Chebyshev
//	IOPT[76] = 10000; 		// LAN1.DAT
//	IOPT[77] = 0; 		// 0 CG precond, 1 [H-EI]^2  , 2 [L_max - [H-EI]^2]^[2n+1]
//	IOPT[78] = 1;
//	IOPT[62] = 262144;		// record length in direct access files, 8192/262144
	IOPT[60] = 0;			// UJ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IOPT[61] = 250000000;	// how big vector to use openmp not blas       
//	IOPT[63] = IOPT[62]*2;	// how big record
//	IOPT[68] = 0;			// test of ortogonality
	IOPT[65] = 10000;		// precond iteration limit
	IOPT[67] = 10000;		// minres iteration limit
	IOPT[68] = 10000;		// cgm iteration limit
	IOPT[1] = 10000;		// Chebyshev iteration limit
	IOPT[70] = 1;			// 1 v pameti;  0 na disku  DIAG, X
	IOPT[71] = 1;			// seq job
	IOPT[2] = 3;			// 1: nincs prekondícionálás, 2: Jacobi-prekondícionálás, 3: ILU-prekondícionálás
	
	mmem = 1;	// scr{mem} is {z,x,dia2,yy,vpot}
	if(IOPT[71] == 1)
	{
		mmem = mmem + 2;	// sequence run in mult
		if(IOPT[70] == 1)
		{
			mmem = mmem + 2;	// Diag,X in memory
		}
	}
	
	double *scr = new double[mmem*prodnDIM];
//	Vec_DP scr(mmem*prodnDIM);
	
	for(i = 0; i < prodnDIM; i++)
	{
		scr[i + 2*prodnDIM] = ud[i];
	}
	
	cout << "kezdődik az új Lánczos" << endl;
//	INITIALIZATION OF VARIABLES AND ARRAYS

	for(i = 0; i < maxiter; i++)
	{
		alpha[i] = 0.0;
		beta[i] = 0.0;
		elozo[i] = 0.0;
	}

	for(i = 0; i < prodnDIM; i++)
	{
		for(j = 0; j < maxiter; j++)
		{
			lanvec[i][j] = 0.0;
		}
	}
	
	if(bsvektor == 1) svec = "Yes";
	if(bsvektor == 0) svec = "No";

	out1 << "This is the Lánczos output file             " << endl;
	out1 << "Basis number for Streching:                 " << ns1 << endl;
	out1 << "Basis number for Bending:                   " << np << endl;
	out1 << "Rotational Quantum number:                  " << jmax << endl;
	out1 << "Matrix dimension:                           " << prodnDIM << endl;
	out1 << "Maximum number of iterations:               " << maxiter << endl;
	out1 << "Number of wanted eigenvalues:               " << neig << endl;
	out1 << "Calculating eigenvectors:                   " << svec << endl;
	out1 << "Convergence condition for Eigenvalues:      " << "10.0E" << -IOPT[96] << endl;
	out1 << "Convergence condition for Eigenvectors:     " << "10.0E" << -IOPT[69] << endl;
	out1 << "Convergence condition for preCGM:           " << "10.0E" << -IOPT[64] << endl;
	
//	FIRST LANCZOS VECTOR
	
//	qtmp = diag; 
	for(i = 0; i < prodnDIM; i++)
	{
		qtmp[i] = diag[i];
	}
	
	double *q1 = new double[prodnDIM];
	double *q2 = new double[prodnDIM];
//	Vec_DP q1(prodnDIM), q2(prodnDIM);
	
//	q1 => qtmp(1:prodnDIM)
	for(i = 0; i < prodnDIM; i++)
	{
		q1[i] = qtmp[i];
	}
	
//	q2 => lanvec(1:prodnDIM,1)
	for(i = 0; i < prodnDIM; i++)
	{
		q2[i] = lanvec[i][0];
	}
	
/*	cout << "előtt q1 = " << endl;
	for(i = 0; i < prodnDIM; i++)
	{
		cout << q1[i] << endl;
	}
	cout << "előtt q2 = " << endl;
	for(i = 0; i < prodnDIM; i++)
	{
		cout << q2[i] << endl;
	}
	cout << "előtt scr = " << endl;
	for(i = 0; i < mmem*prodnDIM; i++)
	{
		cout << scr[i] << endl;
	}	*/
		
/*	for(i = 0; i < prodnDIM; i++)
	{
		q2[i] = 1.0;
	}	*/
	
//	precond_nostore(prodnDIM,scr,q2,q1,val,ncol,nrow,diagind,IOPT);
//	BiCGM(prodnDIM,IOPT[65],pow(10.0,(-(double)IOPT[64])),val,ncol,nrow,ud,q2,q1);
//	minres(prodnDIM,IOPT[67],pow(10.0,(-(double)IOPT[66])),val,ncol,nrow,ud,q2,q1,out2);
//	CGM(prodnDIM,IOPT[65],pow(10.0,(-(double)IOPT[64])),val,ncol,nrow,ud,diagind,q2,q1,out1,IOPT[2],out2);
//	chebyshev(prodnDIM,IOPT[1],pow(10.0,(-(double)IOPT[0])),val,ncol,nrow,ud,q2,q1,lambda_min,lambda_max);
	
/*	cout << "után q1 = " << endl;
	for(i = 0; i < prodnDIM; i++)
	{
		cout << q1[i] << endl;
	}
	cout << "után q2 = " << endl;
	for(i = 0; i < prodnDIM; i++)
	{
		cout << q2[i] << endl;
	}
	cout << "után scr = " << endl;
	for(i = 0; i < mmem*prodnDIM; i++)
	{
		cout << scr[i] << endl;
	}	*/
	
/*	beta1 = 0.0;
	beta1 = norm2(prodnDIM,beta1,q2,IOPT);
	
	if(beta1 == 0.0)
	{
		cout<< "beta1 = " << beta1 <<endl;
		return;
	}
	beta1 = 1.0 / beta1;
	scal(prodnDIM,beta1,q2,IOPT);
	
//	lanvec feltöltése!!!!!!!!!!!!
	for(i = 0; i < prodnDIM; i++)
	{
		qtmp[i] = q1[i];
	}
	for(i = 0; i < prodnDIM; i++)
	{
		lanvec[i][0] = q2[i];
	}	*/
	ifstream wf;
	wf.open("thing.dat");
	wf.setf(ios_base::fixed,ios_base::floatfield);
	for(i = 0; i < prodnDIM; i++)
	{
		wf >> lanvec[i][0];
	}
	wf.close();
	
	if(maxiter <= (neig+25))
	{
		out1 << "Error occured!!" << endl;
		out1 << "Maxiter must be bigger, gt.neig+25" << endl;
		return;
	}
	
//	PARAMETERS FOR THE FIRST LANCZOS CORE
	
	nkezd = 0;
	nstop = neig;
	nstep = 1;
	iterszam = 1;
	ncheckeigvalues = 1;
A:	;

//	LANCZOS CORE
	
	for(j = nkezd; j < nstop; j++)
	{
		out1 << "*********" << nstep << ". main cycle *********" << endl;
		out1 <<  "*********" << j+1 << ". sub cycle *********" << endl;
		out1 << "*********" << nstop << " is the next stop *********" << endl;
		
		iterszam = iterszam + 1;
		
		for(i = 0; i < prodnDIM; i++)
		{
			q1[i] = lanvec[i][j+1];
		}
		for(i = 0; i < prodnDIM; i++)
		{
			q2[i] = lanvec[i][j];
		}
		for(iter = 0; iter < prodnDIM; iter++)
		{
			q2help[iter] = q2[iter];
		}
		
//		precond_nostore(prodnDIM,scr,q1,q2help,val,ncol,nrow,diagind,IOPT);
//		BiCGM(prodnDIM,IOPT[65],pow(10.0,(-(double)IOPT[64])),val,ncol,nrow,ud,q1,q2help);
//		minres(prodnDIM,IOPT[67],pow(10.0,(-(double)IOPT[66])),val,ncol,nrow,ud,q1,q2help,out2);
		CGM(prodnDIM,IOPT[65],pow(10.0,(-(double)IOPT[64])),val,ncol,nrow,ud,diagind,q1,q2help,out1,IOPT[2],out2);
//		chebyshev(prodnDIM,IOPT[1],pow(10.0,(-(double)IOPT[0])),val,ncol,nrow,ud,q1,q2help,lambda_min,lambda_max);
//		ILUprecond(prodnDIM,val,ncol,nrow,ud,diagind,q2help,q1);
		
		alpha1 = 0.0;
		
		alpha1 = proda(prodnDIM,q1,q2,alpha1,IOPT);
		
		alpha[j] = alpha1;
		out1 << "alpha j " << alpha[j] << j << endl;
		
		if(j > (nkezd+1))
		{
			alpha1 = -1.0 * alpha[j];
			axpy(prodnDIM,alpha1,q2,q1,out1,IOPT);
//			q1[j+1] = q1[j+1] - alpha[j] * q2[j]
			beta1 = -1.0 * beta[j-1];
			for(i = 0; i < prodnDIM; i++)
			{
				q2[i] = lanvec[i][j-1];
			}
			axpy(prodnDIM,beta1,q2,q1,out1,IOPT);
//			q1[j+1] = q1[j+1] - beta[j-1] * q2[j-1]
			
		}
		else
		{
			alpha1 = -1.0 * alpha[j];
			axpy(prodnDIM,alpha1,q2,q1,out1,IOPT);
//			q1[j+1] = q1[j+1] - alpha[j] * q2[j]
		}
		
//		FULL RE-ORTHOGONALIZATION
		for(i = 0; i <= j; i++)
		{
			gamma1 = 0.00;
			for(iter = 0; iter < prodnDIM; iter++)
			{
				q2[iter] = lanvec[iter][j-i];
			}
			for(iter = 0; iter < prodnDIM; iter++)
			{
				gamma1 = gamma1 + q1[iter] * q2[iter];
			}
			gamma1 = -1.0 * gamma1;
			for(iter = 0; iter < prodnDIM; iter++)
			{
				q1[iter] = q1[iter] + gamma1 * q2[iter];
			}
		}
		// MÁS ORTOGONALIZÁCIÓ
/*		for(k = 0; k < j+1; k++)
		{
			temp[k] = 0.0;
			for(i = 0; i < prodnDIM; i++)
			{
				temp[k] = temp[k] + lanvec[i][k] * q1[i];
			}
		}
		for(i = 0; i < prodnDIM; i++)
		{
			temp2[i] = 0.0;
			for(k = 0; k < j+1; k++)
			{
				temp2[i] = temp2[i] + lanvec[i][k] * temp[k];
			}
		}
		for(i = 0; i < prodnDIM; i++)
		{
			q1[i] = q1[i] - temp2[i];
		}	*/
		
//		CALCULATING beta[j]
		
		beta1 = 0.0;
		beta1 = norm2(prodnDIM,beta1,q1,IOPT);
		
		beta[j] = beta1;
		out1 << "beta, j " << beta1 << j << endl;
		if(beta1 == 0.0)
		{
			out1 << "error2 dnrm2 == 0" << endl;
			return;
		}
		beta1 = 1.0 / beta1;
		scal(prodnDIM,beta1,q1,IOPT);
		
		for(i = 0; i < prodnDIM; i++)
		{
			lanvec[i][j+1] = q1[i];
		}
		for(i = 0; i < prodnDIM; i++)
		{
			q2[i] = lanvec[i][j+1];
		}
	} // endfor j
	
//	if ncheckeigvalues is 1 this mean that we need eigenvectors and 
//	the requested eigenvalues are converged so we dont not need to check them
	
	if(ncheckeigvalues != 1)
	{
//		construct the T tridiagonal matrix to check the eigenvalues
		hiba = pow(10.0,(-(double)IOPT[96]));
		for(i = 0; i < nstop; i++)
		{
			de[i] = alpha[i];
		}
		for(i = 0; i < nstop-1; i++)
		{
			ee[i] = beta[i];
		}
		
		for(i = 0; i < neig; i++)
		{
			for(k = 0; k < neig; k++)
			{
				if(i == k)
				{
					tv[i][k]=1.00;
				}
				else
				{
					tv[i][k] = 0.00;
				}
			}
		}
		
		tqli2(de,ee,tv,nstop);
		eigsrt(de,tv,nstop);
		
		for(i = 0; i < nstop; i++)
		{
			diff[i] = 0.0;
			de[i] = 1.0 / (de[i]);
			diff[i] = abs(elozo[i]-de[i]);
			elozo[i] = de[i];
		}
		out1 << endl;
		out1 << "Eigenvalues:" << endl;
		out1 << endl;
		out1 << "th       values / hartree       error / hartree        values / cm-1        error / cm-1" << endl;
		ihlpp = min(nstop,neig);
		
		for (i = 0; i < ihlpp; i++)
		{
			out1 << i+1 << "        " << de[i] << "               " << diff[i] << "               " << de[i] * ehcm << "         " << diff[i] * ehcm << endl;
		}
		nconv = 0;
		for(i = 0; i < neig; i++)
		{
			if(diff[i] > hiba)
			{
				nconv++;
			}
		}
		if(nconv > 0)
		{
			out1 << endl;
			out1 << "Eigenvalues are still wrong :), so again!" << endl;
			nkezd = nstop;
			nstop = nstop + nconv;
			if(maxiter <= nstop)
			{
				out1 << "Error occured!!" << endl;
				out1 << "maxiter must be bigger, than " << nstop << endl;
				return;
			}
			nstep++;
			goto A;
		}
		else
		{
			out1 << endl;
			out1 << "Eigenvalues converged!!" << endl;
			out1 << endl;
			ncheckeigvalues = 1;
		}
	}
// 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// 	Calculate and check eigenvectors
// 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(bsvektor == 1)
	{
//	Ritz vector production and
//	checking of vector error, should be add IOPT for this
		out1 << endl;
		vhiba = pow(10.0,(-(double)IOPT[69]));
		for(i = 0; i < maxiter; i++)
		{
			for(j = 0; j < maxiter; j++)
			{
				tv[i][j] = 0.0;
			}
		}
		for(i = 0; i < nstop; i++)
		{
			de[i] = alpha[i];
		}
		for(i = 0; i < nstop-1; i++)
		{
			ee[i] = beta[i];
		}

		for(i = 0; i < maxiter; i++)
		{
			for(k = 0; k < maxiter; k++)
			{
				if(i == k)
				{
					tv[i][k] = 1.00;
				}
				else
				{
					tv[i][k] = 0.00;
				}
			}
		}
	
		tqli2(de,ee,tv,nstop);
		eigsrt(de,tv,nstop);
	
/*		for(i = 0; i < nstop; i++)
		{
			de[i] = 1.0 / de[i];
		}	*/
		for(i = 0; i < nstop; i++)
		{
			diff[i] = 0.0;
			de[i] = 1.0 / (de[i]);
			diff[i] = abs(elozo[i]-de[i]);
			elozo[i] = de[i];
		}
		out1 << endl;
		out1 << "Eigenvalues:" << endl;
		out1 << endl;
		out1 << "th       values / hartree       error / hartree        values / cm-1        error / cm-1" << endl;
		ihlpp = min(nstop,neig);
		
		for (i = 0; i < ihlpp; i++)
		{
			out1 <<setprecision(8)<< i+1 << "        " << de[i] << "               " << diff[i] << "               " << de[i] * ehcm << "         " << diff[i] * ehcm << endl;
		}
		
		for(iter = 0; iter < neig; iter++)
		{
			for(i = 0; i < prodnDIM; i++)
			{
				q2help[i] = 0.0;
				for(k = 0; k < nstop; k++)
				{
					q2help[i] = q2help[i] + lanvec[i][k] * tv[k][iter];
					svmat[i][iter] = q2help[i];
				}
			}
			mult(prodnDIM,q2help,q1,val,ncol,nrow);
			hlpp = -1.0 * de[iter];
			axpy(prodnDIM,hlpp,q2help,q1,out1,IOPT);
//			q1 = q1 - hlpp * q2help = A * q2help - de * q2help
			beta1 = norm2(prodnDIM,beta1,q1,IOPT);
			er[iter] = beta1;
		}
		out1 << endl;
		out1 << "Error of eigenvectors" << endl;
		out1 << endl;
		for(l = 0; l < neig; l++)
		{
			out1 << scientific << er[l] << endl;
			cout << scientific << er[l] << endl;
		}
		nconv = 0;
		for(i = 0; i < neig; i++)
		{
			if(abs(er[i]) > vhiba)
				nconv++;
		}
		if(nconv > 0)
		{
			out1 << endl;
			out1 << "Eigenvectors are still wrong :), so again!" << endl;
			nkezd = nstop;
			nstop = nstop + nconv;
			if(maxiter <= nstop)
			{
				out1 << "Error occured!" << endl;
				out1 << "Maxiter must be bigger, than " << nstop << endl;
				return;
			}
			nstep++;
			goto A;
		}
		else
		{
			out1 << endl;
			out1 << "Eigenvectors converged!!" << endl;
			out1 << endl;
		}
//		for calculating expectation values
//		write out eigenvectors
		eivectype = "Lanczos";
		ofstream eivec;
		eivec.open ("LEIVEC.txt");
		eivec.setf(ios_base::fixed,ios_base::floatfield);
		for(i = 0; i < prodnDIM; i++)
		{
			for(j = 0; j < neig; j++)
			{
				eivec << setprecision(30) << svmat[i][j] << " ";
			}
			eivec << endl;
		}
		eivec.close();
//		write out eigenvectors in binary format
//		ofstream eivecbinary("EIVE.dat", ios::out|ios::binary|ios::ate);
//		eivecbinary.write(reinterpret_cast<const char*>(svmat),prodnDIM*neig);
//		eivecbinary.close();
	}
	out1 << endl;
	cout << "Lánczos vége" << endl;
	t =  clock() - t;
	out1 << "új Lánczos: " << t << " clicks, " << (float) t / CLOCKS_PER_SEC << " seconds." << endl;
	out1.close();
	return;
}

void davidson(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *diag, int *diagind, double *ud, string &eivectype, int vector)
{
	clock_t ttt;
	clock_t t1;
	clock_t t2;
	clock_t t3;
	clock_t t4;
	clock_t t5;
	t1 = clock();
	int nrot = 100;
	int checkeigval = 1;
	const double ehcm = 219474.630;
	double norm = 0.0, lambda = 0.0, elozo = 0.0, sum = 0.0, gamma1;
//	double *q = new double[n];
	double *x = new double[n];
	double *y = new double[n];
	double *v = new double[n];
	double *r = new double[n];
	double *temp = new double[n];
//	double *temp2 = new double[n];
//	double *s = new double[n];
	double *diagtemp = new double[n];
	double *theta = new double[maxiter];
	double *beta = new double[maxiter];
	double **Q = new double *[n];
	for(int i = 0; i < n; i++)
		Q[i] = new double[maxiter];
	double **W = new double *[n];
	for(int i = 0; i < n; i++)
		W[i] = new double[maxiter];
	double **A = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		A[i] = new double[maxiter];
	double **B = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		B[i] = new double[maxiter];
	double **sv = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		sv[i] = new double[maxiter];
//	Mat_IO_DP A(maxiter,maxiter), B(maxiter,maxiter), sv(maxiter,maxiter);

	ofstream file;
	file.open("davidson.txt");
	file.setf(ios_base::fixed,ios_base::floatfield);
	
	cout << "Kezdődik a Davidson" << endl;
	
	ifstream wf;
	wf.open("thing.dat");
	wf.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < n; i++)
	{
		wf >> Q[i][0];
	}
	wf.close();
/*	CGM(n,100000,pow(10.0,-11.0),val,ncol,nrow,ud,diagind,q,diag,file);
	for(int i = 0; i < n; i++)
	{
		Q[i][0] = q[i];
	}	*/
/*	for(int i = 0; i < n; i++)
	{
		Q[i][0] = 1.0;
	}
	for(int i = 0; i < n; i++)
	{
		norm = norm + Q[i][0] * Q[i][0];
	}
	norm = sqrt(norm);
	for(int i = 0; i < n; i++)
	{
		Q[i][0] = (1.0 / norm) * Q[i][0];
	}	*/
	for(int j = 0; j < maxiter; j++)
	{
		t2 = clock();
		file << "ÚJ CIKLUS j = " << j+1 << endl;
/*		for(int i = 0; i < n; i++)
		{
			sum = 0.0;
			for(int k = 0; k < n; k++)
			{
				sum = sum + A[i][k] * Q[k][j];
			}
			W[i][j] = sum;
		}	*/
		ttt = clock();
		multfordavidson2(n,j,Q,W,val,ncol,nrow);
		ttt = clock() - ttt;
//		cout << "ttt mult = " << (float) ttt / CLOCKS_PER_SEC << " seconds." << endl;
		for(int k = 0; k <= (j - 1); k++)
		{
//			cout << "haho" << endl;
			sum = 0.0;
			for(int i = 0; i < n; i++)
			{
				sum = sum + Q[i][k] * W[i][j];
			}
			A[k][j] = A[j][k] = B[k][j] = B[j][k] = sum;
		}
		sum = 0.0;
		for(int i = 0; i < n; i++)
		{
			sum = sum + Q[i][j] * W[i][j];
		}
		A[j][j] = B[j][j] = sum;
		t2 = clock() - t2;
//		cout << "t2 (mátrix felépítése) = " << (float) t2 / CLOCKS_PER_SEC << " seconds." << endl;
		t3 = clock();
//		cout << "A = " << endl;
/*		ofstream bbb;
		bbb.open("B.dat");
		bbb.setf(ios_base::fixed,ios_base::floatfield);
		for(int i = 0; i <= j; i++)
		{
			for(int k = 0; k <= j; k++)
			{
				cout << A[i][k] << " ";
			}
			cout << endl;
		}
		bbb.close();	*/
//		powerfordavidson(B,val,ncol,nrow,s,lambda,j+1,10.0E-7,1000);
		for(int i = 0; i < j+1; i++)
		{
			for(int k = 0; k < j+1; k++)
			{
				if(i == k)
				{
					sv[i][k] = 1.00;
				}
				else
				{
					sv[i][k] = 0.00;
				}
			}
		}
		tred2fordavidson(j+1,B,theta,beta);
/*		cout << "B = " << endl;
		for(int i = 0; i <= j; i++)
		{
			for(int k = 0; k <= j; k++)
			{
				cout << B[i][k] << " ";
			}
				cout << endl;
		}
		for(int i = 0; i < j+1; i++)
		{
			cout << "theta = " << theta[i] << endl
		}	*/
		for(int i = 0; i < j+1; i++)
		{
			beta[i] = beta[i+1];
		}
/*		for(int i = 0; i < j+1; i++)
		{
			cout << "beta = " << beta[i] << endl;
		}	*/
		tqli2fordavidson(theta,beta,sv,j+1);
/*		for(int i = 0; i < j+1; i++)
		{
			cout << "theta = " << theta[i] << endl;
		}	*/
		eigsrtfordavidson(theta,sv,j+1);
//		cout << "theta = " << theta[j] << endl;
		if((abs(theta[j] - elozo) <= eps) && (checkeigval == 1))
		{
			checkeigval = 0;
			file << "Eigenvalues converged in " << j+1 << " iteration!" << endl;
			file << "lambda / hartree = " << theta[j] << " lambda / cm-1 = " << setprecision(6) << theta[j] * ehcm << endl;
			if(vector == 0)
			{
				t1 =  clock() - t1;
				file << "Davidson: " << t1 << " clicks, " << (float) t1 / CLOCKS_PER_SEC << " seconds." << endl;
				file.close();
				cout << "Davidson vége" << endl;
				return;
			}
		}
		elozo = theta[j];
		if(j > 1)
		{
/*			for(int i = 0; i < j+1; i++)
                        {
                             cout << "sv[i][j] = " << sv[i][j] << endl;
                        }	*/
			for(int i = 0; i <= j; i++)
			{
				temp[i] = 0.0;
				for(int k = 0; k <= j; k++)
				{
					temp[i] = temp[i] + B[i][k] * sv[k][j];
				}
			}
/*			for(int i = 0; i < j+1; i++)
			{
				cout << "temp!!!! = " << temp[i] << endl;
			}	*/
			for(int i = 0; i < j+1; i++)
			{
				sv[i][j] = temp[i];
			}
		}
		for(int i = 0; i <= j; i++)
		{
			for(int k = 0; k <= j; k++)
			{
				B[i][k] = A[i][k];
			}
		}
//		jacobifordavidson(j+1,B,theta,sv,nrot);
//		eigsrt(theta,sv,j+1);
		t3 = clock() - t3;
//		cout << "t3 (diagonalizálás) = " << (float) t3 / CLOCKS_PER_SEC << " seconds." << endl;
		t4 = clock();
/*		cout << "theta = " << endl;
		for(int i = 0; i < j+1; i++)
		{
			cout << theta[i] << endl;
		}	*/
/*		cout << "sv = " << endl;
		for(int i = 0; i < j+1; i++)
		{
			for(int k = 0; k < j+1; k++)
			{
				cout << sv[i][k] << " ";
			}
			cout << endl;
		}	*/
		for(int i = 0; i < n; i++)
		{
			y[i] = x[i] = 0.0;
			for(int k = 0; k < (j + 1); k++)
			{
				y[i] = y[i] + Q[i][k] * sv[k][j];
				x[i] = x[i] + W[i][k] * sv[k][j];
			}
		}
/*		if(j == restart)
		{
			file << "restarted" << endl;
			ujra = 1;
			goto A;
		}	*/
/*		for(int i = 0; i < n; i++)
		{
			sum = 0.0;
			for(int ii = 0; ii < n; ii++)
			{
				sum = sum + A[i][ii] * y[ii];
			}
			r[i] = sum;
		}	*/
//		multfordavidson1(n,y,r,val,ncol,nrow);
//		CGM(n,100000,pow(10.0,-11.0),val,ncol,nrow,ud,diagind,r,y,file);
		for(int i = 0; i < n; i++)
		{
			r[i] = x[i] - theta[j] * y[i];
		}
		norm = 0.0;
		for(int i = 0; i < n; i++)
		{
			norm = norm + r[i] * r[i];
		}
		norm = sqrt(norm);
		cout << scientific << norm << endl;
		if(norm < eps)
		{
			file << "Eigenvector is also converged in " << j+1 << " iteration!" << endl;
			file << "The Residual norm is: " << norm << endl;
			file << "Eigenvalues again: " << j+1 << " iteration!" << endl;
			file << "lambda / hartree = " << theta[j] << " lambda / cm-1 = " << setprecision(6) << theta[j] * ehcm << endl;
			eivectype = "Davidson";
			ofstream eivec;
			eivec.open("DEIVEC.txt");
			eivec.setf(ios_base::fixed,ios_base::floatfield);
			for(int i = 0; i < n; i++)
			{
				eivec << setprecision(30) << y[i] << endl;
			}
			eivec.close();
			t1 =  clock() - t1;
			file << "Davidson: " << t1 << " clicks, " << (float) t1 / CLOCKS_PER_SEC << " seconds." << endl;
			file.close();
			cout << "Davidson vége" << endl;
			return;
		}
		t4 = clock() - t4;
//		cout << "t4 (hiba számitás) = " << (float) t4 / CLOCKS_PER_SEC << " seconds." << endl;
		t5 = clock();
		transform2(n,diag,diagtemp,theta[j]);
		for(int i = 0; i < n; i++)
		{
			v[i] = diagtemp[i] * r[i];
		}
//		CGMfordavidson(n,maxiter,eps,C,v,r,theta[0],file);
/*		for(int k = 0; k < j+1; k++)
		{
			temp[k] = 0.0;
			for(int i = 0; i < n; i++)
			{
				temp[k] = temp[k] + Q[i][k] * v[i];
			}
		}
		for(int i = 0; i < n; i++)
		{
			temp2[i] = 0.0;
			for(int k = 0; k < j+1; k++)
			{
				temp2[i] = temp2[i] + Q[i][k] * temp[k];
			}
		}
		for(int i = 0; i < n; i++)
		{
			v[i] = v[i] - temp2[i];
		}	*/
		for(int i = 0; i <= j; i++)
		{
			gamma1 = 0.00;
			for(int iter = 0; iter < n; iter++)
			{
				gamma1 = gamma1 + v[iter] * Q[iter][j-i];
			}
			gamma1 = -1.0 * gamma1;
			for(int iter = 0; iter < n; iter++)
			{
				v[iter] = v[iter] + gamma1 * Q[iter][j-i];
			}
		}
		norm = 0.0;
		for(int i = 0; i < n; i++)
		{
			norm = norm + v[i] * v[i];
		}
		norm = sqrt(norm);
		for(int i = 0; i < n; i++)
		{
			Q[i][j+1] = (1.0 / norm) * v[i];
		}
		t5 = clock() - t5;
//		cout << "t5 (ortogonalizáció) = " << (float) t5 / CLOCKS_PER_SEC << " seconds." << endl;
	}
	return;
}

void jacobidavidson(int n, double m, int maxiter, double eps, double *val, int *ncol, int *nrow, double *diag, int *diagind, string &eivectype, int vector)
{
	clock_t ttt;
	clock_t t1;
	clock_t t2;
	clock_t t3;
	clock_t t4;
	clock_t t5;
	t1 = clock();
	int nrot = 100;
	int checkeigval = 1;
	const double ehcm = 219474.630;
	double norm = 0.0, lambda = 0.0, elozo = 0.0, sum = 0.0, gamma1, tau1 = 0.0, tau2 = 0.0;
//	double *q = new double[n];
	double *x = new double[n];
	double *y = new double[n];
	double *v = new double[n];
	double *r = new double[n];
//	double *w = new double[n];
	double *temp = new double[maxiter];
//	double *temp2 = new double[n];
//	double *s = new double[n];
//	double *diagtemp = new double[n];
	double *theta = new double[maxiter];
	double *beta = new double[maxiter];
	double **Q = new double *[n];
	for(int i = 0; i < n; i++)
		Q[i] = new double[maxiter];
	double **W = new double *[n];
	for(int i = 0; i < n; i++)
		W[i] = new double[maxiter];
	double **A = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		A[i] = new double[maxiter];
	double **B = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		B[i] = new double[maxiter];
	double **sv = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		sv[i] = new double[maxiter];
//	Mat_IO_DP A(maxiter,maxiter), B(maxiter,maxiter), sv(maxiter,maxiter);

	ofstream file;
	file.open("jacobidavidson.txt");
	file.setf(ios_base::fixed,ios_base::floatfield);
	
	cout << "Kezdődik a JacobiDavidson" << endl;
	
	ifstream wf;
	wf.open("thing.dat");
	wf.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < n; i++)
	{
		wf >> Q[i][0];
	}
	wf.close();
//	CGM(n,100000,pow(10.0,-12.0),val,ncol,nrow,ud,diagind,q,diag,file);
//	minres(n,10000,pow(10.0,-12.0),val,ncol,nrow,ud,q,diag);
//	minresforjacobidavidson(n,maxiter,pow(10.0,-12.0),theta,mu1,y,k,u,val,ncol,nrow,diagind,v,r_tilde);
/*	for(int i = 0; i < n; i++)
	{
		Q[i][0] = 1.0;
	}
	for(int i = 0; i < n; i++)
	{
		norm = norm + Q[i][0] * Q[i][0];
	}
	norm = sqrt(norm);
	for(int i = 0; i < n; i++)
	{
		Q[i][0] = (1.0 / norm) * Q[i][0];
	}	*/
	for(int j = 0; j < maxiter; j++)
	{
		t2 = clock();
		file << "ÚJ CIKLUS j = " << j+1 << endl;
/*		for(int i = 0; i < n; i++)
		{
			q[i] = Q[i][j];
		}	*/
		ttt = clock();
//		CGM(n,100000,pow(10.0,-11.0),val,ncol,nrow,ud,diagind,w,q,file);
		multfordavidson2(n,j,Q,W,val,ncol,nrow);
		ttt = clock() - ttt;
//		cout << "ttt mult = " << (float) ttt / CLOCKS_PER_SEC << " seconds." << endl;
/*		for(int i = 0; i < n; i++)
		{
			W[i][j] = w[i];
		}	*/
		for(int k = 0; k <= (j - 1); k++)
		{
			sum = 0.0;
			for(int i = 0; i < n; i++)
			{
				sum = sum + Q[i][k] * W[i][j];
			}
			A[k][j] = A[j][k] = B[k][j] = B[j][k] = sum;
		}
		sum = 0.0;
		for(int i = 0; i < n; i++)
		{
			sum = sum + Q[i][j] * W[i][j];
		}
//		cout << sum << endl;
		A[j][j] = B[j][j] = sum;
		t2 = clock() - t2;
//		cout << "t2 (mátrix felépítése) = " << (float) t2 / CLOCKS_PER_SEC << " seconds." << endl;
		t3 = clock();
		for(int i = 0; i < j+1; i++)
		{
			for(int k = 0; k < j+1; k++)
			{
				if(i == k)
				{
					sv[i][k] = 1.00;
				}
				else
				{
					sv[i][k] = 0.00;
				}
			}
		}
		tred2fordavidson(j+1,B,theta,beta);
		for(int i = 0; i < j+1; i++)
		{
			beta[i] = beta[i+1];
		}
		tqli2fordavidson(theta,beta,sv,j+1);
		eigsrtfordavidson(theta,sv,j+1);
//		cout << "theta = " << theta[j] << endl;
		if((abs(theta[j] - elozo) <= eps) && (checkeigval == 1))
		{
			checkeigval = 0;
			file << "Eigenvalues converged in " << j+1 << " iteration!" << endl;
			file << "lambda / hartree = " << theta[j] << " lambda / cm-1 = " << setprecision(6) << theta[j] * ehcm << endl;
			if(vector == 0)
			{
				t1 =  clock() - t1;
				file << "JacobiDavidson: " << t1 << " clicks, " << (float) t1 / CLOCKS_PER_SEC << " seconds." << endl;
				file.close();
				cout << "JacobiDavidson vége" << endl;
				return;
			}
		}
		elozo = theta[j];
		if(j > 0)
		{
			for(int i = 0; i <= j; i++)
			{
				temp[i] = 0.0;
				for(int k = 0; k <= j; k++)
				{
					temp[i] = temp[i] + B[i][k] * sv[k][j];
				}
			}
			for(int i = 0; i < j+1; i++)
			{
				sv[i][j] = temp[i];
			}
		}
		for(int i = 0; i <= j; i++)
		{
			for(int k = 0; k <= j; k++)
			{
				B[i][k] = A[i][k];
			}
		}
		t3 = clock() - t3;
//		cout << "t3 (diagonalizálás) = " << (float) t3 / CLOCKS_PER_SEC << " seconds." << endl;
		t4 = clock();
		for(int i = 0; i < n; i++)
		{
			y[i] = x[i] = 0.0;
			for(int k = 0; k < (j + 1); k++)
			{
				y[i] = y[i] + Q[i][k] * sv[k][j];
				x[i] = x[i] + W[i][k] * sv[k][j];
			}
		}
//		multfordavidson1(n,y,r,val,ncol,nrow);
//		CGM(n,100000,pow(10.0,-11.0),val,ncol,nrow,ud,diagind,r,y,file);
		for(int i = 0; i < n; i++)
		{
			r[i] = x[i] - theta[j] * y[i];
		}
		norm = 0.0;
		for(int i = 0; i < n; i++)
		{
			norm = norm + r[i] * r[i];
		}
		norm = sqrt(norm);
		cout << scientific << norm << endl;
		if(norm < eps)
		{
			file << "Eigenvector is also converged in " << j+1 << " iteration!" << endl;
			file << "The Residual norm is: " << norm << endl;
			file << "Eigenvalues again: " << j+1 << " iteration!" << endl;
			file << "lambda / hartree = " << theta[j] << " lambda / cm-1 = " << setprecision(6) << theta[j] * ehcm << endl;
			eivectype = "JacobiDavidson";
			ofstream eivec;
			eivec.open("JDEIVEC.txt");
			eivec.setf(ios_base::fixed,ios_base::floatfield);
			for(int i = 0; i < n; i++)
			{
				eivec << setprecision(30) << y[i] << endl;
			}
			eivec.close();
			t1 =  clock() - t1;
			file << "JacobiDavidson: " << t1 << " clicks, " << (float) t1 / CLOCKS_PER_SEC << " seconds." << endl;
			file.close();
			cout << "JacobiDavidson vége" << endl;
			return;
		}
		t4 = clock() - t4;
//		cout << "t4 (hiba számitás) = " << (float) t4 / CLOCKS_PER_SEC << " seconds." << endl;
		t5 = clock();
//		jdcorrectionequation(n,m,maxiter,eps,theta[j],v,r,y,val,ud,ncol,nrow,diagind,file);
		jdcorrectionequation_tuti(n,m,maxiter,eps,theta[j],v,r,y,val,diag,ncol,nrow,diagind,j+1,file);
/*		tau1 = tau2 = 0.0;
		for(int i = 0; i < n; i++)
		{
			tau1 = tau1 + v[i] * v[i];
		}
		for(int k = 0; k < j+1; k++)
		{
			temp[k] = 0.0;
			for(int i = 0; i < n; i++)
			{
				temp[k] = temp[k] + Q[i][k] * v[i];
			}
		}
		for(int i = 0; i < n; i++)
		{
			temp2[i] = 0.0;
			for(int k = 0; k < j+1; k++)
			{
				temp2[i] = temp2[i] + Q[i][k] * temp[k];
			}
		}
		for(int i = 0; i < n; i++)
		{
			v[i] = v[i] - temp2[i];
		}
		for(int i = 0; i < n; i++)
		{
			tau2 = tau2 + v[i] * v[i];
		}
		if((tau2 / tau1) < 0.1)
		{
			cout << "MOST" << endl;
			for(int k = 0; k < j+1; k++)
			{
				temp[k] = 0.0;
				for(int i = 0; i < n; i++)
				{
					temp[k] = temp[k] + Q[i][k] * v[i];
				}
			}
			for(int i = 0; i < n; i++)
			{
				temp2[i] = 0.0;
				for(int k = 0; k < j+1; k++)
				{
					temp2[i] = temp2[i] + Q[i][k] * temp[k];
				}
			}
			for(int i = 0; i < n; i++)
			{
				v[i] = v[i] - temp2[i];
			}
		}	*/
		for(int i = 0; i <= j; i++)
		{
			gamma1 = 0.00;
			for(int iter = 0; iter < n; iter++)
			{
				gamma1 = gamma1 + v[iter] * Q[iter][j-i];
			}
			gamma1 = -1.0 * gamma1;
			for(int iter = 0; iter < n; iter++)
			{
				v[iter] = v[iter] + gamma1 * Q[iter][j-i];
			}
		}
		norm = 0.0;
		for(int i = 0; i < n; i++)
		{
			norm = norm + v[i] * v[i];
		}
		norm = sqrt(norm);
		for(int i = 0; i < n; i++)
		{
			Q[i][j+1] = (1.0 / norm) * v[i];
		}
		t5 = clock() - t5;
//		cout << "t5 (ortogonalizáció) = " << (float) t5 / CLOCKS_PER_SEC << " seconds." << endl;
	}
	return;
}

void jdcorrectionequation_tuti(int n, int nonzero, int maxiter, double eps, double theta, double* v, double* r, double *u, double* val, double* diag, int* ncol, int* nrow, int* diagind, int step, ofstream &file)
{
	double *y = new double[n];
	double *k = new double[n];
	double *r_hat = new double[n];
	double *r_tilde = new double[n];
	double mu1, mu2;
	
	for(int i = 0; i < n; i++)
	{
		k[i] = diag[i] - theta;
		k[i] = 1.0 / k[i];
	}
	for(int i = 0; i < n; i++)
	{
		y[i] = k[i] * u[i];
	}
//	ILUprecondforjacobidavidson(n,theta,val,ncol,nrow,k,diagind,y,u);
	mu1 = 0.0;
	for(int i = 0; i < n; i++)
	{
		mu1 = mu1 + u[i] * y[i];
	}
	for(int i = 0; i < n; i++)
	{
		r_hat[i] = k[i] * r[i];
	}
	mu2 = 0.0;
	for(int i = 0; i < n; i++)
	{
		mu2 = mu2 + u[i] * r_hat[i];
	}
	for(int i = 0; i < n; i++)
	{
		r_tilde[i] = r_hat[i] - (mu2 / mu1) * y[i];
		r_tilde[i] = -1.0 * r_tilde[i];
	}
//	CGMforjacobidavidson_tuti(n,10,pow(10.0,-12.0),theta,mu1,y,k,u,val,ncol,nrow,diagind,v,r_tilde,file);
	minresforjacobidavidson(n,15,pow(10.0,-12.0),theta,mu1,y,k,u,val,ncol,nrow,diagind,v,r_tilde);
	return;
}

void jdcorrectionequation(int n, int nonzero, int maxiter, double eps, double lambda, double* v, double* r, double *y, double* val, double* ud, int* ncol, int* nrow, int* diagind, ofstream &file)
{
	double **Y = new double *[n];
	for(int i = 0; i < n; i++)
		Y[i] = new double[n];
	double **A = new double *[n];
	for(int i = 0; i < n; i++)
		A[i] = new double[n];
	double **AY = new double *[n];
	for(int i = 0; i < n; i++)
		AY[i] = new double[n];
	double sum = 0.0;
	
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			sum = 0.0;
			sum = sum + y[i] * y[j];
			if(i == j)
			{
				Y[i][j] = 1.00 - sum;
			}
			else
			{
				Y[i][j] = -sum;
			}
		}
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			A[i][j] = 0.0;
		}
	}
	comp2mtx(n,nonzero,A,val,ud,ud,ncol,nrow);
	for(int i = 0; i < n; i++)
	{
		A[i][i] = A[i][i] - lambda;
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			sum = 0.0;
			for(int k = 0; k < n; k++)
			{
				sum = sum + A[i][k] * Y[k][j];
			}
			AY[i][j] = sum;
		}
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			sum = 0.0;
			for(int k = 0; k < n; k++)
			{
				sum = sum + Y[i][k] * AY[k][j];
			}
			A[i][j] = sum;
		}
	}
	CGMfordavidson(n,maxiter,eps,A,v,r,lambda,file);
	return;
}

void precond_nostore(int prodnDIM, double *z, double *v, double *r, double *val, int *ncol, int *nrow, int *diagind, int *IOPT)
{
//       z(1, ...) --> scratch,  v-->scratch, later output, r--> input, scratch later

	double delta, c, t, d, tmp, tmp2;
	double wprecond;
	int m, k, i, po1, po2, prodnDIM2;
	double *temp = new double[prodnDIM];
//	Vec_DP temp(prodnDIM);
	 
/*	for(i = 0; i < prodnDIM; i++)
	{
		cout << "v[i] = q2[i] = " << v[i] << endl;
		cout << "r[i] = q1[i] = " << r[i] << endl;
	} */
	 
	delta = pow(10.0,(-(double)IOPT[64]));
	m = IOPT[65];
	prodnDIM2 = 2 * prodnDIM;
	        
	for(i = 0; i < prodnDIM; i++)
	{
		v[i] = 0.00;
	}
	for(i = 0; i < prodnDIM2; i++)
	{
		z[i] = 0.00;
	}
	for(i = 0; i < prodnDIM; i++)
	{
		z[i] = r[i] * z[2*prodnDIM + i]; // z = r / diag read in svl
	}
/*	for(i = 0; i < prodnDIM; i++)
	{
		temp[i] = z[2*prodnDIM + i];
	}
	ILUprecond(prodnDIM,val,ncol,nrow,temp,diagind,z,r);	*/
/*	for(i = 0; i < prodnDIM; i++)
	{
		cout << "z[i] = " << z[i] << endl;
	}	*/
	for(i = 0; i < prodnDIM; i++)
	{
		v[i] = z[i];
	}
/*	for(i = 0; i < prodnDIM; i++)
	{
		cout << "v[i] = " << v[i] << endl;
	}	*/
	c = 0.00;
	for(i = 0; i < prodnDIM; i++)
	{
		  c = c + z[i] * r[i];
	}
	for(k = 1; k <= m; k++)
	{
		tmp = 0.00;
		tmp2 = 0.00;
		for(i = 0; i < prodnDIM; i++)
		{
			tmp2 = tmp2 + v[i] * v[i];
		}
		tmp = sqrt(tmp2);
		if(tmp < delta)
		{
			goto A;
		}
		mult(prodnDIM,v,z,val,ncol,nrow);
		tmp = 0.00;
		for(i = 0; i < prodnDIM; i++)
		{
			tmp = tmp + v[i] * z[i];
		}
		t = c / tmp;
		for(i = 0; i < prodnDIM; i++)
		{
			z[prodnDIM+i] = z[prodnDIM+i] + v[i] * t;
		}
		t = -1.0 * t;   
		for(i = 0; i < prodnDIM; i++)
		{
			r[i] = r[i] + z[i] * t;
		}   
		for(i = 0; i < prodnDIM; i++)
		{
			z[i] = 0.00;
		}
		for(i = 0; i < prodnDIM; i++)
		{
			z[i] = r[i] * z[2*prodnDIM + i]; // z = r / diag read in svl
		}
/*		for(i = 0; i < prodnDIM; i++)
		{
			temp[i] = z[2*prodnDIM + i];
		}
		ILUprecond(prodnDIM,val,ncol,nrow,temp,diagind,z,r);	*/
		d = 0.00;
		for(i = 0; i < prodnDIM; i++)
		{
			d = d + r[i] * z[i];
		}
		tmp = d / c;
		for(i = 0; i < prodnDIM; i++)
		{
			v[i] =  v[i] * tmp;
		}
		tmp=1.00;     
		for(i = 0; i < prodnDIM; i++)
		{
			v[i] = v[i] + z[i];
		}
		c = d;
	}
A:	;
	for(i = 0; i < prodnDIM; i++)
	{
		v[i] = z[prodnDIM+i];
	}
//	cout << k << endl;
	return;
}

void CGM2LANCZOS(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *ud, int *diagind, double *x2, double *b, ofstream &out1)
{
	clock_t t_cgm;
	t_cgm = clock();
	int conv;
	double *x1 = new double[n];
	double *z = new double[n];
	double *r = new double[n];
	double *p = new double[n];
	double *Ap = new double[n];
	double *q = new double[n];
	double **B = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		B[i] = new double[maxiter];
	double **P = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		P[i] = new double[maxiter];
	double **R = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		R[i] = new double[maxiter];
	double **LAMBDA = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		LAMBDA[i] = new double[maxiter];
	double **DELTA = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		DELTA[i] = new double[maxiter];
//	Vec_DP x1(n), z(n), r(n), p(n), q(n);
	double alpha = 0.0, beta = 0.0, rho1 = 0.0, rho2 = 0.0, norm;
	
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if(i == j)
				B[i][j] = 1.0;
			else
				B[i][j] = 0.0;
		}
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			LAMBDA[i][j] = 0.0;
		}
	}
//	initial guess for x
	for(int i = 0; i < n; i++)
	{
		x1[i] = 1.0;
	}
	mult(n,x1,r,val,ncol,nrow);
//	calc r = b - Ax
	for(int i = 0; i < n; i++)
	{
		r[i] = b[i] - r[i];
	}
	for(int i = 1; i <= maxiter; i++)
	{
/*		for(int j = 0; j < n; j++)
		{
			z[j] = r[j];// * ud[j];	// nincs prekond; z[j] = r[j] * ud[j]: Jacobi-prekond
		}	*/
		ILUprecond(n,val,ncol,nrow,ud,diagind,z,r);	// ILU prekond
		rho2 = 0.0;
		for(int j = 0; j < n; j++)
		{
			rho2 = rho2 + r[j] * z[j];
		}
		if(i == 1)
		{
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j];
			}
			for(int j = 0; j < n; j++)
			{
				P[j][i-1] = p[j];
			}
			mult(n,p,Ap,val,ncol,nrow);
			for(int j = 0; j < n; j++)
			{
				LAMBDA[i-1][i-1] = LAMBDA[i-1][i-1] + p[j] * Ap[j];
			}
		}
		else
		{
			beta = rho2 / rho1;
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j] + beta * p[j];
			}
			B[i-1][i] = -beta;
			for(int j = 0; j < n; j++)
			{
				P[j][i-1] = p[j];
			}
		}
		
		mult(n,p,q,val,ncol,nrow);
		alpha = 0.0;
		for(int j = 0; j < n; j++)
		{
			alpha = alpha + p[j] * q[j];
		}
		alpha = rho2 / alpha;
		for(int j = 0; j < n; j++)
		{
			x2[j] = x1[j] + alpha * p[j];
		}
/*		conv = 0;
		for(int j = 0; j < n; j++)
		{
			if(abs(x2[j] - x1[j]) > eps)
			{
				conv++;
			}
		}
		if(conv == 0)
		{
			cout << i << endl;
			t_cgm = clock() - t_cgm;
//			cout << "cgm = " << (float) t_cgm / CLOCKS_PER_SEC << " seconds." << endl;
			return;
		}	*/
		for(int j = 0; j < n; j++)
		{
			r[j] = r[j] - alpha * q[j];
		}
		for(int j = 0; j < n; j++)
		{
			R[i-1][j] = r[j];
		}
		norm = 0.0;
		for(int j = 0; j < n; j++)
		{
			norm = norm + r[j] * r[j];
		}
		norm = sqrt(norm);
		
		if(norm <= eps)
		{
			cout << i << endl;
			t_cgm = clock() - t_cgm;
//			cout << "cgm = " << (float) t_cgm / CLOCKS_PER_SEC << " seconds." << endl;
			return;
		}
		if(i == maxiter)
		{
			cout << i << endl;
			t_cgm = clock() - t_cgm;
//			cout << "cgm = " << (float) t_cgm / CLOCKS_PER_SEC << " seconds." << endl;
			return;
		}
		for(int j = 0; j < n; j++)
		{
			x1[j] = x2[j];
		}
		rho1 = rho2;
	}
	out1 << "Nem konvergált a preCGM " << maxiter << " lépésben ! :(" << endl;
	return;
}

void CGM(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *ud, int *diagind, double *x2, double *b, ofstream &out1, int precondtype, ofstream &out2)
{
	clock_t t_cgm;
	t_cgm = clock();
	int conv;
	double *x1 = new double[n];
	double *z = new double[n];
	double *r = new double[n];
	double *p = new double[n];
	double *q = new double[n];
//	Vec_DP x1(n), z(n), r(n), p(n), q(n);
	double alpha = 0.0, beta = 0.0, rho1 = 0.0, rho2 = 0.0, norm;
	
//	initial guess for x
	for(int i = 0; i < n; i++)
	{
		x1[i] = 0.0;
	}
/*	norm = 0.0;
	for(int i = 0; i < n; i++)
	{
		norm = norm + x1[i] * x1[i];
	}
	norm = sqrt(norm);
	for(int i = 0; i < n; i++)
	{
		x1[i] = (1.0 / norm) * x1[i];
	}	*/
	mult(n,x1,r,val,ncol,nrow);
//	calc r = b - Ax
	for(int i = 0; i < n; i++)
	{
		r[i] = b[i] - r[i];
	}
	for(int i = 1; i <= maxiter; i++)
	{
		if(precondtype == 1)
		{
			for(int j = 0; j < n; j++)
			{
				z[j] = r[j];	// nincs prekond
			}
		}
		if(precondtype == 2)
		{
			for(int j = 0; j < n; j++)
			{
				z[j] = r[j] * ud[j];	// z[j] = r[j] * ud[j]: Jacobi-prekond
			}
		}
		if(precondtype == 3)
		{
			ILUprecond(n,val,ncol,nrow,ud,diagind,z,r);	// ILU prekond
		}
		rho2 = 0.0;
		for(int j = 0; j < n; j++)
		{
			rho2 = rho2 + r[j] * z[j];
		}
		if(i == 1)
		{
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j];
			}
		}
		else
		{
			beta = rho2 / rho1;
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j] + beta * p[j];
			}
		}
		mult(n,p,q,val,ncol,nrow);
		alpha = 0.0;
		for(int j = 0; j < n; j++)
		{
			alpha = alpha + p[j] * q[j];
		}
		alpha = rho2 / alpha;
		for(int j = 0; j < n; j++)
		{
			x2[j] = x1[j] + alpha * p[j];
		}
/*		conv = 0;
		for(int j = 0; j < n; j++)
		{
			if(abs(x2[j] - x1[j]) > eps)
			{
				conv++;
			}
		}
		if(conv == 0)
		{
//			cout << i << endl;
			t_cgm = clock() - t_cgm;
//			cout << "cgm = " << (float) t_cgm / CLOCKS_PER_SEC << " seconds." << endl;
			return;
		}	*/
		for(int j = 0; j < n; j++)
		{
			r[j] = r[j] - alpha * q[j];
		}
		norm = 0.0;
		for(int j = 0; j < n; j++)
		{
			norm = norm + r[j] * r[j];
		}
		norm = sqrt(norm);
		if(norm <= eps)
		{
//			cout << i << endl;
			t_cgm = clock() - t_cgm;
//			cout << "cgm = " << (float) t_cgm / CLOCKS_PER_SEC << " seconds." << endl;
			return;
		}
		if(i == maxiter)
		{
			out2 << i << endl;
			t_cgm = clock() - t_cgm;
//			cout << "cgm = " << (float) t_cgm / CLOCKS_PER_SEC << " seconds." << endl;
			return;
		}
		for(int j = 0; j < n; j++)
		{
			x1[j] = x2[j];
		}
		rho1 = rho2;
	}
	out1 << "Nem konvergált a preCGM " << maxiter << " lépésben ! :(" << endl;
	return;
}

void ILUprecond(int n, double *val, int *ncol, int *nrow, double *ud, int *diagind, double *y, double *x)
{
	double sum;
	double *z = new double[n];
//	Vec_DP z(n);
	
	for(int i = 0; i < n; i++)
	{
		z[i] = y[i] = 0.0;
	}
	for(int i = 0; i < n; i++)
	{
		sum =  0.0;
		for(int j = nrow[i]; j <= diagind[i]; j++)
		{
			sum = sum + val[j] * z[ncol[j]];
		}
		z[i] = ud[i] * (x[i] - sum);
	}
	for(int i = n - 1; i >= 0; i--)
	{
		sum = 0.0;
		for(int j = diagind[i]; j < nrow[i+1]; j++)
		{
			sum = sum + val[j] * y[ncol[j]];
			y[i] = z[i] - ud[i] * sum;
		}
	}
	return;
}

void ILUprecondforjacobidavidson(int n, double theta, double *val, int *ncol, int *nrow, double *ud, int *diagind, double *y, double *x)
{
	double sum;
	double *z = new double[n];
//	Vec_DP z(n);
	
	for(int i = 0; i < n; i++)
	{
		z[i] = y[i] = 0.0;
	}
	for(int i = 0; i < n; i++)
	{
		sum =  0.0;
		for(int j = nrow[i]; j <= diagind[i]; j++)
		{
			if(j == diagind[i])
			{
				sum = sum + (val[j] - theta) * z[ncol[j]];
			}
			else
			{
				sum = sum + val[j] * z[ncol[j]];
			}
		}
		z[i] = ud[i] * (x[i] - sum);
	}
	for(int i = n - 1; i >= 0; i--)
	{
		sum = 0.0;
		for(int j = diagind[i]; j < nrow[i+1]; j++)
		{
			if(j == diagind[i])
			{
				sum = sum + (val[j] - theta) * y[ncol[j]];
			}
			else
			{
				sum = sum + val[j] * y[ncol[j]];
			}
			y[i] = z[i] - ud[i] * sum;
		}
	}
	return;
}

void ILUprecondfordavidson(int n, double **A, double *x, double *y)
{
	double **L = new double *[n];
	for(int i = 0; i < n; i++)
		L[i] = new double[n];
	double **U = new double *[n];
	for(int i = 0; i < n; i++)
		U[i] = new double[n];
	double *z = new double[n];
	double *temp = new double[n];
//	Mat_DP L(n,n), U(n,n);
//	Vec_DP z(n), temp(n);
	
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if(j <= i - 1)
				L[i][j] = A[i][j];
			else
				L[i][j] = 0.0;
			if(j >= i + 1)
				U[i][j] = A[i][j];
			else
				U[i][j] = 0.0;
		}
	}
	for(int i = 0; i < n; i++)
	{
		temp[i] = 0.0;
		if(i == 0)
		{
			z[i] = (1.0 / A[i][i]) * y[i];
		}
		else
		{
			for(int j = 0; j < i; j++)
			{
				temp[i] = temp[i] + L[i][j] * z[j];
			}
			z[i] = (1.0 / A[i][i]) * (y[i] - temp[i]);
		}
	}
	for(int i = n - 1; i >= 0; i--)
	{
		temp[i] = 0.0;
		if(i == n - 1)
		{
			x[i] = z[i];
		}
		else
		{
			for(int j = n - 1; j > i; j--)
			{
				temp[i] = temp[i] + U[i][j] * x[j];
			}
			x[i] = z[i] - (1.0 / A[i][i]) * temp[i];
		}
	}
	return;
}

void CGMfordavidson(int n, int maxiter, double eps, double **A, double *x2, double *b, double theta, ofstream &out1)
{
	int conv;
	double *x1 = new double[n];
	double *z = new double[n];
	double *r = new double[n];
	double *p = new double[n];
	double *q = new double[n];
	double **B = new double *[n];
	for(int i = 0; i < n; i++)
		B[i] = new double[n];
//	Vec_DP x1(n), z(n), r(n), p(n), q(n);
	double alpha = 0.0, beta = 0.0, rho1 = 0.0, rho2 = 0.0;
	
//	initial guess for x
	for(int i = 0; i < n; i++)
	{
		x1[i] = 1.0;
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			B[i][j] = A[i][j];
		}
	}
	for(int i = 0; i < n; i++)
	{
//		A[i][i] = A[i][i] - theta;
	}
	for(int i = 0; i < n; i++)
	{
		r[i] = 0.0;
		for(int j = 0; j < n; j++)
		{
			r[i] = r[i] + A[i][j] * x1[j];
		}
	}
//	calc r = b - Ax
	for(int i = 0; i < n; i++)
	{
		r[i] = b[i] - r[i];
	}
	for(int i = 1; i <= maxiter; i++)
	{
		for(int j = 0; j < n; j++)
		{
			z[j] = r[j] / A[j][j];
		}
//		ILUprecondfordavidson(n,A,z,r);
		rho2 = 0.0;
		for(int j = 0; j < n; j++)
		{
			rho2 = rho2 + r[j] * z[j];
		}
		if(i == 1)
		{
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j];
			}
		}
		else
		{
			beta = rho2 / rho1;
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j] + beta * p[j];
			}
		}
		for(int j = 0; j < n; j++)
		{
			q[j] = 0.0;
			for(int k = 0; k < n; k++)
			{
				q[j] = q[j] + A[j][k] * p[k];
			}
		}
		alpha = 0.0;
		for(int j = 0; j < n; j++)
		{
			alpha = alpha + p[j] * q[j];
		}
		alpha = rho2 / alpha;
		for(int j = 0; j < n; j++)
		{
			x2[j] = x1[j] + alpha * p[j];
		}
/*		for(int j = 0; j < n; j++)
		{
			cout << "x2[" << j << "] = " << x2[j] << "     x1[" << j << "] = " << x1[j] << endl;
		}
		for(int j = 0; j < n; j++)
		{
			cout << "előtt r[" << j << "] = " << r[j] << endl;
		}
		cout << "előtt alpha = " << alpha << endl;
		for(int j = 0; j < n; j++)
		{
			cout << "előtt q[" << j << "] = " << q[j] << endl;
		}	*/
		for(int j = 0; j < n; j++)
		{
			r[j] = r[j] - alpha * q[j];
		}
/*		for(int j = 0; j < n; j++)
		{
			cout << "után r[" << j << "] = " << r[j] << endl;
		}
		cout << "után alpha = " << alpha << endl;
		for(int j = 0; j < n; j++)
		{
			cout << "után q[" << j << "] = " << q[j] << endl;
		}	*/
		conv = 0;
		for(int j = 0; j < n; j++)
		{
			if(abs(x2[j] - x1[j]) > eps)
			{
				conv++;
			}
		}
		if(conv == 0)
		{
			out1 << "i = " << i << endl;
			for(int k = 0; k < n; k++)
			{
				for(int j = 0; j < n; j++)
				{
					A[k][j] = B[k][j];
				}
			}
			return;
		}
		for(int j = 0; j < n; j++)
		{
			x1[j] = x2[j];
		}
		rho1 = rho2;
	}
	cout << "Nem konvergált CGMfordavidson! :(" << endl;
	return;
}

void chebyshev(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *ud, double *x2, double *b, double lambda_min, double lambda_max)
{
	int conv;
	double *x1 = new double[n];
	double *r1 = new double[n];
	double *r2 = new double[n];
	double *z = new double[n];
	double *p = new double[n];
	double *q = new double[n];
//	Vec_DP x1(n), r(n), z(n), p(n);
	double alpha, beta, c, d, norm;
//	initial guess for x
	for(int i = 0; i < n; i++)
	{
		x1[i] = 1.0;
	}
	mult(n,x1,r1,val,ncol,nrow);
	for(int i = 0; i < n; i++)
	{
		r1[i] = b[i] - r1[i];
	}
	d = (lambda_max + lambda_min) / 2.0;
	c = (lambda_max - lambda_min) / 2.0;
	for(int i = 1; i <= maxiter; i++)
	{
		for(int j = 0; j < n; j++)
		{
			z[j] = r1[j];// * ud[j];
		}
		if(i == 1)
		{
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j];
			}
			alpha = 2.0 / d;
		}
		else
		{
			beta = 0.0;
			beta = ((c * alpha) / 2.0) * ((c * alpha) / 2.0);
			alpha = 0.0;
			alpha = 1.0 / (d - beta);
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j] + beta * p[j];
			}
		}
		for(int j = 0; j < n; j++)
		{
			x2[j] = x1[j] + alpha * p[j];
		}
		conv = 0;
		for(int j = 0; j < n; j++)
		{
			if(abs(x2[j] - x1[j]) > eps)
			{
				conv++;
			}
		}
		if(conv == 0)
		{
			cout << i << endl;
			return;
		}
		mult(n,p,q,val,ncol,nrow);
		for(int j = 0; j < n; j++)
		{
			r2[j] = r1[j] - alpha * q[j];
		}
/*		norm = 0.0;
		for(int j = 0; j < n; j++)
		{
			norm = norm + r2[j] * r2[j];
		}
		norm = sqrt(norm);
		if(norm <= eps)
		{
			cout << i << endl;
			return;
		}	*/
		if(i == maxiter)
		{
			cout << i << endl;
			return;
		}
		for(int j = 0; j < n; j++)
		{
			x1[j] = x2[j];
			r1[j] = r2[j];
		}
	}
	return;
}

void BiCGM(int n, int maxiter, double eps, double *val2, int *ncol, int *nrow, double *ud, double *x2, double *b)
{
	int conv;
	double *x1 = new double[n];
	double *r = new double[n];
	double *rtilde = new double[n];
	double *z = new double[n];
	double *ztilde = new double[n];
	double *p = new double[n];
	double *ptilde = new double[n];
	double *q = new double[n];
	double *qtilde = new double[n];
//	Vec_DP x1(n), r(n), rtilde(n), z(n), ztilde(n);
//	Vec_DP p(n), ptilde(n), q(n), qtilde(n);
	double alpha, beta, rho1, rho2, norm, normtilde;
	
	for(int i = 0; i < n; i++)
	{
		x1[i] = 1.0;
	}
	mult(n,x1,r,val2,ncol,nrow);
	for(int i = 0; i < n; i++)
	{
		r[i] = b[i] - r[i];
	}
	for(int i = 0; i < n; i++)
	{
		rtilde[i] = r[i];
	}
	for(int i = 1; i <= maxiter; i++)
	{
		for(int j = 0; j < n; j++)
		{
			z[j] = r[j];// * ud[j];
			ztilde[j] = rtilde[j];// * ud[j];
		}
		rho2 = 0.0;
		for(int j = 0; j < n; j++)
		{
			rho2 = rho2 + z[j] * rtilde[j];
		}
		if(i == 1)
		{
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j];
				ptilde[j] = ztilde[j];
			}
		}
		else
		{
			beta = rho2 / rho1;
			for(int j = 0; j < n; j++)
			{
				p[j] = z[j] + beta * p[j];
				ptilde[j] = ztilde[j] + beta * ptilde[j];
			}
		}
		mult(n,p,q,val2,ncol,nrow);
		mult2(n,ptilde,qtilde,val2,ncol,nrow);
		alpha = 0.0;
		for(int j = 0; j < n; j++)
		{
			alpha = alpha + ptilde[j] * q[j];
		}
		alpha = rho2 / alpha;
		for(int j = 0; j < n; j++)
		{
			x2[j] = x1[j] + alpha * p[j];
		}
		conv = 0;
		for(int j = 0; j < n; j++)
		{
			if(abs(x2[j] - x1[j]) > eps)
			{
				conv++;
			}
		}
		if(conv == 0)
		{
			cout << "i = " << i << endl;
			return;
		}
		for(int j = 0; j < n; j++)
		{
			r[j] = r[j] - alpha * q[j];
			rtilde[j] = rtilde[j] - alpha * qtilde[j];
		}
/*		norm = normtilde = 0.0;
		for(int j = 0; j < n; j++)
		{
			norm = norm + r[j] * r[j];
			normtilde = normtilde + rtilde[j] * rtilde[j];
		}
		norm = sqrt(norm);
		normtilde = sqrt(normtilde);
		if((norm <= eps) && (normtilde <= eps))
		{
			cout << i << endl;
			return;
		}	*/
		if(i == maxiter)
		{
			cout << "i = " << i << endl;
			return;
		}
		for(int j = 0; j < n; j++)
		{
			x1[j] = x2[j];
		}
		rho1 = rho2;
	}
	return;
}

void minres(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double *ud, double *x2, double *b, ofstream &out2)
{
	double *x1 = new double[n];
	double *v1 = new double[n];
	double *v2 = new double[n];
	double *v3 = new double[n];
	double *Hv2 = new double[n];
	double *w1 = new double[n];
	double *w2 = new double[n];
	double *w3 = new double[n];
//	Vec_DP x1(n), v1(n), v2(n), v3(n), Hv2(n);
//	Vec_DP w1(n), w2(n), w3(n);
	double alpha, beta1 = 0.0, beta2 = 0.0;
	double gamma1 = 1.0, gamma2 = 1.0, gamma3;
	double delta, eta, sigma1 = 0.0, sigma2 = 0.0, sigma3;
	double rho1, rho2, rho3, r1norm, r2norm;
	int conv;
	
//	initial guess for x; v_i-1 = 0; w_i-2 = w_i-1 = w_i = 0
	for(int i = 0; i < n; i++)
	{
		x1[i] = 1.0;
		v1[i] = 0.0;
		w1[i] = w2[i] = 0.0;
	}
/*	for(int i = 0; i < n; i++)
	{
		r1norm = r1norm + x1[i] * x1[i];
	}
	r1norm = sqrt(r1norm);
	for(int i = 0; i < n; i++)
	{
		x1[i] = (1.0 / r1norm) * x1[i];
	}	*/
	mult(n,x1,v2,val,ncol,nrow);
/*	for(int i = 0; i < n; i++)
	{
		v2[i] = v2[i] * ud[i];
	}	*/
	for(int i = 0; i < n; i++)
	{
		v2[i] = b[i] - v2[i];
	}
	for(int i = 0; i < n; i++)
	{
		beta1 = beta1 + v2[i] * v2[i];
	}
	beta1 = sqrt(beta1);
	r1norm = beta1;
	eta = beta1;
	for(int i = 1; i <= maxiter; i++)
	{
//		cout << "i = " << i << endl;
		for(int j = 0; j < n; j++)
		{
			v2[j] = v2[j] / beta1;
		}
		mult(n,v2,Hv2,val,ncol,nrow);	// Hv_i = H * v_i
/*		for(int j = 0; j < n; j++)
		{
			Hv2[j] = Hv2[j] * ud[j];
		}	*/
		alpha = 0.0;
		for(int j = 0; j < n; j++)
		{
			alpha = alpha + v2[j] * Hv2[j];
		}
		for(int j = 0; j < n; j++)
		{
			v3[j] = Hv2[j] - alpha * v2[j] - beta1 * v1[j];
		}
		beta2 = 0.0;
		for(int j = 0; j < n; j++)
		{
			beta2 = beta2 + v3[j] * v3[j];
		}
		beta2 = sqrt(beta2);
		
		delta = gamma2 * alpha - gamma1 * sigma2 * beta1;
		rho1 = sqrt(delta * delta + beta2 * beta2);
		rho2 = sigma2 * alpha + gamma1 * gamma2 * beta1;
		rho3 = sigma1 * beta1;
		
		gamma3 = delta / rho1;
		sigma3 = beta2 / rho1;
		
		for(int j = 0; j < n; j++)
		{
			w3[j] = v2[j] - rho3 * w1[j] - rho2 * w2[j];
			w3[j] = w3[j] / rho1;
		}
		for(int j = 0; j < n; j++)
		{
			x2[j] = x1[j] + gamma3 * eta * w3[j];
//			cout << "x2[" << j << "] = " << x2[j] << "     x1[" << j << "] = " << x1[j] << endl;
		}
//		check convergence
/*		conv = 0;
		for(int j = 0; j < n; j++)
		{
			if(abs(x2[j] - x1[j]) > eps)
			{
				conv++;
			}
		}
		if(conv == 0)
		{
			cout << i << endl;
			return;
		}	*/
		r2norm = abs(sigma3) * r1norm;
		if(r2norm <= eps)
		{
			cout << i << endl;
			return;
		}
		if(i == maxiter)
		{
			cout << i << endl;
			return;
		}
		eta = - sigma3 * eta;
		
		for(int j = 0; j < n; j++)
		{
			v1[j] = v2[j];
			v2[j] = v3[j];
			w1[j] = w2[j];
			w2[j] = w3[j];
			x1[j] = x2[j];
		}
		beta1 = beta2;
		gamma1 = gamma2;
		gamma2 = gamma3;
		sigma1 = sigma2;
		sigma2 = sigma3;
		r1norm = r2norm;
	}
}

void minresforjacobidavidson(int n, int maxiter, double eps, double theta, double mu1, double *y, double *k, double *u, double *val, int *ncol, int *nrow, int *diagind, double *x2, double *b)
{
	double *x1 = new double[n];
	double *v1 = new double[n];
	double *v2 = new double[n];
	double *v3 = new double[n];
	double *Hv2 = new double[n];
	double *w1 = new double[n];
	double *w2 = new double[n];
	double *w3 = new double[n];
	double *temp = new double[n];
//	Vec_DP x1(n), v1(n), v2(n), v3(n), Hv2(n);
//	Vec_DP w1(n), w2(n), w3(n);
	double alpha, beta1 = 0.0, beta2 = 0.0;
	double gamma1 = 1.0, gamma2 = 1.0, gamma3;
	double delta, eta, sigma1 = 0.0, sigma2 = 0.0, sigma3;
	double rho1, rho2, rho3, r1norm, r2norm;
	double mu2 = 0.0;
	int conv;
	
//	initial guess for x; v_i-1 = 0; w_i-2 = w_i-1 = w_i = 0
	for(int i = 0; i < n; i++)
	{
		x1[i] = 0.0;
		v1[i] = 0.0;
		w1[i] = w2[i] = 0.0;
	}
//	mult(n,x1,v2,val,ncol,nrow);
	mult_tuti(n,x1,v2,val,ncol,nrow,diagind,theta);
	for(int i = 0; i < n; i++)
	{
		temp[i] = k[i] * v2[i];
	}
//	ILUprecond(n,val,ncol,nrow,k,diagind,temp,v2);
//	ILUprecondforjacobidavidson(n,theta,val,ncol,nrow,k,diagind,temp,v2);
	mu2 = 0.0;
	for(int i = 0; i < n; i++)
	{
		mu2 = mu2 + u[i] * temp[i];
	}
	for(int i = 0; i < n; i++)
	{
		v2[i] = temp[i] - (mu2 / mu1) * y[i];
	}
	for(int i = 0; i < n; i++)
	{
		v2[i] = b[i] - v2[i];
	}
	for(int i = 0; i < n; i++)
	{
		beta1 = beta1 + v2[i] * v2[i];
	}
	beta1 = sqrt(beta1);
	r1norm = beta1;
	eta = beta1;
	for(int i = 1; i <= maxiter; i++)
	{
//		cout << "i = " << i << endl;
		for(int j = 0; j < n; j++)
		{
			v2[j] = v2[j] / beta1;
		}
//		mult(n,v2,Hv2,val,ncol,nrow);	// Hv_i = H * v_i
		mult_tuti(n,v2,Hv2,val,ncol,nrow,diagind,theta);
		for(int j = 0; j < n; j++)
		{
			temp[j] = k[j] * Hv2[j];
		}
//		ILUprecond(n,val,ncol,nrow,k,diagind,temp,Hv2);
//		ILUprecondforjacobidavidson(n,theta,val,ncol,nrow,k,diagind,temp,Hv2);
		mu2 = 0.0;
		for(int j = 0; j < n; j++)
		{
			mu2 = mu2 + u[j] * temp[j];
		}
		for(int j = 0; j < n; j++)
		{
			Hv2[j] = temp[j] - (mu2 / mu1) * y[j];
		}
		alpha = 0.0;
		for(int j = 0; j < n; j++)
		{
			alpha = alpha + v2[j] * Hv2[j];
		}
		for(int j = 0; j < n; j++)
		{
			v3[j] = Hv2[j] - alpha * v2[j] - beta1 * v1[j];
		}
		beta2 = 0.0;
		for(int j = 0; j < n; j++)
		{
			beta2 = beta2 + v3[j] * v3[j];
		}
		beta2 = sqrt(beta2);
		
		delta = gamma2 * alpha - gamma1 * sigma2 * beta1;
		rho1 = sqrt(delta * delta + beta2 * beta2);
		rho2 = sigma2 * alpha + gamma1 * gamma2 * beta1;
		rho3 = sigma1 * beta1;
		
		gamma3 = delta / rho1;
		sigma3 = beta2 / rho1;
		
		for(int j = 0; j < n; j++)
		{
			w3[j] = v2[j] - rho3 * w1[j] - rho2 * w2[j];
			w3[j] = w3[j] / rho1;
		}
		for(int j = 0; j < n; j++)
		{
			x2[j] = x1[j] + gamma3 * eta * w3[j];
//			cout << "x2[" << j << "] = " << x2[j] << "     x1[" << j << "] = " << x1[j] << endl;
		}
//		check convergence
/*		conv = 0;
		for(int j = 0; j < n; j++)
		{
			if(abs(x2[j] - x1[j]) > eps)
			{
				conv++;
			}
		}
		if(conv == 0)
		{
			cout << i << endl;
			return;
		}	*/
		r2norm = abs(sigma3) * r1norm;
		if(r2norm <= eps)
		{
			cout << i << endl;
			return;
		}
		if(i == maxiter)
		{
//			cout << "maximum iteration number " << i << endl;
			return;
		}
		eta = - sigma3 * eta;
		
		for(int j = 0; j < n; j++)
		{
			v1[j] = v2[j];
			v2[j] = v3[j];
			w1[j] = w2[j];
			w2[j] = w3[j];
			x1[j] = x2[j];
		}
		beta1 = beta2;
		gamma1 = gamma2;
		gamma2 = gamma3;
		sigma1 = sigma2;
		sigma2 = sigma3;
		r1norm = r2norm;
	}
}

void minresfordavidson(int n, int maxiter, double eps, double *val, int *ncol, int *nrow, double **A, double *x2, double *b, int param)
{
	double *x1 = new double[n];
	double *v1 = new double[n];
	double *v2 = new double[n];
	double *v3 = new double[n];
	double *Hv2 = new double[n];
	double *w1 = new double[n];
	double *w2 = new double[n];
	double *w3 = new double[n];
	double *temp = new double[n];
//	Vec_DP x1(n), v1(n), v2(n), v3(n), Hv2(n);
//	Vec_DP w1(n), w2(n), w3(n);
//	Vec_DP temp(n);
	double alpha, beta1 = 0.0, beta2 = 0.0;
	double gamma1 = 1.0, gamma2 = 1.0, gamma3;
	double delta, eta, sigma1 = 0.0, sigma2 = 0.0, sigma3;
	double rho1, rho2, rho3;
	int conv;
	
//	initial guess for x; v_i-1 = 0; w_i-2 = w_i-1 = w_i = 0
	for(int i = 0; i < n; i++)
	{
		x1[i] = 1.0;
		v1[i] = v2[i] = 0.0;
		w1[i] = w2[i] = 0.0;
	}
	if(param == 0)
	{
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				v2[i] = v2[i] + A[i][j] * x1[j];
			}
		}
	}
/*	if(param == 1)
	{
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				v2[i] = v2[i] + Y[i][j] * x1[j];
			}
		}
//		mult(n,v2,temp,val,ncol,nrow);
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				temp[i] = temp[i] + A[i][j] * v2[j];
			}
		}
		for(int i = 0; i < n; i++)
		{
			v2[i] = 0.0;
			for(int j = 0; j < n; j++)
			{
				v2[i] = v2[i] + Y[i][j] * temp[j];
			}
		}
	}	*/
	for(int i = 0; i < n; i++)
	{
		v2[i] = b[i] - v2[i];
	}
	for(int i = 0; i < n; i++)
	{
		beta1 = beta1 + v2[i] * v2[i];
	}
	beta1 = sqrt(beta1);
	eta = beta1;
	for(int i = 1; i <= maxiter; i++)
	{
//		cout << "i = " << i << endl;
		for(int j = 0; j < n; j++)
		{
			v2[j] = v2[j] / beta1;
		}
		if(param == 0)
		{
			for(int j = 0; j < n; j++)
			{
				Hv2[j] = 0.0;
			}
			for(int j = 0; j < n; j++)
			{
				for(int k = 0; k < n; k++)
				{
					Hv2[j] = Hv2[j] + A[j][k] * v2[k];
				}
			}
		}
/*		if(param == 1)
		{
			for(int j = 0; j < n; j++)
			{
				Hv2[j] = 0.0;
				for(int k = 0; k < n; k++)
				{
					Hv2[j] = Hv2[j] + Y[j][k] * v2[k];
				}
			}
//			mult(n,Hv2,temp,val,ncol,nrow);	// Hv_i = H * v_i
			for(int j = 0; j < n; j++)
			{
				temp[j] = 0.0;
				for(int k = 0; k < n; k++)
				{
					temp[j] = temp[j] + A[j][k] * Hv2[k];
				}
			}
			for(int j = 0; j < n; j++)
			{
				Hv2[j] = 0.0;
				for(int k = 0; k < n; k++)
				{
					Hv2[j] = Hv2[j] + Y[j][k] * temp[k];
				}
			}
		}	*/
		alpha = 0.0;
		for(int j = 0; j < n; j++)
		{
			alpha = alpha + v2[j] * Hv2[j];
		}
		for(int j = 0; j < n; j++)
		{
			v3[j] = Hv2[j] - alpha * v2[j] - beta1 * v1[j];
		}
		beta2 = 0.0;
		for(int i = 0; i < n; i++)
		{
			beta2 = beta2 + v3[i] * v3[i];
		}
		beta2 = sqrt(beta2);
		
		delta = gamma2 * alpha - gamma1 * sigma2 * beta1;
		rho1 = sqrt(delta * delta + beta2 * beta2);
		rho2 = sigma2 * alpha + gamma1 * gamma2 * beta1;
		rho3 = sigma1 * beta1;
		
		gamma3 = delta / rho1;
		sigma3 = beta2 / rho1;
		
		for(int j = 0; j < n; j++)
		{
			w3[j] = v2[j] - rho3 * w1[j] - rho2 * w2[j];
			w3[j] = w3[j] / rho1;
		}
		for(int j = 0; j < n; j++)
		{
			x2[j] = x1[j] + gamma3 * eta * w3[j];
		}
//		check convergence
		conv = 0;
		for(int j = 0; j < n; j++)
		{
			if(abs(x2[j] - x1[j]) > eps)
			{
				conv++;
			}
		}
		if(conv == 0)
		{
//			cout << "minresfordavidson iter = " << i << endl;
			return;
		}
		eta = - sigma3 * eta;
		
		for(int j = 0; j < n; j++)
		{
			v1[j] = v2[j];
			v2[j] = v3[j];
			w1[j] = w2[j];
			w2[j] = w3[j];
			x1[j] = x2[j];
		}
		beta1 = beta2;
		gamma1 = gamma2;
		gamma2 = gamma3;
		sigma1 = sigma2;
		sigma2 = sigma3;
	}
	cout << "Nem konvergált! (minresfordavidson) :(" << endl;
	return;
}

//***********************************************************************
//
// axpy is subrutine, which making y(1-n) =y(1-n) + da * x(1-n)  
//
// IF IOPT(60) eq 0 then is standard BLAS rutine daxpy used 
//
//***********************************************************************
void axpy(int n, double da, double *x, double *y, ofstream &out1, int *IOPT)
{
//	local variables
	int i, m, mp1;
//	******************************************************************
	if((IOPT[60] == 0) || (n < IOPT[61]))
	{
		daxpy(n,da,x,1,y,1); //Standard BLAS routins
		return;
	}
//	******************************************************************
	if(IOPT[60] == 1)
	{
		for (i = 0; i < n; i++)
		{
			y[i] = y[i] + x[i] * da;
		}
		return;
	}
//	******************************************************************
	if(IOPT[60] == 2)
	{
		m = n % 4;
		if(m == 0)
			goto A;
		for (i = 0; i < m; i++)
		{
			y[i] = y[i] = da * x[i];
		}
		if(n < 4)
			return;
A:		mp1 = m + 1;
		for (i = mp1 - 1; i < n; i = i + 4)
		{
			y[i] = y[i] = da * x[i];
			y[i + 1] = y[i + 1] + da * x[i + 1];
			y[i + 2] = y[i + 2] + da * x[i + 2];
			y[i + 3] = y[i + 3] + da * x[i + 3];
		}
		return;
	}
//	******************************************************************
	out1 << "Error occured, wrong number in IOPT(60) !!" << endl;
	return;
}

void daxpy(int n, double da, double *dx, int incx, double *dy, int incy)
{
	int i, ix, iy, m, mp1;
	
	if(n < 0)
		return;
	if(da == 0.0)
		return;
	if((incx == 1) && (incy == 1))
	{
		m = n % 4;
		if(m != 0)
		{
			for(i = 0; i < m; i++)
			{
				dy[i] = dy[i] + da * dx[i];
			}
		}
		if(n < 4)
			return;
		mp1 = m + 1;
		for(i = (mp1 - 1); i < n; i = i + 4)
		{
			dy[i] = dy[i] + da * dx[i];
			dy[i + 1] = dy[i + 1] + (da * dx[i + 1]);
			dy[i + 2] = dy[i + 2] + (da * dx[i + 2]);
			dy[i + 3] = dy[i + 3] + (da * dx[i + 3]);
		}
	}
	else
	{
		ix = 1;
		iy = 1;
		if(incx < 0)
			ix = ((-n + 1) * incx) + 1;
		if(incy < 0)
			iy = ((-n + 1) * incy) + 1;
		for(i = 0; i < n; i++)
		{
			dy[iy] = dy[iy] + da * dx[ix];
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return;
}


//*************************************************************************
//         SAJATERTEK ES SAJATVEKTOR RENDZES
//*************************************************************************
void eigsrt(double *d, double **v, int n)
{
//	local variables
	int i, j, k;
	double p;
	for(i = 0; i < n-1; i++)
	{
		k = i;
		p = d[i];
		for(j = i; j < n; j++)
		{
			if(d[j] >= p)
			{   
				k = j;
				p = d[j];
			}
		}
		if(k != i)
		{
			d[k] = d[i];
			d[i] = p;
			for(j = 0; j < n; j++)
			{
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
	return; 
}

void eigsrtfordavidson(double *d, double **v, int n)
{
//	local variables
	int i, j, k;
	double p;
	for(i = 0; i < n-1; i++)
	{
		k = i;
		p = d[i];
		for(j = i; j < n; j++)
		{
			if(d[j] >= p)
			{   
				k = j;
				p = d[j];
			}
		}
		if(k != i)
		{
			d[k] = d[i];
			d[i] = p;
			for(j = 0; j < n; j++)
			{
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
	return; 
}

void mult(int n, double *x, double *y, double *val, int *ncol, int *nrow)
{
	int i, j;

//	x is the input
//	y is the output (y=Hx)
      
	for(i = 0; i < n; i++)
	{
		y[i] = 0.00;
		for(j = nrow[i]; j < nrow[i + 1]; j++) 
		{
			y[i] = y[i] + val[j] * x[ncol[j]];
		}
	}
	return;
}

void multfordavidson1(int n, double *x, double *y, double *val, int *ncol, int *nrow)
{
	int i, j;

//	x is the input
//	y is the output (y=Hx)
      
	for(i = 0; i < n; i++)
	{
		y[i] = 0.00;
		for(j = nrow[i]; j < nrow[i + 1]; j++) 
		{
			y[i] = y[i] + val[j] * x[ncol[j]];
		}
	}
	return;
}

void multfordavidson2(int n, int m, double **X, double **Y, double *val, int *ncol, int *nrow)
{
	int i, j;

//	x is the input
//	y is the output (y=Hx)
      
	for(i = 0; i < n; i++)
	{
		Y[i][m] = 0.00;
		for(j = nrow[i]; j < nrow[i + 1]; j++) 
		{
			Y[i][m] = Y[i][m] + val[j] * X[ncol[j]][m];
		}
	}
	return;
}

void multfordavidson3(int n, int m, double **X, double **Y, double *val, int *ncol, int *nrow)
{
	int i, j, k;

//	x is the input
//	y is the output (y=Hx)

	for(k = 0; k < m; k++)
	{
		for(i = 0; i < n; i++)
		{
			Y[i][k] = 0.00;
			for(j = nrow[i]; j < nrow[i + 1]; j++) 
			{
				Y[i][k] = Y[i][k] + val[j] * X[ncol[j]][k];
			}
		}
	}
	return;
}

void mult2(int n, double *x, double *y, double *val, int *ncol, int *nrow)
{
	int i, j;
	
	for(i = 0; i < n; i++)
	{
		y[i] = 0.00;
	}
	for(j = 0; j < n; j++)
	{
		for(i = nrow[j]; i < nrow[j + 1]; i++)
		{
			y[ncol[i]] = y[ncol[i]] + val[i] * x[j];
		}
	}
	return;
}

void mult_tuti(int n, double *x, double *y, double *val, int *ncol, int *nrow, int *diagind, double theta)
{
	int i, j, k = 0;
	
	for(i = 0; i < n; i++)
	{
		y[i] = 0.00;
		for(j = nrow[i]; j < nrow[i + 1]; j++) 
		{
			if(j == diagind[k])
			{
				y[i] = y[i] + (val[j] - theta) * x[ncol[j]];
				k++;
			}
			else
			{
				y[i] = y[i] + val[j] * x[ncol[j]];
			}
		}
	}
	return;
}

void transform_tuti(int q, double *val, double *val2, int *ncol, int *nrow, double *ud, int *diagind, double theta)
{
	int i, j, k = 0;

	for(i = 0; i < q; i++)
	{
		if(i == diagind[k])
		{
//			cout << "diagind = " << diagind[k] << endl;
			val2[i] = val[i] - theta;
			ud[k] = (1.0 / val2[i]);
			k++;
		}
		else
		{
			val2[i] = val[i];
		}
	}
	return;
}

void transform(int q, double *val, double *val2, int *ncol, int *nrow, double *ud, int *diagind, double theta)
{
	int i, j, k = 0;

	for(i = 0; i < q; i++)
	{
		if(i == diagind[k])
		{
//			cout << "diagind = " << diagind[k] << endl;
			val2[i] = val[i] - theta;
			ud[k] = 1.0 / val2[i];
			k++;
		}
		else
		{
			val2[i] = val[i];
		}
	}
	return;
}

void transform2(int n, double *diag, double *diagtemp, double theta)
{
	for(int i = 0; i < n; i++)
	{
		diagtemp[i] = diag[i] - theta;
	}
	for(int i = 0; i < n; i++)
	{
		diagtemp[i] = (1.0 / diagtemp[i]);
	}
}

void MinvA(int n, double *val, double *val2, int *ncol, int *nrow, double *diag)
{
	int i, j;
      
	for(i = 0; i < n; i++)
	{
		for(j = nrow[i]; j < nrow[i + 1]; j++) 
		{
			val2[j] = (1.0 / diag[i]) * val[j];
		}
	}
	return;
}


//***********************************************************************
//
// norm2 is subrutine, which making norma2 =dsqrt( sum(   v(1-n) * v(1-n) )  )  
//
// IF IOPT(60) eq 0 then is standard BLAS function dnrm2 used
//
//***********************************************************************
double norm2(int n, double norma2, double *v, int *IOPT)
{
//	local variables
	double vxv;
	int i;

//**********************************************************************
	if((IOPT[60] == 0) || (n < IOPT[61]))
	{
		norma2 = dnrm2(n,v,1); // Standard BLAS routins
		return norma2;
	}
//**********************************************************************
	if((IOPT[60] == 1) || (IOPT[60] == 2))
	{
		norma2 = 0.0;
		vxv    = 0.0;
		for(i = 0; i < n; i++)
		{
			vxv = vxv + v[i] * v[i];
		}
		norma2 = sqrt(vxv);
		return norma2;
	}
//**********************************************************************
	return norma2;
}

double dnrm2(int n, double *x, int incx)
{
	const double one = 1.0;
	const double zero = 0.0;
	double absxi, norm, scale, ssq;
	int ix;
	
	if((n < 1) || (incx < 1))
		norm = zero;
	else if(n == 1)
		norm = abs(x[0]);
	else
	{
		scale = zero;
		ssq = one;
		for(ix = 0; ix < (1+(n-1)*incx); ix = ix + incx)
		{
			if(x[ix] != zero)
			{
				absxi = abs(x[ix]);
				if(scale < absxi)
				{
					ssq = one + ssq * (scale/absxi) * (scale/absxi);
					scale = absxi;
				}
				else
				{
					ssq = ssq + (absxi/scale) * (absxi/scale);
				}
			}
		}
		norm = scale * sqrt(ssq);
//		cout << "norm = " << norm << endl;
	}
	return norm;
}

//***********************************************************************
//
// Scal is subrutine, which making x(1-n) = da * x(1-n)  
//
// IF IOPT(60) eq 0 then is standard BLAS rutine dscal used
//
//***********************************************************************
void scal(int n, double da, double *x, int *IOPT)
{
//	local variables
	int i,m,mp1;
//**********************************************************************
	if ((IOPT[60] == 0) || (n < IOPT[61]))
	{
		dscal(n,da,x,1); // Standard BLAS1 routins
		return;
	}
//**********************************************************************
	if(IOPT[60] == 1)
	{
		for(i = 0; i < n; i++)
		{
			x[i] =  x[i] * da;
		}
	return;
	}
//**********************************************************************
	if(IOPT[60] == 2)
	{
		m = n % 5;
		if(m == 0)
			goto A;
		for(i = 0; i < m; i++)
		{
			x[i] = da * x[i];
		}
		if(n < 5)
			return;
A:		mp1 = m + 1;
		for(i = mp1 -1; i < n; i = i + 5)
		{
			x[i] = da * x[i];
			x[i + 1] = da * x[i + 1];
			x[i + 2] = da * x[i + 2];
			x[i + 3] = da * x[i + 3];
			x[i + 4] = da * x[i + 4];
		}
		return;
	}
	return;
}

void dscal(int n, double da, double *dx, int incx)
{
	int i,m,mp1,nincx;
	
	if((n < 0) || (incx < 0))
		return;
	if(incx == 1)
	{
		m = n % 5;
		if(m != 0)
		{
			for(i = 0; i < m; i++)
			{
				dx[i] = da * dx[i];
			}
			if(n < 5)
				return;
		}
		mp1 = m + 1;
		for(i = mp1 - 1; i < n; i = i + 5)
		{
			dx[i] = da * dx[i];
			dx[i + 1] = da * dx[i + 1];
			dx[i + 2] = da * dx[i + 2];
			dx[i + 3] = da * dx[i + 3];
			dx[i + 4] = da * dx[i + 4];
		}
	}
	else
	{
		nincx = n * incx;
		for(i = 0; i < nincx; i = i + incx)
		{
			dx[i] = da * dx[i];
		}
	}
	return;
}

//***********************************************************************
//
// proda is subrutine, which making prod =sum(  x(1-n) * y(1-n)  )
//
// IF IOPT[60] eq 0 then is standard BLAS function ddot used
//
//***********************************************************************
double proda(int n, double *x, double *y, double prod, int *IOPT)
{
//	local variables
	int i, m, mp1;
	double szum;
//**********************************************************************
	if((IOPT[60] == 0) || (n < IOPT[61]))
	{
		prod = ddot(n,x,1,y,1); // Standard BLAS routins
		return prod;
	}
//**********************************************************************
	if(IOPT[60] == 1)
	{
		szum=0.0;
		for(i = 0; i < n; i++)
		{
			szum = szum + x[i] * y[i];
		}
		prod = szum;
		return prod;
	}
//**********************************************************************
	if(IOPT[60] == 2)
	{  
		szum=0.0;
		m = n % 5;
		if(m == 0)
			goto A;
		for(i = 0; i < m; i++)
		{
			szum = szum + x[i]*y[i];
		}
		if(n < 5)
			goto B;
A:		mp1 = m + 1;
		for(i = mp1 - 1; i < n; i = i + 5)
		{
			szum=szum+x[i]*y[i]+x[i+1]*y[i+1]+x[i+2]*y[i+2]+x[i+3]*y[i+3]+x[i+4]*y[i+4];
		}
B:		prod = szum;
		return prod;
	}
//**********************************************************************
	return prod;
}

double ddot(int n, double *dx, int incx, double *dy, int incy)
{
	double dtemp = 0.0;
	int i, ix, iy, m, mp1;
	
	if(n < 0)
		return dtemp;
	if((incx == 1) && (incy == 1))
	{
		m = n % 5;
		if(m != 0)
		{
			for(i = 0; i < m; i++)
			{
				dtemp = dtemp + dx[i] * dy[i];
			}
			if(n < 5)
				return dtemp;
		}
		mp1 = m + 1;
		for(i = mp1 - 1; i < n; i = i + 5)
		{
			dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2] * dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i + 4];
		}
	}
	else
	{
		ix = 1;
		iy = 1;
		if(incx < 0)
			ix = (-n + 1) * incx + 1;
		if(incy < 0)
			iy = (-n + 1) * incy + 1;
		for(i = 0; i < n; i++)
		{
			dtemp = dtemp + dx[ix] * dy[iy];
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return dtemp;
}

void trialbasis(int n, int bsize, double **G, double *val, int *ncol, int *nrow)
{
	int l = 0;
	double alpha;
	
	for(int i = 0; i < bsize; i++)
	{
		for(int j = 0; j < n; j++) 
		{
			if(j == ncol[l])
			{
				G[j][i] = val[l];
				l++;
			}
			else
			{
				G[j][i] = 0.00;
			}
			if(j > 0 && ncol[l+1] < ncol[l])
				continue;
		}
	}
	double norm;
	for(int i = 0; i < bsize; i++)
	{
		norm = 0.00;
		for(int j = 0; j < n; j++) 
		{
			norm = norm + G[j][i] * G[j][i];
		}
		norm = sqrt(norm);
		for(int j = 0; j < n; j++) 
		{
			G[j][i] = (1.0 / norm) * G[j][i];
		}
	}
	for(int k = bsize - 1; k >= 0; k--)
	{
		for(int i = 0; i < bsize; i++)
		{
			if(k != i)
			{
				alpha = 0.00;
				for(int j = 0; j < n; j++)
				{
					alpha = alpha + G[j][k] * G[j][i];
				}
				alpha = -1.0 * alpha;
				for(int j = 0; j < n; j++)
				{
					G[j][k] = G[j][k] + alpha * G[j][i];
				}
			}
		}
	}
}

void comp2mtx(int nt, int nonzero, double **A, double *val, double *diag, double *ud, int *ncol, int *nrow)
{
	int i, j, k, l;
	i = j = k = l = 0;
	ofstream mtx;
	mtx.setf(ios_base::fixed,ios_base::floatfield);
	
	while(k < nonzero)
	{
		if(k + 1 < nonzero)
		{
			if(ncol[k] < ncol[k+1])
			{
				j = ncol[k];
				A[i][j] = val[k];
				j++;
			}
			else
			{
				j = ncol[k];
				A[i][j] = val[k];
				i++;
			}
		}
		k++;
	}
	A[nt-1][nt-1] = val[nonzero-1];
	
/*	for(i = 0; i < nt; i++)
	{
		for(j = 0; j < nt; j++)
		{
			if(i == j)
				A[i][j] = A[i][j];
			else if(i == j + 1 || i == j - 1)
				A[i][j] = A[i][j];
			else if(i == j + 2 || i == j - 2)
				A[i][j] = A[i][j];
			else
				A[i][j] = 0.00;
		}
	}	*/
	
//	cout << "A = " << endl;
	mtx.open("mtx.dat");
	for(i = 0; i < nt; i++)
	{
		for(j = 0; j < nt; j++)
		{
//			cout << A[i][j] << " ";
			mtx << setprecision(20) << A[i][j] << " ";
		}
//		cout << endl;
		mtx << endl;
	}
	mtx.close();
	return;
}

void mtx2comp(int nt, int nonzero, double **A, double *val, double *diag, int *diagind, double *ud, int *ncol, int *nrow)
{
	int i, j, k, l, m;
	i = j = k = l = m = 0;
	bool w = true;
	int sum = 0;
	ofstream mtx;
	mtx.open("mtx.dat");
	
	cout << "A = " << endl;
	for(i = 0; i < nt; i++)
	{
		for(j = 0; j < nt; j++)
		{
			mtx << A[i][j] << " ";
		}
		mtx << endl;
	}
	for(i = 0; i < nt; i++)
	{
		diag[i] = A[i][i];
	}
	for(i = 0; i < nt; i++)
	{
		w = true;
		for(j = 0; j < nt; j++)
		{
			if(A[i][j] != 0.00)
			{
				val[k] = A[i][j];
				if(i == j)
				{
					diagind[m] = k;
					m++;
				}
				ncol[k] = j;
				if(w == true)
				{
					nrow[l] = k;
					l++;
					w = false;
				}
				k++;
			}
		}
	}
	for(i = 0; i < nt; i++)
	{
		if(A[nt-1][i] != 0.0)
			sum++;
	}
	nrow[nt] = nrow[nt-1] + sum;
	for(i = 0; i < nt; i++)
	{
		if(abs(diag[i]) >= 1.0E-10)
		{
			ud[i] = 1.0 / diag[i];
		}
		else
		{
			ud[i] = 1.0;
//			ud[i] = 1.0E10;
		}
	}
		
	ofstream d;
	d.open("diag.txt");
	for(i = 0; i < nt; i++)
	{
		cout << setprecision(20) << "diag[" << i << "] = " << diag[i] << endl;
//		d << setprecision(20) << diag[i] << endl;
	}
	ofstream u;
	u.open("ud.txt");
	for(i = 0; i < nt; i++)
	{
		cout << setprecision(20) << "ud[" << i << "] = " << ud[i] << endl;
//		u << setprecision(20) << ud[i] << endl;
	}
	ofstream v;
	v.open("val.txt");
	for(i = 0; i < nonzero; i++)
	{
		cout << setprecision(20) << "val[" << i << "] = " << val[i] << endl;
//		v << setprecision(20) << val[i] << endl;
	}
	ofstream nc;
	nc.open("ncol.txt");
	for(i = 0; i < nonzero; i++)
	{
		cout << setprecision(20) << "ncol[" << i << "] = " << ncol[i] << endl;
//		nc << ncol[i] << endl;
	}
	ofstream nr;
	nr.open("nrow.txt");
	for(i = 0; i < nt+1; i++)
	{
		cout << setprecision(20) << "nrow[" << i << "] = " << nrow[i] << endl;
//		nr << nrow[i] << endl;
	}
	d.close();
	u.close();
	v.close();
	nc.close();
	nr.close();
		
	cout<<"DONE!"<<endl;
	return;
}

void expectval(int nmax, int ns1, int np, int nt, double Re, double Rmax, double tm2atom, double tm3atom, string eivectype)
{
	double var1, var2, var3;
	double *allgrid = new double[ns1+np];
	double *de = new double[ns1];
	double *q = new double[np];
	double *R1 = new double[nmax];
	double *R1_2 = new double[nmax];
	double *R1_3 = new double[nmax];
	double *R1_M1 = new double[nmax];
	double *R1_M2 = new double[nmax];
	double *R1_M3 = new double[nmax];
	double *R2 = new double[nmax];
	double *R2_2 = new double[nmax];
	double *R2_3 = new double[nmax];
	double *R2_M1 = new double[nmax];
	double *R2_M2 = new double[nmax];
	double *R2_M3 = new double[nmax];
	double *Theta = new double[nmax];
	double **svmat = new double *[nt];
	for(int i = 0; i < nt; i++)
		svmat[i] = new double[nmax];
/*	Vec_DP allgrid(ns1+np), de(ns1), q(np);
	Vec_DP R1(nmax), R1_2(nmax), R1_3(nmax), R1_M1(nmax), R1_M2(nmax), R1_M3(nmax);
	Vec_DP R2(nmax), R2_2(nmax), R2_3(nmax), R2_M1(nmax), R2_M2(nmax), R2_M3(nmax);
	Vec_DP Theta(nmax);
	Mat_DP svmat(nt,nmax);	*/
//	******várható érték számítás************
//	svmat olvasás fájlból
	ifstream eivec;
	if(eivectype == "Lanczos")
		eivec.open ("LEIVEC.txt");
	if(eivectype == "Davidson")
		eivec.open ("DEIVEC.txt");
	if(eivectype == "JacobiDavidson")
		eivec.open ("JDEIVEC.txt");
	if(eivectype == "RQ-iteration")
		eivec.open ("RQEIVEC.txt");
	if(eivectype == "Inverse-Power-iteration")
		eivec.open ("INVPOWEIVEC.txt");
	if(eivectype == "Power-iteration")
		eivec.open ("POWEIVEC.txt");
	eivec.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < nt; i++)
	{
		for(int j = 0; j < nmax; j++)
		{
			eivec >> svmat[i][j];
		}
	}
	eivec.close();
	cout << eivectype << " állapotfüggvényre várható érték számítás kezdődik" << endl;
//	grid pontok olvasása fájlból
	ifstream grid;
	grid.open ("gridfile.txt");
	grid.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < ns1+np; i++)
	{
		grid >> allgrid[i];
	}
	grid.close();
	for(int i = 0; i < ns1; i++)
	{
		de[i] = allgrid[i];
	}
	for(int i = 0; i < np; i++)
	{
		q[i] = allgrid[i+ns1];
	}
	double allk = de[ns1-1]/Rmax;
	for(int i = 0; i < ns1; i++)
	{
		de[i] = de[i] / allk + Re;
	}
	for(int p = 0; p < nmax; p++)
	{
		R1[p] = 0.00;
		R1_2[p] = 0.00;
		R1_3[p] = 0.00;
		R1_M1[p] = 0.00;
		R1_M2[p] = 0.00;
		R1_M3[p] = 0.00;
		
		R2[p] = 0.00;
		R2_2[p] = 0.00;
		R2_3[p] = 0.00;
		R2_M1[p] = 0.00;
		R2_M2[p] = 0.00;
		R2_M3[p] = 0.00;
		Theta[p] = 0.00;
	}
	int l = 0;
	for(int p = 0; p < nmax; p++)
	{
		for(int i = 0; i < ns1; i++)
		{
			for(int j = 0; j < ns1; j++)
			{
				for(int k = 0; k < np; k++)
				{
					jacbolbelso2(tm3atom,tm2atom,de[i],de[j],q[k],var1,var2,var3);
					R1[p] = R1[p] + svmat[l][p] * var1 * svmat[l][p];
					R1_2[p] = R1_2[p] + svmat[l][p] * (var1 * var1) * svmat[l][p];
					R1_3[p] = R1_3[p] + svmat[l][p] * (var1 * var1 * var1) * svmat[l][p];
					R1_M1[p] = R1_M1[p] + svmat[l][p] * pow(var1,-1.0) * svmat[l][p];
					R1_M2[p] = R1_M2[p] + svmat[l][p] * pow(var1,-2.0) * svmat[l][p];
					R1_M3[p] = R1_M3[p] + svmat[l][p] * pow(var1,-3.0) * svmat[l][p];
					
					R2[p] = R2[p] + svmat[l][p] * var2 * svmat[l][p];
					R2_2[p] = R2_2[p] + svmat[l][p] * (var2 * var2) * svmat[l][p];
					R2_3[p] = R2_3[p] + svmat[l][p] * (var2 * var2 * var2) * svmat[l][p];
					R2_M1[p] = R2_M1[p] + svmat[l][p] * pow(var2,-1.0) * svmat[l][p];
					R2_M2[p] = R2_M2[p] + svmat[l][p] * pow(var2,-2.0) * svmat[l][p];
					R2_M3[p] = R2_M3[p] + svmat[l][p] * pow(var2,-3.0) * svmat[l][p];
					
					Theta[p] = Theta[p] + svmat[l][p] * acos(var3) * svmat[l][p];
					if(l != nt - 1)
					{
						l++;
					}
					else
					{
						l = 0;
					}
				}
			}
		}
	}
	for(int p = 0; p < nmax; p++)
	{
		cout << "R1[" << p << "] = " << setprecision(5) << R1[p] * 0.529177249 <<
		" || sqrt(R1_2[" << p << "]) = " << sqrt(R1_2[p]) * 0.529177249 <<
		" || sqrt(R1_3[" << p << "]) = " << pow(R1_3[p],1.0/3.0) * 0.529177249 <<
		" || sqrt(R1_M1[" << p << "]) = " << pow(R1_M1[p],-1.0) * 0.529177249 <<
		" || sqrt(R1_M2[" << p << "]) = " << pow(R1_M2[p],(-1.0/2.0)) * 0.529177249 <<
		" || sqrt(R1_M3[" << p << "]) = " << pow(R1_M3[p],(-1.0/3.0)) * 0.529177249 << endl;
		
		cout << "R2[" << p << "] = " << setprecision(5) << R2[p] * 0.529177249 <<
		" || sqrt(R2_2[" << p << "]) = " << sqrt(R2_2[p]) * 0.529177249 <<
		" || sqrt(R2_3[" << p << "]) = " << pow(R2_3[p],1.0/3.0) * 0.529177249 <<
		" || sqrt(R2_M1[" << p << "]) = " << pow(R2_M1[p],-1.0) * 0.529177249 <<
		" || sqrt(R2_M2[" << p << "]) = " << pow(R2_M2[p],(-1.0/2.0)) * 0.529177249 <<
		" || sqrt(R2_M3[" << p << "]) = " << pow(R2_M3[p],(-1.0/3.0)) * 0.529177249 << endl;
		
		cout << "Theta[" << p << "] = " << setprecision(6) << Theta[p] * 180.0/3.14159265359 << endl;
	}
	cout << eivectype << " állapotfüggvényre várható érték számítás vége" << endl;
}

void expectval_rot(int nmax, int ns1, int np, int jmin, int jmax, int nt, double Re, double Rmax, double tm2atom, double tm3atom, string eivectype)
{
	double var1, var2, var3;
	double *allgrid = new double[ns1+np];
	double *de = new double[ns1];
	double *q = new double[np];
	double *R1 = new double[nmax];
	double *R1_2 = new double[nmax];
	double *R1_3 = new double[nmax];
	double *R1_M1 = new double[nmax];
	double *R1_M2 = new double[nmax];
	double *R1_M3 = new double[nmax];
	double *R2 = new double[nmax];
	double *R2_2 = new double[nmax];
	double *R2_3 = new double[nmax];
	double *R2_M1 = new double[nmax];
	double *R2_M2 = new double[nmax];
	double *R2_M3 = new double[nmax];
	double *Theta = new double[nmax];
	double **svmat = new double *[nt];
	for(int i = 0; i < nt; i++)
		svmat[i] = new double[nmax];
/*	Vec_DP allgrid(ns1+np), de(ns1), q(np);
	Vec_DP R1(nmax), R1_2(nmax), R1_3(nmax), R1_M1(nmax), R1_M2(nmax), R1_M3(nmax);
	Vec_DP R2(nmax), R2_2(nmax), R2_3(nmax), R2_M1(nmax), R2_M2(nmax), R2_M3(nmax);
	Vec_DP Theta(nmax);
	Mat_DP svmat(nt,nmax);	*/
//	******várható érték számítás************
//	svmat olvasás fájlból
	ifstream eivec;
	if(eivectype == "Lanczos")
		eivec.open ("LEIVEC.txt");
	if(eivectype == "Davidson")
		eivec.open ("DEIVEC.txt");
	if(eivectype == "JacobiDavidson")
		eivec.open ("JDEIVEC.txt");
	if(eivectype == "RQ-iteration")
		eivec.open ("RQEIVEC.txt");
	if(eivectype == "Inverse-Power-iteration")
		eivec.open ("INVPOWEIVEC.txt");
	if(eivectype == "Power-iteration")
		eivec.open ("POWEIVEC.txt");
	eivec.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < nt; i++)
	{
		for(int j = 0; j < nmax; j++)
		{
			eivec >> svmat[i][j];
		}
	}
	eivec.close();
	cout << eivectype << " állapotfüggvényre várható érték számítás kezdődik" << endl;
//	grid pontok olvasása fájlból
	ifstream grid;
	grid.open ("gridfile.txt");
	grid.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < ns1+np; i++)
	{
		grid >> allgrid[i];
	}
	grid.close();
	for(int i = 0; i < ns1; i++)
	{
		de[i] = allgrid[i];
	}
	for(int i = 0; i < np; i++)
	{
		q[i] = allgrid[i+ns1];
	}
	double allk = de[ns1-1]/Rmax;
	for(int i = 0; i < ns1; i++)
	{
		de[i] = de[i] / allk + Re;
	}
	for(int p = 0; p < nmax; p++)
	{
		R1[p] = 0.00;
		R1_2[p] = 0.00;
		R1_3[p] = 0.00;
		R1_M1[p] = 0.00;
		R1_M2[p] = 0.00;
		R1_M3[p] = 0.00;
		
		R2[p] = 0.00;
		R2_2[p] = 0.00;
		R2_3[p] = 0.00;
		R2_M1[p] = 0.00;
		R2_M2[p] = 0.00;
		R2_M3[p] = 0.00;
		Theta[p] = 0.00;
	}
	int l = 0;
	for(int p = 0; p < nmax; p++)
	{
		for(int km = jmin; km <= jmax; km++)
		{
			for(int km1 = jmin; km1 <= jmax; km1++)
			{
				if(km == km1)
				{
					for(int i = 0; i < ns1; i++)
					{
						for(int j = 0; j < ns1; j++)
						{
							for(int k = 0; k < np; k++)
							{
								jacbolbelso2(tm3atom,tm2atom,de[i],de[j],q[k],var1,var2,var3);
								R1[p] = R1[p] + svmat[l][p] * var1 * svmat[l][p];
								R1_2[p] = R1_2[p] + svmat[l][p] * (var1 * var1) * svmat[l][p];
								R1_3[p] = R1_3[p] + svmat[l][p] * (var1 * var1 * var1) * svmat[l][p];
								R1_M1[p] = R1_M1[p] + svmat[l][p] * pow(var1,-1.0) * svmat[l][p];
								R1_M2[p] = R1_M2[p] + svmat[l][p] * pow(var1,-2.0) * svmat[l][p];
								R1_M3[p] = R1_M3[p] + svmat[l][p] * pow(var1,-3.0) * svmat[l][p];
								
								R2[p] = R2[p] + svmat[l][p] * var2 * svmat[l][p];
								R2_2[p] = R2_2[p] + svmat[l][p] * (var2 * var2) * svmat[l][p];
								R2_3[p] = R2_3[p] + svmat[l][p] * (var2 * var2 * var2) * svmat[l][p];
								R2_M1[p] = R2_M1[p] + svmat[l][p] * pow(var2,-1.0) * svmat[l][p];
								R2_M2[p] = R2_M2[p] + svmat[l][p] * pow(var2,-2.0) * svmat[l][p];
								R2_M3[p] = R2_M3[p] + svmat[l][p] * pow(var2,-3.0) * svmat[l][p];
								
								Theta[p] = Theta[p] + svmat[l][p] * acos(var3) * svmat[l][p];
								if(l != nt - 1)
								{
									l++;
								}
								else
								{
									l = 0;
								}
							}
						}
					}
				}
			}
		}
	}
	for(int p = 0; p < nmax; p++)
	{
		cout << "R1[" << p << "] = " << setprecision(5) << R1[p] * 0.529177249 <<
		" || sqrt(R1_2[" << p << "]) = " << sqrt(R1_2[p]) * 0.529177249 <<
		" || sqrt(R1_3[" << p << "]) = " << pow(R1_3[p],1.0/3.0) * 0.529177249 <<
		" || sqrt(R1_M1[" << p << "]) = " << pow(R1_M1[p],-1.0) * 0.529177249 <<
		" || sqrt(R1_M2[" << p << "]) = " << pow(R1_M2[p],(-1.0/2.0)) * 0.529177249 <<
		" || sqrt(R1_M3[" << p << "]) = " << pow(R1_M3[p],(-1.0/3.0)) * 0.529177249 << endl;
		
		cout << "R2[" << p << "] = " << setprecision(5) << R2[p] * 0.529177249 <<
		" || sqrt(R2_2[" << p << "]) = " << sqrt(R2_2[p]) * 0.529177249 <<
		" || sqrt(R2_3[" << p << "]) = " << pow(R2_3[p],1.0/3.0) * 0.529177249 <<
		" || sqrt(R2_M1[" << p << "]) = " << pow(R2_M1[p],-1.0) * 0.529177249 <<
		" || sqrt(R2_M2[" << p << "]) = " << pow(R2_M2[p],(-1.0/2.0)) * 0.529177249 <<
		" || sqrt(R2_M3[" << p << "]) = " << pow(R2_M3[p],(-1.0/3.0)) * 0.529177249 << endl;
		
		cout << "Theta[" << p << "] = " << setprecision(6) << Theta[p] * 180.0/3.14159265359 << endl;
	}
	cout << eivectype << " állapotfüggvényre várható érték számítás vége" << endl;
/*	for(int p = 0; p < nmax; p++)
	{
		cout << fixed << setprecision(5) << R1[p] * 0.529177249 << " & " << sqrt(R1_2[p]) * 0.529177249 << " & " << pow(R1_3[p],1.0/3.0) * 0.529177249 << " & " << pow(R1_M1[p],-1.0) * 0.529177249 << " & " << pow(R1_M2[p],(-1.0/2.0)) * 0.529177249 << " & " << pow(R1_M3[p],(-1.0/3.0)) * 0.529177249 << endl;
		cout << fixed << setprecision(5) << R2[p] * 0.529177249 << " & " << sqrt(R2_2[p]) * 0.529177249 << " & " << pow(R2_3[p],1.0/3.0) * 0.529177249 << " & " << pow(R2_M1[p],-1.0) * 0.529177249 << " & " << pow(R2_M2[p],(-1.0/2.0)) * 0.529177249 << " & " << pow(R2_M3[p],(-1.0/3.0)) * 0.529177249 << endl;
		cout << fixed << setprecision(3) << Theta[p] * 180.0/3.14159265359 << endl;
		cout << endl;
	}	*/
	ofstream expectvalues;
	expectvalues.open ("expectvalues.txt");
	expectvalues.setf(ios_base::fixed,ios_base::floatfield);
	expectvalues << setprecision(10) << R1[0] << endl;
	expectvalues << setprecision(10) << R2[0] << endl;
	expectvalues << setprecision(10) << Theta[0] << endl;
	expectvalues.close();
}

void varhatoertektolelteres(int nmax, int ns1, int np, int jmin, int jmax, int nt, double Re, double Rmax, double tm2atom, double tm3atom, string eivectype)
{
	double var1, var2, var3;
	double *allgrid = new double[ns1+np];
	double *expvals = new double[3];
	double *de = new double[ns1];
	double *q = new double[np];
	double *R1 = new double[nmax];
	double *R1_2 = new double[nmax];
	double *R1_3 = new double[nmax];
	double *R1_M1 = new double[nmax];
	double *R1_M2 = new double[nmax];
	double *R1_M3 = new double[nmax];
	double *R2 = new double[nmax];
	double *R2_2 = new double[nmax];
	double *R2_3 = new double[nmax];
	double *R2_M1 = new double[nmax];
	double *R2_M2 = new double[nmax];
	double *R2_M3 = new double[nmax];
	double *Theta = new double[nmax];
	double *Theta2 = new double[nmax];
	double **svmat = new double *[nt];
	for(int i = 0; i < nt; i++)
		svmat[i] = new double[nmax];
/*	Vec_DP allgrid(ns1+np), de(ns1), q(np);
	Vec_DP R1(nmax), R1_2(nmax), R1_3(nmax), R1_M1(nmax), R1_M2(nmax), R1_M3(nmax);
	Vec_DP R2(nmax), R2_2(nmax), R2_3(nmax), R2_M1(nmax), R2_M2(nmax), R2_M3(nmax);
	Vec_DP Theta(nmax);
	Mat_DP svmat(nt,nmax);	*/
//	******várható érték számítás************
//	svmat olvasás fájlból
	ifstream eivec;
	if(eivectype == "Lanczos")
		eivec.open ("LEIVEC.txt");
	if(eivectype == "Davidson")
		eivec.open ("DEIVEC.txt");
//	if(eivectype == "JacobiDavidson")
		eivec.open ("JDEIVEC.txt");
	if(eivectype == "RQ-iteration")
		eivec.open ("RQEIVEC.txt");
	if(eivectype == "Inverse-Power-iteration")
		eivec.open ("INVPOWEIVEC.txt");
	if(eivectype == "Power-iteration")
		eivec.open ("POWEIVEC.txt");
	eivec.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < nt; i++)
	{
		for(int j = 0; j < nmax; j++)
		{
			eivec >> svmat[i][j];
		}
	}
	eivec.close();
//	cout << eivectype << " állapotfüggvényre várható érték számítás kezdődik" << endl;
//	grid pontok olvasása fájlból
	ifstream grid;
	grid.open ("gridfile.txt");
	grid.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < ns1+np; i++)
	{
		grid >> allgrid[i];
	}
	grid.close();
	ifstream expectvalues;
	expectvalues.open ("expectvalues.txt");
	grid.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < 3; i++)
	{
		expectvalues >> expvals[i];
	}
	expectvalues.close();
	for(int i = 0; i < ns1; i++)
	{
		de[i] = allgrid[i];
	}
	for(int i = 0; i < np; i++)
	{
		q[i] = allgrid[i+ns1];
	}
	double allk = de[ns1-1]/Rmax;
	for(int i = 0; i < ns1; i++)
	{
		de[i] = de[i] / allk + Re;
	}
	for(int p = 0; p < nmax; p++)
	{
		R1[p] = 0.00;
		R1_2[p] = 0.00;
		R1_3[p] = 0.00;
		R1_M1[p] = 0.00;
		R1_M2[p] = 0.00;
		R1_M3[p] = 0.00;
		
		R2[p] = 0.00;
		R2_2[p] = 0.00;
		R2_3[p] = 0.00;
		R2_M1[p] = 0.00;
		R2_M2[p] = 0.00;
		R2_M3[p] = 0.00;
		Theta[p] = 0.00;
	}
	int l = 0;
	for(int p = 0; p < nmax; p++)
	{
		for(int km = jmin; km <= jmax; km++)
		{
			for(int km1 = jmin; km1 <= jmax; km1++)
			{
				if(km == km1)
				{
					for(int i = 0; i < ns1; i++)
					{
						for(int j = 0; j < ns1; j++)
						{
							for(int k = 0; k < np; k++)
							{
								jacbolbelso2(tm3atom,tm2atom,de[i],de[j],q[k],var1,var2,var3);
								R1[p] = R1[p] + svmat[l][p] * (var1 - expvals[0]) * svmat[l][p];
								R1_2[p] = R1_2[p] + svmat[l][p] * (var1 - expvals[0]) * (var1 - expvals[0]) * svmat[l][p];
								R1_3[p] = R1_3[p] + svmat[l][p] * (var1 * var1 * var1) * svmat[l][p];
								R1_M1[p] = R1_M1[p] + svmat[l][p] * pow(var1,-1.0) * svmat[l][p];
								R1_M2[p] = R1_M2[p] + svmat[l][p] * pow(var1,-2.0) * svmat[l][p];
								R1_M3[p] = R1_M3[p] + svmat[l][p] * pow(var1,-3.0) * svmat[l][p];
								
								R2[p] = R2[p] + svmat[l][p] * (var2 - expvals[1]) * svmat[l][p];
								R2_2[p] = R2_2[p] + svmat[l][p] * (var2 - expvals[1]) * (var2 - expvals[1]) * svmat[l][p];
								R2_3[p] = R2_3[p] + svmat[l][p] * (var2 * var2 * var2) * svmat[l][p];
								R2_M1[p] = R2_M1[p] + svmat[l][p] * pow(var2,-1.0) * svmat[l][p];
								R2_M2[p] = R2_M2[p] + svmat[l][p] * pow(var2,-2.0) * svmat[l][p];
								R2_M3[p] = R2_M3[p] + svmat[l][p] * pow(var2,-3.0) * svmat[l][p];
								
								Theta[p] = Theta[p] + svmat[l][p] * (acos(var3) - expvals[2]) * svmat[l][p];
								Theta2[p] = Theta2[p] + svmat[l][p] * (acos(var3) - expvals[2]) * (acos(var3) - expvals[2]) * svmat[l][p];
								if(l != nt - 1)
								{
									l++;
								}
								else
								{
									l = 0;
								}
							}
						}
					}
				}
			}
		}
	}
/*	for(int p = 0; p < nmax; p++)
	{
		cout << "R1[" << p << "] = " << setprecision(5) << R1[p] * 0.529177249 <<
		" || sqrt(R1_2[" << p << "]) = " << R1_2[p] * 0.529177249 <<
		" || sqrt(R1_3[" << p << "]) = " << pow(R1_3[p],1.0/3.0) * 0.529177249 <<
		" || sqrt(R1_M1[" << p << "]) = " << pow(R1_M1[p],-1.0) * 0.529177249 <<
		" || sqrt(R1_M2[" << p << "]) = " << pow(R1_M2[p],(-1.0/2.0)) * 0.529177249 <<
		" || sqrt(R1_M3[" << p << "]) = " << pow(R1_M3[p],(-1.0/3.0)) * 0.529177249 << endl;
		
		cout << "R2[" << p << "] = " << setprecision(5) << R2[p] * 0.529177249 <<
		" || sqrt(R2_2[" << p << "]) = " << R2_2[p] * 0.529177249 <<
		" || sqrt(R2_3[" << p << "]) = " << pow(R2_3[p],1.0/3.0) * 0.529177249 <<
		" || sqrt(R2_M1[" << p << "]) = " << pow(R2_M1[p],-1.0) * 0.529177249 <<
		" || sqrt(R2_M2[" << p << "]) = " << pow(R2_M2[p],(-1.0/2.0)) * 0.529177249 <<
		" || sqrt(R2_M3[" << p << "]) = " << pow(R2_M3[p],(-1.0/3.0)) * 0.529177249 << endl;
		
		cout << "Theta[" << p << "] = " << setprecision(6) << Theta[p] * 180.0/3.14159265359 << endl;
	}
	cout << eivectype << " állapotfüggvényre várható érték számítás vége" << endl;	*/
	ofstream gausssigma;
	gausssigma.open ("gausssigma.txt");
	gausssigma.setf(ios_base::fixed,ios_base::floatfield);
	gausssigma << setprecision(10) << R1_2[0] - (R1[0] * R1[0]) << endl;
	gausssigma << setprecision(10) << R2_2[0] - (R2[0] * R2[0]) << endl;
	gausssigma << setprecision(10) << Theta2[0] - (Theta[0] * Theta[0]) << endl;
	gausssigma.close();
}

double RQiteration(double *val, double *val2, int *ncol, int *nrow, double *ud, int *diagind, int n, double eps, int maxiter, int q, string &eivectype)
{
	clock_t t;
	t = clock();
	cout << "Rayleigh-hányados iteráció kezdődik" << endl;
	int i, j, k, iter;
	double s, lambda = 0.0, norm = 0.0;
	double *x = new double[n];
	double *y = new double[n];
	double *r = new double[n];
//	Vec_DP x(n), y(n);
	const double ehcm = 219474.630;
	ofstream file1;
	file1.open("RQiter.txt");
	file1.setf(ios_base::fixed,ios_base::floatfield);
	ofstream file2;
	file2.open("RQiterminresiter.txt");
	file2.setf(ios_base::fixed,ios_base::floatfield);
	ifstream wf;
	wf.open ("thing.dat");
	wf.setf(ios_base::fixed,ios_base::floatfield);
	for(i = 0; i < n; i++)
	{
		wf >> x[i];
	}
	wf.close();
/*	for(i = 0; i < n; i++)
	{
		x[i] = 1.0;
	}
	norm = 0.0;
	for(i = 0; i < n; i++)
	{
		norm = norm + x[i] * x[i];
	}
	norm = sqrt(norm);
	for(i = 0; i < n; i++)
	{
		x[i] = (1.0 / norm) * x[i];
	}	*/
	mult(n,x,y,val,ncol,nrow);
	for(i = 0; i < n; i++)
	{
		lambda = lambda + x[i] * y[i];
	}
	for(k = 1; k <= maxiter; k++)
	{
		iter = k;
		transform(q,val,val2,ncol,nrow,ud,diagind,lambda);
/*		for(i = 0; i < q; i++)
		{
//			cout << "val2 = " << val2[i] << endl;
		}
		for(i = 0; i < n; i++)
		{
//			cout << "diag = " << diag[i] << endl;
		}	*/
		minres(n,maxiter,pow(10.0,-12.0),val,ncol,nrow,ud,y,x,file2);
//		CGM(n,maxiter,pow(10.0,-12.0),val,ncol,nrow,ud,diagind,y,x,file,3,file2);
//		AZ RÉGI RQ-ITERÁCIÓBAN EZ ITT NINCS
		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + y[i] * y[i];
		}
		norm = sqrt(norm);
		for(i = 0; i < n; i++)
		{
			x[i] = (1.0 / norm) * y[i];
		}
		mult(n,x,y,val,ncol,nrow);
//		************************************************
		lambda = 0.0;
		for(i = 0; i < n; i++)
		{
			lambda = lambda + x[i] * y[i];
		}
		cout << 1.0 / lambda << endl;
		for(i = 0; i < n; i++)
		{
			r[i] = y[i] - lambda * x[i];
		}
		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + r[i] * r[i];
		}
		norm = sqrt(norm);
//		konvergencia ellenőrzése:
		if(norm <= eps)
		{
			eivectype = "RQ-iteration";
			file1 << "Konvergált!\n";
			file1 << "iter = " << iter << endl;
			file1 << "Konvergált!\n";
			file1 << "lambda / hartree = " << setprecision(30) << 1.0/ lambda << " lambda / cm-1 = " << (1.0 / lambda) * ehcm << endl;
			ofstream eivec;
			eivec.open("RQEIVEC.txt");
			eivec.setf(ios_base::fixed,ios_base::floatfield);
			for(i = 0; i < n; i++)
			{
				eivec << setprecision(30) << x[i] << endl;
			}
			eivec.close();
			t =  clock() - t;
			file1 << "Rayleigh-hányados iteráció vége" << endl;
        		file1 << "Rayleigh-hányados iteráció: " << t << " clicks, " << (float) t / CLOCKS_PER_SEC << " seconds." << endl;
			file1.close();
			cout << "Rayleigh-hányados iteráció vége" << endl;
			return lambda;
		}
//		A RÉGI RQ-ITERÁCIÓBAN EZ ITT VAN
/*		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + y[i] * y[i];
		}
		norm = sqrt(norm);
		for(i = 0; i < n; i++)
		{
			x[i] = (1.0 / norm) * y[i];
		}	*/
//		************************************************
		if(k == maxiter)
		{
			eivectype = "RQ-iteration";
			file1 << "Elértem a maximális lépésszámot!\n";
			file1 << "iter = " << iter << endl;
			file1 << "Elértem a maximális lépésszámot!\n";
			file1 << "lambda / hartree = " << setprecision(30) << 1.0/ lambda << " lambda / cm-1 = " << (1.0 / lambda) * ehcm << endl;
			ofstream eivec;
			eivec.open("RQEIVEC.txt");
			eivec.setf(ios_base::fixed,ios_base::floatfield);
			for(i = 0; i < n; i++)
			{
				eivec << setprecision(30) << x[i] << endl;
			}
			eivec.close();
		}
	}
	file1 << endl;
	t =  clock() - t;
	file1 << "Rayleigh-hányados iteráció vége" << endl;
	file1 << "Rayleigh-hányados iteráció: " << t << " clicks, " << (float) t / CLOCKS_PER_SEC << " seconds." << endl;
	file1.close();
	file2.close();
	cout << "Rayleigh-hányados iteráció vége" << endl;
	return lambda;
}

double invpower(double *val, int *ncol, int *nrow, double *ud, int *diagind, int n, double eps, int maxiter, string &eivectype)
{
	clock_t t;
	t = clock();
	cout << "Inverz hatványiteráció kezdődik" << endl;
	int i, j, k, iter;
	double s, lambda = 0.0, norm = 0.0;
	double *x = new double[n];
	double *y = new double[n];
	double *r = new double[n];
//	Vec_DP x(n), y(n);
	const double ehcm = 219474.630;
	ofstream file1;
	file1.open("invhatvany.txt");
	file1.setf(ios_base::fixed,ios_base::floatfield);
	ofstream file2;
	file2.open("invhatvanycgmiter.txt");
	file2.setf(ios_base::fixed,ios_base::floatfield);
	ifstream wf;
	wf.open("thing.dat");
	wf.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < n; i++)
	{
		wf >> x[i];
	}
	wf.close();
/*	for(i = 0; i < n; i++)
	{
		x[i] = 1.0;
	}
	norm = 0.0;
	for(i = 0; i < n; i++)
	{
		norm = norm + x[i] * x[i];
	}
	norm = sqrt(norm);
	for(i = 0; i < n; i++)
	{
		x[i] = (1.0 / norm) * x[i];
	}	*/
	for(k = 1; k <= maxiter; k++)
	{
		iter = k;
//		BiCGM(n,maxiter,eps,val2,ncol,nrow,ud,y,x);
		CGM(n,10000,pow(10.0,-12.0),val,ncol,nrow,ud,diagind,y,x,file1,3,file2);
//		AZ RÉGI INVEZITERÁCIÓBAN EZ ITT NINCS
		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + y[i] * y[i];
		}
		norm = sqrt(norm);
		for(i = 0; i < n; i++)
		{
			x[i] = (1.0 / norm) * y[i];
		}
		mult(n,x,y,val,ncol,nrow);
//		************************************************
		lambda = 0.0;
		for(i = 0; i < n; i++)
		{
			lambda = lambda + x[i] * y[i];
		}
		for(i = 0; i < n; i++)
		{
			r[i] = y[i] - lambda * x[i];
		}
		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + r[i] * r[i];
		}
		norm = sqrt(norm);
//		konvergencia ellenőrzése:
		if(norm <= eps)
		{
			eivectype = "Inverse-Power-iteration";
			file1 << "Konvergált!\n";
			file1 << "iter = " << iter << endl;
			file1 << "Konvergált!\n";
			file1 << endl;
			file1 << "lambda / hartree = " << setprecision(30) << lambda << " lambda / cm-1 = " << lambda * ehcm << endl;
			ofstream eivec;
			eivec.open("INVPOWEIVEC.txt");
			eivec.setf(ios_base::fixed,ios_base::floatfield);
			for(i = 0; i < n; i++)
			{
				eivec << setprecision(30) << x[i] << endl;
			}
			eivec.close();
			file1 << endl;
			t =  clock() - t;
			file1 << "Inverz hatványiteráció vége" << endl;
			file1 << "Inverz hatványiteráció: " << t << " clicks, " << (float) t / CLOCKS_PER_SEC << " seconds." << endl;
			file1.close();
			return lambda;
		}
//		A RÉGI INVERZITERÁCIÓBAN EZ ITT VAN
/*		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + y[i] * y[i];
		}
		norm = sqrt(norm);
		for(i = 0; i < n; i++)
		{
			x[i] = (1.0 / norm) * y[i];
		}	*/
//		************************************************
		if(k == maxiter)
		{
			eivectype = "Inverse-Power-iteration";
			file1 << "Elértem a maximális lépésszámot!\n";
			file1 << "iter = " << iter << endl;
			file1 << "Elértem a maximális lépésszámot!\n";
			file1 << endl;
			file1 << "lambda / hartree = " << setprecision(30) << lambda << " lambda / cm-1 = " << lambda * ehcm << endl;
			ofstream eivec;
			eivec.open("INVPOWEIVEC.txt");
			eivec.setf(ios_base::fixed,ios_base::floatfield);
			for(i = 0; i < n; i++)
			{
				eivec << setprecision(30) << x[i] << endl;
			}
			eivec.close();
		}
	}
	file1 << endl;
	t =  clock() - t;
	file1 << "Inverz hatványiteráció vége" << endl;
	file1 << "Inverz hatványiteráció: " << t << " clicks, " << (float) t / CLOCKS_PER_SEC << " seconds." << endl;
	file1.close();
	file2.close();
	cout << "Inverz hatványiteráció vége" << endl;
	return lambda;
}

double power(double *val, int *ncol, int *nrow, int *diagind, int n, double eps, int maxiter, int nonzero, string &eivectype)
{
	clock_t t;
	t = clock();
	cout << "Hatványiteráció kezdődik" << endl;
	int i, j, k, iter, shift = 1;
	double s, lambda = 0.0, store, norm = 0.0;
	double *x = new double[n];
	double *y = new double[n];
	double *r = new double[n];
	double *val2 = new double[nonzero];
//	Vec_DP x(n), y(n);
	const double ehcm = 219474.630;
	ofstream file;
	file.open("hatvany.txt");
	file.setf(ios_base::fixed,ios_base::floatfield);
	if(shift == 1)
	{
		for(int i = 0; i < nonzero; i++)
		{
			val2[i] = val[i];
		}
	}
A:	if(shift == 0)
	{
		cout << "Shift kezdődik" << endl;
		for(int i = 0; i < n; i++)
		{
			val2[diagind[i]] = val[diagind[i]] - lambda;
		}
	}
	ifstream wf;
	wf.open ("thing.dat");
	wf.setf(ios_base::fixed,ios_base::floatfield);
	for(i = 0; i < n; i++)
	{
		wf >> x[i];
	}
	wf.close();
/*	for(i = 0; i < n; i++)
	{
		x[i] = 1.0;
	}
	norm = 0.0;
	for(i = 0; i < n; i++)
	{
		norm = norm + x[i] * x[i];
	}
	norm = sqrt(norm);
	for(i = 0; i < n; i++)
	{
		x[i] = (1.0 / norm) * x[i];
	}	*/
	for(k = 1; k <= maxiter; k++)
	{
		iter = k;
		mult(n,x,y,val2,ncol,nrow);
//		AZ RÉGI HATVÁNYITERÁCIÓBAN EZ ITT NINCS
		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + y[i] * y[i];
		}
		norm = sqrt(norm);
		for(i = 0; i < n; i++)
		{
			x[i] = (1.0 / norm) * y[i];
		}
		mult(n,x,y,val2,ncol,nrow);
//		************************************************
		lambda = 0.0;
		for(i = 0; i < n; i++)
		{
			lambda = lambda + x[i] * y[i];
		}
		for(i = 0; i < n; i++)
		{
			r[i] = y[i] - lambda * x[i];
		}
		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + r[i] * r[i];
		}
		norm = sqrt(norm);
//		konvergencia ellenőrzése:
		if(norm <= eps)
		{
			eivectype = "Power-iteration";
			file << "Konvergált!\n";
			file << "iter = " << iter << endl;
			file << "Konvergált!\n";
			file << endl;
			if(shift == 1)
			{
				file << "lambda / hartree = " << setprecision(30) << lambda << " lambda / cm-1 = " << lambda * ehcm << endl;
				store = lambda;
			}
			if(shift == 0)
			{
				file << "lambda / hartree = " << setprecision(30) << store + lambda << " lambda / cm-1 = " << (store + lambda) * ehcm << endl;
				ofstream eivec;
				eivec.open("POWEIVEC.txt");
				eivec.setf(ios_base::fixed,ios_base::floatfield);
				for(i = 0; i < n; i++)
				{
					eivec << setprecision(30) << x[i] << endl;
				}
				eivec.close();
			}
			file << endl;
			t =  clock() - t;
			file << "Hatványiteráció vége" << endl;
        	file << "Hatványiteráció: " << t << " clicks, " << (float) t / CLOCKS_PER_SEC << " seconds." << endl;
			if(shift == 1)
			{
				shift = 0;
				goto A;
			}
			file.close();
			cout << "Hatványiteráció vége" << endl;
			return lambda;
		}
//		A RÉGI HATVÁNYITERÁCIÓBAN EZ ITT VAN
/*		norm = 0.0;
		for(i = 0; i < n; i++)
		{
			norm = norm + y[i] * y[i];
		}
		norm = sqrt(norm);
		for(i = 0; i < n; i++)
		{
			x[i] = (1.0 / norm) * y[i];
		}	*/
//		************************************************
		if(k == maxiter)
		{
			eivectype = "Power-iteration";
			file << "Elértem a maximális lépésszámot!\n";
			file << "iter = " << iter << endl;
			file << "Elértem a maximális lépésszámot!\n";
			file << endl;
			if(shift == 1)
			{
				file << "lambda / hartree = " << setprecision(30) << lambda << " lambda / cm-1 = " << lambda * ehcm << endl;
			}
			if(shift == 0)
			{
				file << "lambda / hartree = " << setprecision(30) << store + lambda << " lambda / cm-1 = " << (store + lambda) * ehcm << endl;
				ofstream eivec;
				eivec.open("POWEIVEC.txt");
				eivec.setf(ios_base::fixed,ios_base::floatfield);
				for(i = 0; i < n; i++)
				{
					eivec << setprecision(30) << x[i] << endl;
				}
				eivec.close();
			}
		}
	}
	file << endl;
	t =  clock() - t;
	file << "Hatványiteráció vége" << endl;
	file << "Hatványiteráció: " << t << " clicks, " << (float) t / CLOCKS_PER_SEC << " seconds." << endl;
	file.close();
	cout << "Hatványiteráció vége" << endl;
	return lambda;
}

void powerfordavidson(double **A, double *val, int *ncol, int *nrow, double *y, double &lambda2, int n, double eps, int maxiter)
{
	int i, j, k, iter, conv;
	double s, lambda1 = 0.0;
	lambda2 = 0.0;
	double *x = new double[n];
//	Vec_DP x(n);
	for(i = 0; i < n; i++)
	{
		x[i] = 1.0;
	}
	for(k = 1; k <= maxiter; k++)
	{
		iter = k;
		for(i = 0; i < n; i++)
		{
			y[i] = 0.0;
		}
//		minresfordavidson(n,maxiter,eps,val,ncol,nrow,A,y,x,0);
//		CGMfordavidson(n,maxiter,eps,A,y,x);
/*		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n; j++)
			{
				y[i] = y[i] + A[i][j] * x[j];
			}
		}	*/
		s = 0.0;
		for(i = 0; i < n; i++)
		{
			s = s + y[i] * y[i];
		}
		s = sqrt(s);
		for(i = 0; i < n; i++)
		{
			y[i] = (1.0 / s) * y[i];
		}
		conv = 0;
		for(i = 0; i < n; i++)
		{
			if(abs(y[i] - x[i]) > eps)
				conv++;
		}
		for(i = 0; i < n; i++)
		{
			x[i] = 0.0;
			for(j = 0; j < n; j++)
			{
				x[i] = x[i] + A[i][j] * y[j];
			}
		}
		lambda2 = 0.0;
		for(i = 0; i < n; i++)
		{
			lambda2 = lambda2 + y[i] * x[i];
		}
//		konvergencia ellenőrzése:
		if((abs(lambda2 - lambda1) <= eps) && (conv == 0))
		{
//			cout << "powerfordavidson konvergált" << endl;
			return;
		}
		for(i = 0; i < n; i++)
		{
			x[i] = y[i];
		}
		lambda1 = lambda2;
		if(k == maxiter)
		{
//			cout << "Nem konvergált! (powerfordavidson) :( !\n";
		}
	}
	return;
}

// eredeti fgv-k

void SVLprecgm(int n, int nmax, int maxiter, double hiba, double *de, int ndont,int *nrow, int *ncol, double *val, double *diag)
{
	clock_t t;
	t = clock();
	ofstream out;
	out.open("res.txt");
	out.setf(ios_base::fixed,ios_base::floatfield);

	double **q = new double *[n];
	for(int i = 0; i < n; i++)
		q[i] = new double[maxiter];
	double **tv = new double *[maxiter];
	for(int i = 0; i < maxiter; i++)
		tv[i] = new double[maxiter];
	double *beta = new double[maxiter];
	double *elozo = new double[maxiter];
	double *alpha = new double[maxiter];
	double *z = new double[n];
	double *ee = new double[maxiter];
	double *decm = new double[maxiter];
	double *x = new double[n];
	double *diff = new double[maxiter];
	double *er = new double[maxiter];
/*	Mat_DP q(n,maxiter), tv(maxiter,maxiter);
	Vec_DP beta(maxiter),elozo(maxiter), alpha(maxiter);
	Vec_DP z(n),ee(maxiter), decm(maxiter);
	Vec_DP x(n),diff(maxiter),er(maxiter);	*/
	double beta1,vorma,gamma1;

	int neig = nmax;
	int nstep = 1;
	int nkezd = 1;
	
	cout << "kezdődik a régi Lánczos" << endl;
         
	for(int i = 0; i < maxiter; i++)
	{
		alpha[i] = 0.00;
		beta[i] = 0.00;
		elozo[i] = 0.00;
	}

	for(int i = 0; i < n; i++)
	{
		q[i][0] = diag[i];
	}

	for(int i = 0; i < n; i++)
	{
		z[i] = q[i][0];
	}

	precond(n,x,z, nrow, ncol, val, diag);

	for(int i = 0; i < n; i++)
	{
		z[i] = x[i];
		//cout<<z[i]<<endl;
	}

	vorma = 0.00;

	norman(q, vorma, n, 0);

	if(vorma == 0.00)
	{
		cout<<"vorma = 0.00"<<endl;
		exit(-1);
	}

	for(int i = 0; i < n; i++)
	{
		q[i][0] = q[i][0]/vorma;
	}

// ********** ELSO OSZLOP *************

	for(int i = 0; i < n; i++)
	{
		z[i] = q[i][0];
	}

	precond(n,x,z, nrow, ncol, val, diag);

	for(int i = 0; i < n; i++)
	{
		z[i] = x[i];
	}

	for(int i = 0; i < n; i++)
	{
		alpha[0] = alpha[0] + z[i]*q[i][0];
	}

	for(int i = 0; i < n; i++)
	{
		z[i] =z[i] - alpha[0]*q[i][0];
	}

	norm(z,n, beta1);

	if(beta1 == 0.00)
	{
		cout<<"beta1 = 0.00"<<endl;
		exit(-1);
	}

	beta[0] = beta1;


	for(int i = 0; i < n; i++)
	{
		q[i][1] = z[i] / beta1;
	}

	out<<"***************"<<nstep<<". CIKLUS **************"<<endl;
	out<<"nmax = "<<nmax<<endl;
	out<<"nkezd"<<nkezd<<endl;

	if(nmax > maxiter)
	{
		out<<"FIGYELEM! NINCS TELJES KONVERGENCIA!"<<endl;
		exit(-1);
	}

//
//	FOCIKLUS
//  eloszor nkezd = 1
//

	bool doit = true;

	while(doit)
	{
		for(int j = nkezd; j < nmax; j++)
		{
			for(int i = 0; i < n; i++)
			{
				z[i] = q[i][j];
			}

			precond(n,x,z, nrow, ncol, val, diag);

			for(int i = 0; i < n; i++)
			{
				z[i] = x[i];
			}

			alpha[j] = 0.00;
			beta[j] = 0.00;

			
			for(int i = 0; i < n; i++)
			{
				alpha[j] = alpha[j] + z[i]*q[i][j];
			}

			for(int i = 0; i < n; i++)
			{
				z[i] =z[i] - alpha[j]*q[i][j] - beta[j-1]*q[i][j-1];
			}
	  
			// ortogonalizacio

			for(int i = 0; i < j; i++)
			{
				gamma1 = 0.00;

				for(int k = 0; k < n; k++)
				{
					gamma1 = gamma1 + q[k][j-i+1]*z[k];
				}

				for(int k = 0; k < n; k++)
				{
					z[k] = z[k] - gamma1*q[k][j-i+1];
				}
			}

			beta1 = 0.00;

			norm(z, n, beta1);

			if(beta1 == 0.00)
			{
				cout<<"beta1 = 0.00"<<endl;
				exit(-1);
			}

			beta[j] = beta1;

			for(int i = 0; i < n; i++)
			{
				q[i][j+1] = z[i]/beta[j];
			}
		} // fociklus

		//cout<<"nmax = "<<nmax<<endl;;
		for(int i = 0; i< nmax; i++)
		{
			de[i] = alpha[i];
			ee[i] = beta[i];
			//cout<<alpha[i]<<"  "<<beta[i]<<endl;
		}

		for(int i = 0; i< nmax; i++)
		{
			for(int k = 0; k < nmax; k++)
			{
				if(i == k)
				{
					tv[i][k]=1.00;
				}
				else
				{
					tv[i][k] = 0.00;
				}
			}
		}

		tqli2(de, ee, tv, nmax);
		eigsrt(de,tv,nmax);

		//cout<<"sajaterteke:"<<endl;
		for(int j = 0; j < nmax; j++)
		{
			
			decm[j] = 0.00;
			diff[j] = 0.00;

			de[j] = 1.00/ de[j];
			decm[j] = de[j]*219474.63;
			//cout<<decm[j]<<endl;
			diff[j] = abs(elozo[j]-decm[j]);
			elozo[j] = decm[j];
		}

		out<<"*****"<<endl;
		out<<"Az "<<nstep<<". CIKLUS SAJATERTEKEI: "<<endl;

		nstep++;

		for(int i = 0; i < neig; i++)
		{
			out<<i<<"  "<<decm[i]<<"  "<<diff[i]<<endl;
		}

		int nconv =0;

		for(int i = 0; i < neig; i++)
		{
			if(diff[i] > hiba)
				nconv++;
		}

		if(nconv > 0)
		{
			out<<""<<endl;
			out<<"NEM KONVERGALT ADATOK! UJ CIKLUS KEZDODIK!"<<endl;
			nkezd = nmax;
			nmax = nmax+1;
			
		} else
		{

			out<<"A Lanczos konvergalt!"<<endl;
			doit = false;
		}


	}// while



	t = clock() - t;
	out << "régi Lánczos: " << t << " clicks, " << (float) t / CLOCKS_PER_SEC << " seconds." << endl;

}

void norm(double *v, int n, double &orma)
{
	orma = 0.00;
	double vxv = 0.00;

	for(int i = 0; i < n; i++)
	{
		vxv = vxv + v[i]*v[i];
	}

	orma = sqrt(vxv);
}

void norman(double **r, double &beta, int n, int j)
{
	beta = 0.00;
	double sum = 0.00;

	for(int i = 0; i < n; i++)
	{
		sum = sum + r[i][j]*r[i][j];
	}

	beta = sqrt(sum);
}

void precond(int n, double *x, double *b ,int *nrow, int *ncol, double *val, double *diag)
{
	double *r = new double[n];
	double *z = new double[n];
	double *v = new double[n];
//	Vec_DP r(n),z(n),v(n);
	double delta, c,t,d;

	delta = 1E-9;
	int m = 150;


	for(int i =0; i < n; i++)
	{
		r[i] = 0.00;
		for(int j = nrow[i]; j < nrow[i+1]; j++)
		{
			r[i] = r[i] + val[j]*x[ncol[j]];
		}
	}

	for(int i =0; i < n; i++)
	{
		r[i] = b[i]-r[i];
	}

	for(int i =0; i < n; i++)
	{
		z[i] = r[i]/diag[i];
		v[i] = z[i];
	}

	c=0.00;

	for(int i =0; i < n; i++)
	{
		c = c+z[i]*r[i];
	}

	double vv, vz;
	for(int k = 0; k < m; k++)
	{
		vv = 0.00;
		for(int i =0; i < n; i++)
		{
			vv = vv + v[i]*v[i];
		}

		if(sqrt(vv) < delta)
			break;

		for(int i =0; i < n; i++)
		{
			z[i] = 0.00;
			for(int j = nrow[i]; j < nrow[i+1]; j++)
			{
				z[i] = z[i] + val[j]*v[ncol[j]];
			}
		}

		vz = 0.00;

		for(int i =0; i < n; i++)
		{
			vz = vz + z[i]*v[i];
		}

		t=c/vz;

		for(int i =0; i < n; i++)
		{
			x[i] = x[i] + t*v[i];
		}

		for(int i =0; i < n; i++)
		{
			r[i] = r[i] - t*z[i];
 			z[i] = r[i]/diag[i];
		}

		d=0.00;

		for(int i =0; i < n; i++)
		{
			d = d+ z[i]*r[i];
		}

		for(int i =0; i < n; i++)
		{
			v[i] = z[i] + (d/c)*v[i];
		}

		c=d;
	}
}

void potv(double &vpot, double x1, double x2, double x3)
{
	vpot = 0.00;
}

void legmat2(int np, int km, int km1, double ajm, double **h4, double **ht3, int np_cmat)
{
	double *q = new double[np];
	double **h3 = new double *[np];
	for(int i = 0; i < np; i++)
		h3[i] = new double[np];
	double **tv = new double *[np];
	for(int i = 0; i < np; i++)
		tv[i] = new double[np];
	double **e3 = new double *[np];
	for(int i = 0; i < np; i++)
		e3[i] = new double[np];
	double **tv1 = new double *[np];
	for(int i = 0; i < np; i++)
		tv1[i] = new double[np];
	double **tv3 = new double *[np];
	for(int i = 0; i < np; i++)
		tv3[i] = new double[np];
/*	Vec_DP q(np);
	Mat_DP h3(np,np),tv(np,np) ;
	Mat_DP e3(np,np),tv1(np,np),tv3(np,np);	*/ 

	double akm= (double)km;

	for(int i = 0; i < np; i++)
	{
		for(int j = 0; j < np; j++)
		{
			h3[i][j] = 0.00;
		}
	}

	double ai;
	for(int i = km; i <= np+km-2; i++)
	{
		ai = (double)i;
		h3[i-km][i-km+1] = sqrt((ai+1.00-akm)*(ai+1.00+akm)/(4.00*(ai+1.00)*(ai+1.00)-1.00));
	}

	for(int i = 0; i < np; i++ ) {
		for(int j = 0; j < np; j++ ) {
			h3[j][i] = h3[i][j];
		}
	}

	

	int nrot = 100;
	jacobi(np, h3, q, tv, nrot);
	eigsrt2(np, q, tv);

	double akm1= (double)km1;

	for(int i = 0; i < np; i++)
	{
		for(int j = 0; j < np; j++)
		{
			h3[i][j] = 0.00;
		}
	}

	for(int i = km1; i <= np+km1-2; i++)
	{
		ai = (double)i;
		h3[i-km1][i-km1+1] = sqrt((ai+1.00-akm1)*(ai+1.00+akm1)/(4.00*(ai+1.00)*(ai+1.00)-1.00));
	}

	for(int i = 0; i < np; i++ ) {
		q[i] = 0.00;
		for(int j = 0; j < np; j++ ) {
			h3[j][i] = h3[i][j];
		}
	}
	
	nrot = 100;
	jacobi(np, h3, q, tv1, nrot);
	eigsrt2(np, q, tv1);

	for(int i = 0; i < np; i++)
	{
		for(int j = 0; j < np; j++)
		{
			e3[i][j] = 0.00;
		}
	}

	double coef=-1.000;
    
	if(km == 0)
		coef = -1.00*sqrt(2.00);

	for(int i = 1; i < np; i++)
	{
		ai=(double)(i+1+km);
	    e3[i][i-1]= coef*sqrt((ai-1.00)*ai-akm*(akm+1.00));
	}

	//cout<<"km = "<<km<<endl;
	//e3.Print();

	for(int i = 0; i < np; i++)
	{
		for(int j = 0; j < np; j++)
		{
			h3[i][j] = 0.00;
			h4[i][j] = 0.00;
		}
	}

	for(int i = 0; i < np; i++)
	{
		for(int j = 0; j < np; j++)
		{
			tv3[i][j] = tv[j][i];
		}
	}

	MMUL(tv3,e3,h3,np,np,np);

	MMUL(h3,tv1,h4,np,np,np);

	for(int i = 0; i < np; i++)
	{
		for(int j = 0; j < np; j++)
		{
			h4[i][j]=h4[i][j]*sqrt(ajm*(ajm+1.00)-akm*(akm+1.00));
		}
	}

	for(int i = 0; i < np; i++)
	{
		for(int j = 0; j < np; j++)
		{
			ht3[i][j] = h4[j][i];
		}
	}

	//cout<<"km = "<<km<<endl;
	//ht3.Print();

	
}

void MMUL(double **A, double **B, double **C, int m, int n, int nm)
{
	double temp;

	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < nm; j++)
		{
			temp = 0.00;
			for(int k = 0; k < n; k++)
			{
				 temp = temp  +A[i][k]*B[k][j];
			}

			C[i][j] = temp;
		}
	}
}

void legmat(int np, int km, double **h3, double *q, int np_lmat, ofstream &grid)
{
	double **tv = new double *[np];
	for(int i = 0; i < np; i++)
		tv[i] = new double[np];
//	Mat_DP tv(np,np);
	double akm = (double)km;

	for(int i = 0; i < np; i++ ) {
		for(int j = 0; j < np; j++ ) {
			h3[i][j] = 0.00;
		}
	}
	
	double ai;
	for(int i = km; i <= np+km-2; i++)
	{
		ai = (double)i;
		h3[i-km][i-km+1] = sqrt((ai+1.00-akm)*(ai+1.00+akm)/(4.00*(ai+1.00)*(ai+1.00)-1.00));
	}

	for(int i = 0; i < np; i++ ) {
		for(int j = 0; j < np; j++ ) {
			h3[j][i] = h3[i][j];
		}
	}

	//h3.Print();

	int nrot = 100;
	jacobi(np, h3, q, tv, nrot);
	eigsrt2(np, q, tv);

	//for(int i = 0; i < np; i++)
	//	cout<<q[i]<<" ";

	//cout<<endl;
	
	for(int i = 0; i < np; i++)
	{
		grid << setprecision(30) << q[i] << endl;
	}
	grid.close();

//******Az L2 matrix transzformalasa DVR-be*******************	

	for(int k = 0; k < np; k++)
	{
		for(int l = 0; l < np; l++)
		{
			h3[k][l] = 0.0;
			for(int i = 0; i < np; i++)
			{
				ai = (double)(i+km);
//				cout << "ai = " << ai << endl;
				h3[k][l] = h3[k][l] + tv[i][k]*((ai+1.0)*ai)*tv[i][l];
			}
		}
	}

/*	for(int i = 0; i < np; i++)
	{
		for(int j = 0; j < np; j++)
		{
			cout << setprecision(10) << h3[i][j] << " ";
		}
		cout << endl;
	}	*/
	//h3.Print();

}


void jacbolbelso2(double tm3,double tm2,double R1j,double R2j,double xcosj,double &R1b,double &R2b,double &xcosb)
{
	R1b=R1j;
	R2b = sqrt((tm2/(tm3+tm2))*(tm2/(tm3+tm2))*R1j*R1j+R2j*R2j+2*tm2/(tm3+tm2)*R1j*R2j*xcosj);
	xcosb = R2j/R2b*xcosj+tm2/(tm3+tm2)*R1j/R2b;
}

/*
double norm(double *v)
{
	int n = v.size();
	double vxv = 0.0;
	for(int i = 0; i < n; i++)
		vxv=vxv+v[i]*v[i];

	return sqrt(vxv);
}
*/

inline void rot(double **a, const double s, const double tau, const int i, const int j, const int k, const int l)
{
	double g,h;
	g=a[i][j];
	h=a[k][l];
	a[i][j]=g-s*(h+g*tau);
	a[k][l]=h+s*(g-h*tau);
}

void jacobi(int n, double **a, double *d, double **v, int &nrot)
{
	int i,j,ip,iq;
	double tresh,theta,tau,t,sm,s,h,g,c;
//	DP tresh,theta,tau,t,sm,s,h,g,c;

//	int n=d.size();
	double *b = new double[n];
	double *z = new double[n];
//	Vec_DP b(n),z(n);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0)
			return;
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
						a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1.0+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip;j++)
						rot(a,s,tau,j,ip,j,iq);
					for (j=ip+1;j<iq;j++)
						rot(a,s,tau,ip,j,j,iq);
					for (j=iq+1;j<n;j++)
						rot(a,s,tau,ip,j,iq,j);
					for (j=0;j<n;j++)
						rot(v,s,tau,j,ip,j,iq);
					++nrot;
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	cout<<"Too many iterations in routine jacobi";
}


void jacobifordavidson(int n, double **a, double *d, double **v, int &nrot)
{
	int i,j,ip,iq;
	double tresh,theta,tau,t,sm,s,h,g,c;
//	DP tresh,theta,tau,t,sm,s,h,g,c;

//	int n=d.size();
	double *b = new double[n];
	double *z = new double[n];
//	Vec_DP b(n),z(n);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0)
			return;
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
						a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1.0+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip;j++)
						rot(a,s,tau,j,ip,j,iq);
					for (j=ip+1;j<iq;j++)
						rot(a,s,tau,ip,j,j,iq);
					for (j=iq+1;j<n;j++)
						rot(a,s,tau,ip,j,iq,j);
					for (j=0;j<n;j++)
						rot(v,s,tau,j,ip,j,iq);
					++nrot;
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	cout<<"Too many iterations in routine jacobi";
}

/*void eigsrt(double *d, double **v, int n)
{
	int i,j,k;
	double p;

	//int n=d.size();
	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}	*/

void eigsrt2(int n, double *d, double **v)
{
	int i,j,k;
	double p;

//	int n=d.size();
	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i;j<n;j++)
			if (d[j] <= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

void tred2(int n, double **a, double *d, double *e)
{
	int l,k,j,i;
	double scale,hh,h,g,f;
//	DP scale,hh,h,g,f;

//	int n=d.size();
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<l+1;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<l+1;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<l+1;j++) {
				// Next statement can be omitted if eigenvectors not wanted
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<l+1;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<l+1;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	// Next statement can be omitted if eigenvectors not wanted
	//d[0]=0.0;
	//e[0]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	//	wanted except for statement d[i]=a[i][i];
	for (i=0;i<n;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (j=0;j<l;j++) {
				g=0.0;
				for (k=0;k<l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
	}
}

void tred2fordavidson(int n, double **a, double *d, double *e)
{
	int l,k,j,i;
	double scale,hh,h,g,f;
//	DP scale,hh,h,g,f;

//	int n=d.size();
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<l+1;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<l+1;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<l+1;j++) {
				// Next statement can be omitted if eigenvectors not wanted
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<l+1;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<l+1;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	// Next statement can be omitted if eigenvectors not wanted
	//d[0]=0.0;
	//e[0]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	//	wanted except for statement d[i]=a[i][i];
	for (i=0;i<n;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (j=0;j<l;j++) {
				g=0.0;
				for (k=0;k<l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
	}
}

double pythag(const double a, const double b)
{
	double absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void tqli2(double *d, double *e, double **z, int n)
{
   register int   m,l,iter,i,k;
   double         s,r,p,g,f,dd,c,b;

   //for(i = 1; i < n; i++) e[i-1] = e[i];
     e[n] = 0.0;
   for(l = 0; l < n; l++) {
      iter = 0;
      do {
         for(m = l; m < n-1; m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
            if((double)(fabs(e[m])+dd) == dd) break;
         }
         if(m != l) {
            if(iter++ == 30) {
               printf("\n\nToo many iterations in tqli.\n");
               exit(1);
            }
            g = (d[l+1] - d[l])/(2.0 * e[l]);
            r = pythag(g,1.0);
            g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
            s = c = 1.0;
            p = 0.0;
            for(i = m-1; i >= l; i--) {
               f      = s * e[i];
               b      = c*e[i];
               e[i+1] = (r=pythag(f,g));
               if(r == 0.0) {
                  d[i+1] -= p;
                  e[m]    = 0.0;
                  break;
               }
               s      = f/r;
               c      = g/r;
               g      = d[i+1] - p;
               r      = (d[i] - g) * s + 2.0 * c * b;
               d[i+1] = g + (p = s * r);
               g      = c * r - b;
               for(k = 0; k < n; k++) {
                  f         = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i]   = c * z[k][i] - s * f;
               } /* end k-loop */
            } /* end i-loop */
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
         } /* end if-loop for m != 1 */
      } while(m != l);
   } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */

void tqli2fordavidson(double *d, double *e, double **z, int n)
{
   register int   m,l,iter,i,k;
   double         s,r,p,g,f,dd,c,b;

   //for(i = 1; i < n; i++) e[i-1] = e[i];
     e[n] = 0.0;
   for(l = 0; l < n; l++) {
      iter = 0;
      do {
         for(m = l; m < n-1; m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
            if((double)(fabs(e[m])+dd) == dd) break;
         }
         if(m != l) {
            if(iter++ == 30) {
               printf("\n\nToo many iterations in tqli.\n");
               exit(1);
            }
            g = (d[l+1] - d[l])/(2.0 * e[l]);
            r = pythag(g,1.0);
            g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
            s = c = 1.0;
            p = 0.0;
            for(i = m-1; i >= l; i--) {
               f      = s * e[i];
               b      = c*e[i];
               e[i+1] = (r=pythag(f,g));
               if(r == 0.0) {
                  d[i+1] -= p;
                  e[m]    = 0.0;
                  break;
               }
               s      = f/r;
               c      = g/r;
               g      = d[i+1] - p;
               r      = (d[i] - g) * s + 2.0 * c * b;
               d[i+1] = g + (p = s * r);
               g      = c * r - b;
               for(k = 0; k < n; k++) {
                  f         = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i]   = c * z[k][i] - s * f;
               } /* end k-loop */
            } /* end i-loop */
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
         } /* end if-loop for m != 1 */
      } while(m != l);
   } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */

void tqli(double *d, double *e, double **z, int n)
{
	int m,l,iter,i, k;
	double s,r,p,g,f,dd,c,b;
//	DP s,r,p,g,f,dd,c,b;

	//int n=d.size();
	//for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) cout<<"Too many iterations in tqli";
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					// Next loop can be omitted if eigenvectors not wanted
					for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

void psikiind(int nt, int nsp1nt, int ns1, int np, double Re, double Rmax, double tm1atom, double tm2atom, double tm3atom, double *val, double *diag)
{
	double x, y, z, norm = 0.0, scal = 1.0;
	double mu_q = 0.0, mu_de = 0.0;
	double sum_mu_r1 = 1.84371823805856660527524581994, sum_mu_r2 = 1.8437182375102556530777064836, sum_mu_theta = cos(1.822640511752868563988272399);
	double sum_mu_r1_sq = 3.41627648866375688285756950791, sum_mu_r2_sq = 3.41627648643212822321402200032, sum_mu_theta_sq = cos(3.34606286253576046263447096862);
	double sigma_r1 = sqrt(sum_mu_r1_sq - sum_mu_r1 * sum_mu_r1), sigma_r2 = sqrt(sum_mu_r2_sq - sum_mu_r2 * sum_mu_r2);
	double sigma_theta = sqrt(sum_mu_theta_sq - sum_mu_theta * sum_mu_theta);
	
	double *f_x = new double[nt];
	double *f_y = new double[nt];
	double *f_z = new double[nt];
	double *g = new double[nt];
	double *stdev = new double[3];
	double *allgrid = new double[ns1+np];
	double *de = new double[ns1];
	double *q = new double[np];
	double *gscal1 = new double[nt];
	double *gscal2 = new double[nt];
	double *gscal3 = new double[nt];
	double *var1 = new double[nt];
	double *var2 = new double[nt];
	double *var3 = new double[nt];
	double *pot = new double[nt];
	double *mu_r1 = new double[nt];
	double *mu_r2 = new double[nt];
	double *mu_theta = new double[nt];
	
	double **hr1 = new double *[ns1];
	for(int i = 0; i < ns1; i++)
		hr1[i] = new double[ns1];
	double **tv2 = new double *[ns1];
	for(int i = 0; i < ns1; i++)
		tv2[i] = new double[ns1];
	double **h3 = new double *[np];
	for(int i = 0; i < np; i++)
		h3[i] = new double[np];

//	potenciál olvasása fájlból
/*	ifstream potencial;
	potencial.open ("POTENCIAL.txt");
	potencial.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < nt; i++)
	{
		potencial >> pot[i];
	}
	potencial.close();	*/
//	grid pontok olvasása fájlból: Hermite polinomok gyökei ns1-ig utána a Legendre polinomok gyökei
	ifstream grid;
	grid.open ("gridfile.txt");
	grid.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < ns1+np; i++)
	{
		grid >> allgrid[i];
	}
	grid.close();
	ifstream gausssigma;
	gausssigma.open ("gausssigma.txt");
	gausssigma.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < 3; i++)
	{
		gausssigma >> stdev[i];
	}
	gausssigma.close();
// itt az lett volna a terv, hogy kiveszek bizonyos grid pontokat és helyette mást teszek be, de ez hülyeség
/*	if(ns1 % 2 == 0)	//páros eset
	{
		ns1 = ns1 + 1;
		for(int i = 0; i < ns1; i++ ) {
			for(int j = 0; j < ns1; j++ ) {
				hr1[i][j] =0.00;
			}
		}
		for(int i = 0; i < ns1-1; i++ ) {
			hr1[i][i+1] = sqrt(((double)i+1.0)/2.0);
		}
		for(int i = 1; i < ns1; i++ ) {
			hr1[i][i-1] = sqrt(((double)i-1.0)/2.0);
			//	valamiért anno erre atirtam: hr1[i][i-1] = sqrt(((double)i)/2.0);
		}
		int nrot = 100;
		jacobi(ns1, hr1, de, tv2, nrot);
		eigsrt(de, tv2, ns1);
	}
	if(ns1 % 2 != 0)	//páratlan eset
	{
		for(int i = 0; i < ns1; i++ ) {
			for(int j = 0; j < ns1; j++ ) {
				hr1[i][j] =0.00;
			}
		}
		for(int i = 0; i < ns1-1; i++ ) {
			hr1[i][i+1] = sqrt(((double)i+1.0)/2.0);
		}
		for(int i = 1; i < ns1; i++ ) {
			hr1[i][i-1] = sqrt(((double)i-1.0)/2.0);
			//	valamiért anno erre atirtam: hr1[i][i-1] = sqrt(((double)i)/2.0);
		}
		int nrot = 100;
		jacobi(ns1, hr1, de, tv2, nrot);
		eigsrt(de, tv2, ns1);
		for(int i = 0; i < ns1/2; i++)
			de[i] = de[i+1];
		de[ns1/2] = 0.00;
		for(int i = ns1/2+1; i < ns1; i++)
			de[i] = de[i+1];
	}
	if(np % 2 == 0)	//páros eset
	{
		np = np + 1;
		int km = 0;
		double akm = (double)km;
		for(int i = 0; i < np; i++ ) {
			for(int j = 0; j < np; j++ ) {
				h3[i][j] = 0.00;
			}
		}
		double ai;
		for(int i = km; i <= np+km-2; i++)
		{
			ai = (double)i;
			h3[i-km][i-km+1] = sqrt((ai+1.00-akm)*(ai+1.00+akm)/(4.00*(ai+1.00)*(ai+1.00)-1.00));
		}

		for(int i = 0; i < np; i++ ) {
			for(int j = 0; j < np; j++ ) {
				h3[j][i] = h3[i][j];
			}
		}
		int nrot = 100;
		jacobi(np, h3, q, tv2, nrot);
		eigsrt2(np, q, tv2);
	}
	if(np % 2 != 0)	//páratlan eset
	{
		int km = 0;
		double akm = (double)km;
		for(int i = 0; i < np; i++ ) {
			for(int j = 0; j < np; j++ ) {
				h3[i][j] = 0.00;
			}
		}
		double ai;
		for(int i = km; i <= np+km-2; i++)
		{
			ai = (double)i;
			h3[i-km][i-km+1] = sqrt((ai+1.00-akm)*(ai+1.00+akm)/(4.00*(ai+1.00)*(ai+1.00)-1.00));
		}

		for(int i = 0; i < np; i++ ) {
			for(int j = 0; j < np; j++ ) {
				h3[j][i] = h3[i][j];
			}
		}
		int nrot = 100;
		jacobi(np, h3, q, tv2, nrot);
		eigsrt2(np, q, tv2);
	}	*/
/*	else
	{
		for(int i = 0; i < ns1; i++ ) {
			for(int j = 0; j < ns1; j++ ) {
				hr1[i][j] =0.00;
			}
		}
		for(int i = 0; i < ns1-1; i++ ) {
			hr1[i][i+1] = sqrt(((double)i+1.0)/2.0);
		}
		for(int i = 1; i < ns1; i++ ) {
			hr1[i][i-1] = sqrt(((double)i-1.0)/2.0);
			//	valamiért anno erre atirtam: hr1[i][i-1] = sqrt(((double)i)/2.0);
		}
		int nrot = 100;
		jacobi(ns1, hr1, de, tv2, nrot);
		eigsrt(de, tv2, ns1);
	
	
		int km = 0;
		double akm = (double)km;
		for(int i = 0; i < np; i++ ) {
			for(int j = 0; j < np; j++ ) {
				h3[i][j] = 0.00;
			}
		}
		double ai;
		for(int i = km; i <= np+km-2; i++)
		{
			ai = (double)i;
			h3[i-km][i-km+1] = sqrt((ai+1.00-akm)*(ai+1.00+akm)/(4.00*(ai+1.00)*(ai+1.00)-1.00));
		}

		for(int i = 0; i < np; i++ ) {
			for(int j = 0; j < np; j++ ) {
				h3[j][i] = h3[i][j];
			}
		}
		jacobi(np, h3, q, tv2, nrot);
		eigsrt2(np, q, tv2);
	}	*/
//	ahogy fent írtam, Hermite polinomok gyökei ns1-ig utána a Legendre polinomok gyökei
	for(int i = 0; i < ns1; i++)
	{
		de[i] = allgrid[i];
	}
	for(int i = 0; i < np; i++)
	{
		q[i] = allgrid[i+ns1];
	}
	if(np % 2 == 0)	//páros eset
	{
		for(int i = 0; i < np / 2 - 1; i++)
			q[i] = q[i+1];
		for(int i = np - 1; i > np / 2; i--)
			q[i] = q[i-1];
		q[np/2-1] = -0.00000;
		q[np/2] = +0.00000;
//		for(int i = 0; i < np; i++)
//			cout << setprecision(30) << q[i] << endl;
	}
/*	if(ns1 % 2 == 0)	//páros eset
	{
		de[ns1/2] = 0.00;
	}
	if(np % 2 == 0)	//páros eset
	{
		q[np/2] = 0.00;
	}	*/
/*	for(int i = 0; i < np; i++)
	{
		mu_q = mu_q + q[i];
	}
	mu_q = mu_q / np;	*/
//	bűvészkedés, a kódban Rmax = 0.800, ez paraméterként adódik át a függvénynek. (Re = 2.000)
	double allk = de[ns1-1]/Rmax;
/*	for(int i = 0; i < ns1; i++)
	{
		mu_de = mu_de + de[i];
	}
	mu_de = mu_de / ns1;	*/
//	ezzel nem csinálok semmit
	double ppot = 0.0;
	for(int i = 0; i < nt; i++)
	{
		ppot = ppot + pot[i];
	}
//	ppot = sqrt(ppot);
//	ami itt ki van kommentelve az nem lényeges, ami ez után következik teljesen ugyanazt csinálja.
/*	int m = 0;
	for(int i = 0; i < ns1; i++)
	{
		for(int j = 0; j < ns1; j++)
		{
			for(int k = 0; k < np; k++)
			{
//				jacbolbelso2_r1(tm3atom,tm2atom,de[i]/allk+Re,de[j]/allk+Re,q[k],var1[l],var2[l],var3[l]);
//				jacbolbelso2_r2(tm3atom,tm2atom,de[i]/allk+Re,de[j]/allk+Re,q[k],var1[l],var2[l],var3[l]);
//				jacbolbelso2_theta(tm3atom,tm2atom,de[i]/allk+Re,de[j]/allk+Re,q[k],var1[l],var2[l],var3[l]);
				jacbolbelso2(tm3atom,tm2atom,de[i]/allk+Re,de[j]/allk+Re,q[k],var1[m],var2[m],var3[m]);
				mu_r1 = mu_r1 + var1[m] - pot[m];
				mu_r2 = mu_r2 + var2[m] - pot[m];
				mu_theta = mu_theta + var3[m] - pot[m];
				if(m != nt - 1)
				{
					m++;
				}
				else
				{
					m = 0;
				}
			}
		}
	}	*/
//	a Hamilton-mátrix felépítése ugyanígy néz ki. azért írtam át erre, hogy esetleg itt-ott figyelembe lehessen venni a mátrix bizonyos részeit (ahol itt vannak kikommentelve dolgok ott)
	int m = 0, ki = 0;
	for(int km1 = 0; km1 <= 0; km1++)
	{
		for(int i = 0; i < ns1; i++)
		{
			for(int j = 0; j < ns1; j++)
			{
				for(int l = 0; l < np; l++)
				{
					for(int iv = 0; iv < ns1; iv++)
					{
						if(iv == i)
						{
							for(int jv = 0; jv < ns1; jv++)
							{
								if(jv == j)
								{
									for(int lv = 0; lv < np; lv++)
									{
										if(lv == l)
										{
											jacbolbelso2(tm3atom,tm2atom,de[i]/allk+Re,de[j]/allk+Re,q[l],var1[m],var2[m],var3[m]);
//											a - pot[m]-et ki lehet hagyni
//											mu_r1[m] = var1[m] - pot[m];
//											mu_r2[m] = var2[m] - pot[m];
//											mu_theta[m] = var3[m] - pot[m];
											if(m != nt - 1)
											{
												m++;
											}
										} 
										else
										{   
//											gscal3[m] = gscal3[m] + val[ki];
//											mu_r1 = mu_r1 - val[ki];
//											mu_r2 = mu_r2 - val[ki];
//											mu_theta = mu_theta - val[ki];
										}
										ki++;
									} // endfor lv
								}
								else
								{  
//									gscal2[m] = gscal2[m] + val[ki];
//									mu_r1 = mu_r1 - val[ki];
//									mu_r2 = mu_r2 - val[ki];
//									mu_theta = mu_theta - val[ki];
									ki++;
								}
							} // endfor jv
						} 
						else
						{
//							gscal1[m] = gscal1[m] + val[ki];
//							mu_r1 = mu_r1 - val[ki];
//							mu_r2 = mu_r2 - val[ki];
//							mu_theta = mu_theta - val[ki];
							ki++;
						}
					}// endfor iv
				} //endfor l
			}//endfor j
		}// endfor i
	}
//	itt azzal probalkoztam, hogy a lenti számtani közép helyett inkább a mediánt lenne-e jobb használni.
/*	bubblesort(nt,var1,gscal1);
	bubblesort(nt,var2,gscal2);
	bubblesort(nt,var3,gscal3);
	if(ns1 % 2 == 0)
	{
		mu_r1 = (gscal1[ns1/2] + gscal1[ns1/2+1]) / 2.0;
		mu_r2 = (gscal2[ns1/2] + gscal2[ns1/2+1]) / 2.0;
	}
	if(ns1 % 2 == 0)
	{
		mu_theta = (gscal3[np/2] + gscal3[np/2+1]) / 2.0;
	}
	if(ns1 % 2 != 0)
	{
		double p1 = (double)ns1/2.0 + 0.5;
		cout << "p1 = " << (int)p1 << endl;
		mu_r1 = gscal1[(int)p1];
		mu_r2 = gscal2[(int)p1];
	}
	if(np % 2 != 0)
	{
		double p2 = (double)np/2.0 + 0.5;
		cout << "p2 = " << (int)p2 << endl;
		mu_theta = gscal3[(int)p2];
	}	*/
//	átlagolások	
/*	for(int i = 0; i < nt; i++)
	{
		sum_mu_r1 = sum_mu_r1 + mu_r1[i];
	}
	sum_mu_r1 = sum_mu_r1 / nt;	*/
//	a sigma = 0.2 magic paraméter helyett (lásd fent) lehet számítani értelmesen is a sigma-t, de ha ezt csinálom akkor nem valami jó a kiinduláűsi vektor
	sigma_r1 = 0.0;
	for(int i = 0; i < nt; i++)
	{
		sigma_r1 = sigma_r1 + (var1[i] - sum_mu_r1) * (var1[i] - sum_mu_r1);
	}
	sigma_r1 = sqrt(sigma_r1 / (nt - 1));
//	cout << sigma_r1 << endl;
/*	for(int i = 0; i < nt; i++)
	{
		sum_mu_r2 = sum_mu_r2 + mu_r2[i];
	}
	sum_mu_r2 = sum_mu_r2 / nt;	*/
	sigma_r2 = 0.0;
	for(int i = 0; i < nt; i++)
	{
		sigma_r2 = sigma_r2 + (var2[i] - sum_mu_r2) * (var2[i] - sum_mu_r2);
	}
	sigma_r2 = sqrt(sigma_r2 / (nt - 1));
//	cout << sigma_r2 << endl;
/*	for(int i = 0; i < nt; i++)
	{
		sum_mu_theta = sum_mu_theta + mu_theta[i];
	}
	sum_mu_theta = sum_mu_theta / nt;	*/
	sigma_theta = 0.0;
	for(int i = 0; i < nt; i++)
	{
		sigma_theta = sigma_theta + (var3[i] - sum_mu_theta) * (var3[i] - sum_mu_theta);
	}
	sigma_theta = sqrt(sigma_theta / (nt - 1));
//	cout << sigma_theta << endl;
//	randomnum(nt,gscal1,mu_r1,sigma,pot);
//	randomnum(nt,gscal2,mu_r2,sigma,pot);
//	randomnum(nt,gscal3,mu_theta,sigma,pot);
//	normalRandom(nt,scal1);
//	GAUSSOK SZÁMOLÁSA ÁTLAGOLT KOORDINÁTA ÉRTÉKEK ESETÉN
	m = 0;
	double sign = 1.0;
	for(int i = 0; i < ns1; i++)
	{
//		sign = -1.0 * sign;
		for(int j = 0; j < ns1; j++)
		{
//			sign = -1.0 * sign;
			for(int k = 0; k < np; k++)
			{
				sign = -1.0 * sign;
//				jacbolbelso2(tm3atom,tm2atom,de[i],de[j],q[k],var1,var2,var3);
//				sigma_r1 = sqrt((var1[m] - sum_mu_r1) * (var1[m] - sum_mu_r1));
				f_x[m] = gauss1(sign,scal,3.543049794,var1[m],sum_mu_r1,0.2,pot[m]);
//				sigma_r2 = sqrt((var2[m] - sum_mu_r2) * (var2[m] - sum_mu_r2));
				f_y[m] = gauss1(sign,scal,3.543049794,var2[m],sum_mu_r2,0.2,pot[m]);
//				sigma_theta = sqrt((var3[m] - sum_mu_theta) * (var3[m] - sum_mu_theta));
				f_z[m] = gauss1(sign,scal,3.875686461,var3[m],sum_mu_theta,0.2,pot[m]);
				g[m] = f_x[m] * f_y[m] * f_z[m];
				if(m != nt - 1)
				{
					m++;
				}
				else
				{
					m = 0;
				}
			}
		}
	}
//	GAUSSOK SZÁMOLÁSA AZ EGYES KOORDINÁTA ÉRTÉKEK FIGYELEMBEVÉTELÉVEL, HA KELL AZ ELŐZŐ RÉSZ UNCOMMENT, EZ MEG VISSZA-COMMENT
/*	m = 0;
	double sign = 1.0;
	for(int i = 0; i < ns1; i++)
	{
//		sign = -1.0 * sign;
		for(int j = 0; j < ns1; j++)
		{
//			sign = -1.0 * sign;
			for(int k = 0; k < np; k++)
			{
				sign = -1.0 * sign;
//				jacbolbelso2(tm3atom,tm2atom,de[i],de[j],q[k],var1,var2,var3);
				sigma_r1 = sqrt((var1[m] - mu_r1[m]) * (var1[m] - mu_r1[m]));
				f_x[m] = gauss2(sign,scal,1.0/gscal1[m],var1[m],mu_r1[m],sigma_r1,pot[m]);
				sigma_r2 = sqrt((var2[m] - mu_r2[m]) * (var2[m] - mu_r2[m]));
				f_y[m] = gauss2(sign,scal,1.0/gscal2[m],var2[m],mu_r2[m],sigma_r2,pot[m]);
				sigma_theta = sqrt((var3[m] - mu_theta[m]) * (var3[m] - mu_theta[m]));
				f_z[m] = gauss2(sign,scal,1.0/gscal3[m],var3[m],mu_theta[m],sigma_theta,pot[m]);
				g[m] = f_x[m] * f_y[m] * f_z[m];
				if(m != nt - 1)
				{
					m++;
				}
				else
				{
					m = 0;
				}
			}
		}
	}	*/
//	normálás, majd kiiratás a thing.dat fájlba.
	norm = 0.0;
	for(int i = 0; i < nt; i++)
	{
		norm = norm + g[i] * g[i];
	}
	norm = sqrt(norm);
	for(int i = 0; i < nt; i++)
	{
		g[i] = (1.0 / norm) * g[i];
	}
/*	ifstream egzaktwf;
	egzaktwf.open ("sv.txt");
	egzaktwf.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < nt; i++)
	{
		egzaktwf >> pot[i];
	}
	egzaktwf.close();
	double lambda = 0.0;
	for(int i = 0; i < nt; i++)
	{
		lambda = lambda + g[i] * pot[i];
	}
	cout << lambda << endl;	*/
	ofstream ujgauss;
//	ofstream ujgauss2;
	ujgauss.open ("thing.dat");
	ujgauss.setf(ios_base::fixed,ios_base::floatfield);
//	ujgauss2.open ("thing.txt");
//	ujgauss2.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < nt; i++)
	{
		ujgauss << setprecision(30) << g[i] << endl;
//		ujgauss2 << setprecision(30) << g[i] << endl;
	}
	ujgauss.close();
//	ujgauss2.close();
	return;
}

double gauss1(double sign, double scal, double gscal, double var, double mu, double sigma, double pot)
{
	double coeff = 0.0, kitevo = 0.0, y = 0.0;
	const double pi = 3.14159265358979323846;
	
	coeff = 1.0 / (sqrt(2.0 * pi) * sigma);
	kitevo = -((var - mu) * (var - mu));
	kitevo = kitevo / (2.0 * sigma * sigma);
	y = sign * scal * coeff * exp(kitevo);
//	y = y * pot;
	
	return y;
}

double gauss2(double sign, double scal, double gscal, double var, double mu, double sigma, double pot)
{
	double coeff = 0.0, kitevo = 0.0, y = 0.0;
	const double pi = 3.14159265358979323846;
	
	coeff = gscal / sqrt(pi);
	kitevo = -(gscal * (var - mu) * gscal * (var - mu));
	kitevo = kitevo / 2.0;
	y = sign * scal * coeff * exp(kitevo);
//	y = y * pot;
	
	return y;
}

/*
void jacbolbelso2_r1(double tm3,double tm2,double R1j,double R2j,double xcosj,double &R1b,double &R2b,double &xcosb)
{
	R1b=R1j;
}

void jacbolbelso2_r2(double tm3,double tm2,double R1j,double R2j,double xcosj,double &R1b,double &R2b,double &xcosb)
{
	R2b = sqrt((tm2/(tm3+tm2))*(tm2/(tm3+tm2))*R1j*R1j+R2j*R2j+2*tm2/(tm3+tm2)*R1j*R2j*xcosj);
}

void jacbolbelso2_theta(double tm3,double tm2,double R1j,double R2j,double xcosj,double &R1b,double &R2b,double &xcosb)
{
	xcosb = R2j/R2b*xcosj+tm2/(tm3+tm2)*R1j/R2b;
}
*/

void randomnum(int nt, double *scal, double mu, double sigma, double *pot)
{
	default_random_engine generator;
	normal_distribution<double> distribution(mu,0.00001);
	
	for(int i = 0; i < nt; i++)
	{
		scal[i] = distribution(generator);
		cout << scal[i] << endl;
	}
}

double uniformRandom()
{
	return ( (double)(rand()) + 1.0 )/( (double)(RAND_MAX) + 1.0 );
}
 // return a normally distributed random number
void normalRandom(int nt, double *scal)
{
	for(int i = 0; i < nt; i++)
	{
		double u1=uniformRandom();
		double u2=uniformRandom();
		scal[i] = cos(8.0*atan(1.0)*u2)*sqrt(-2.0*log(u1));
	}
}

void DirectIteration(int n, int nonzero, int maxiter, double *val, int *ncol, int *nrow, double *ud, int *diagind)
{
	int k;
	double lambda;
	
	double *psi = new double[n];
	double *Hpsi = new double[n];
	double *val2 = new double[nonzero];
	double *zero = new double[n];
	for(int i = 1; i <= maxiter; i++)
	{
		zero[i] = 0.00;
	}
	
	ifstream wf;
	wf.open("thing.dat");
	wf.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 0; i < n; i++)
	{
		wf >> psi[i];
	}
	wf.close();
	
	ofstream file;
	ofstream file2;		// CSAK, HOGY MEGFELELŐ SZÁMÚ CGM PARAMÉTER LEGYEN ÁTADVA
	file.open("shiet.txt");
	file.setf(ios_base::fixed,ios_base::floatfield);
	for(int i = 1; i <= maxiter; i++)
	{
		mult(n,psi,Hpsi,val,ncol,nrow);
		lambda = 0.0;
		for(int i = 0; i < n; i++)
		{
			lambda = lambda + psi[i] * Hpsi[i];
		}
		cout << "lambda = " << lambda << endl;
		k = 0;
		for(int j = 0; j < nonzero; j++)
		{
			if(j == diagind[k])
			{
				val2[j] = val[j] - lambda;
				k++;
			}
			else
			{
				val2[j] = val[j];
			}
		}
		CGM(n,100000,pow(10.0,-5.0),val2,ncol,nrow,ud,diagind,psi,zero,file,3,file2);
	}
	file.close();
}
