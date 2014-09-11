#include <iostream>
#include <fstream>
#include <complex>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <omp.h>

using namespace std;

#define PI 3.1415927
#define QE 1.602177e-19
#define E0 8.85419e-12
#define U0 1.2566e-6
#define ME 9.10939e-31
#define MI 5.0e-26
#define KB 1.38066e-23
#define C 2.9979246e8
#define DATATYPE double

// above MI is approximate mass of NO+ ion

//////////////////////

// EMP 1D version: created from 2D code on February 27, 2014


///////////////////////////////////////////////////////////////////////////////////

int main()
{

  // time the whole thing
  double tStart = omp_get_wtime();

  //------------------------------------------------------------------------//
  //-------READ INPUTS FROM FILE -------------------------------------------//

  FILE * inputFile;
  inputFile = fopen("inputs.dat","rb");

  int doionosphere;
  int doioniz;
  double maxalt;
  double stepalt;
  double dx1;
  double dx2;
  double dt;
  int tsteps;
  int numfiles;
  int nprobes;

  fread(&doionosphere,sizeof(int),1,inputFile);
  fread(&doioniz,sizeof(int),1,inputFile);
  fread(&maxalt,sizeof(double),1,inputFile);
  fread(&stepalt,sizeof(double),1,inputFile);
  fread(&dx1,sizeof(double),1,inputFile);
  fread(&dx2,sizeof(double),1,inputFile);
  fread(&dt,sizeof(double),1,inputFile);
  fread(&tsteps,sizeof(int),1,inputFile);

  double Ezsource [tsteps];
  fread(&Ezsource,sizeof(double),tsteps,inputFile);
  fread(&numfiles,sizeof(int),1,inputFile);
  fread(&nprobes,sizeof(int),1,inputFile);
  int probex [nprobes];
  fread(&probex,sizeof(int),nprobes,inputFile);
  fclose(inputFile);

  cout << "made it to here!\n";

  //------------------------------------------------------------------------//
  //-------GRID PARAMETERS -------------------------------------------------//

  // set up radial array

  int xx = stepalt/dx1 + (maxalt - stepalt)/dx2 + 1;
  DATATYPE x [xx];
  DATATYPE dx [xx-1];
  x[0] = 0;
  for (int i = 1; i < xx; i++) {
    if (x[i-1] < stepalt) { 
      x[i] = x[i-1] + dx1;
    } else { 
      x[i] = x[i-1] + dx2; 
    }
    dx[i-1] = x[i] - x[i-1];
  }
  
  //------------------------------------------------------------------------//
  //-------MULTIPLY FACTORS FOR FDTD ---------------------------------------//
  
  double sigm = 0;
  double sig = 0;

  DATATYPE c1h = (2*U0 - sigm*dt)/(2*U0 + sigm*dt);
  DATATYPE c2h = 2*dt/(2*U0 + sigm*dt);
  
  DATATYPE c1e = (2*E0 - sig*dt)/(2*E0 + sig*dt);
  DATATYPE c2e = 2*dt/(2*E0 + sig*dt);
  
  //------------------------------------------------------------------------//
  //-------MAGNETIC FIELD --------------------------------------------------//

  // read magnetic field: varying in theta direction
  double Bx, By, Bz;
  FILE * BFile;
  BFile = fopen("B0.dat","rb");
  fread(&Bx,sizeof(double),1,BFile);
  fread(&By,sizeof(double),1,BFile);
  fread(&Bz,sizeof(double),1,BFile);
  fclose(BFile);

  // compute gyrofrequencies

  DATATYPE wcex = -QE * Bx / ME;
  DATATYPE wcey = -QE * By / ME;
  DATATYPE wcez = -QE * Bz / ME;
  DATATYPE wce0 = sqrt(wcex*wcex + wcey*wcey + wcez*wcez);
  DATATYPE wcix = QE * Bx / MI;
  DATATYPE wciy = QE * By / MI;
  DATATYPE wciz = QE * Bz / MI;
  DATATYPE wci0 = sqrt(wcix*wcix + wciy*wciy + wciz*wciz);

  //------------------------------------------------------------------------//
  //-------OUTPUT TO LOG FILE ----------------------------------------------//

  ofstream logfile;
  logfile.open("log.txt");

  logfile << "\n";
  logfile << "------------------------------------------------------------------\n";
  logfile << " ------------------- 1D EMP Simulation ---------------------------\n";
  logfile << "DoIoniz = " << doioniz << "\n";
  logfile << "Grid is " << xx << " (in x) cells\n";
  logfile << "Bx = " << Bx << "; By = " << By << "; Bz = " << Bz << "\n";
  logfile << "Time step is " << (dt*1e9) << " ns\n";
  logfile << "Simulation will run " << tsteps << " time steps.\n";

  // write space parameters to sferic.dat file

  FILE * sfericFile;

  sfericFile = fopen("sferic.dat","wb");
  fwrite(&tsteps,sizeof(int),1,sfericFile);
  fwrite(&xx,sizeof(int),1,sfericFile);
  fwrite(&numfiles,sizeof(int),1,sfericFile);
  fwrite(&dt,sizeof(DATATYPE),1,sfericFile);
  fwrite(x,sizeof(DATATYPE),xx,sfericFile);
  fclose(sfericFile);


  //------------------------------------------------------------------------//
  //------- FIELD INITIALIZATION --------------------------------------------//

  // electric field vector
  DATATYPE Ex [xx];
  DATATYPE Ey [xx];
  DATATYPE Ez [xx];
  // magnetic field vector
  DATATYPE Hx [xx];
  DATATYPE Hy [xx];
  DATATYPE Hz [xx];
  // spatial-averaged electric field
  DATATYPE Exm, Hym, Hzm;

  // electron current
  DATATYPE Jex [xx];
  DATATYPE Jey [xx];
  DATATYPE Jez [xx];
  DATATYPE Jex0 = 0;
  DATATYPE Jey0 = 0;

  // ion current
  DATATYPE Jix [xx];
  DATATYPE Jiy [xx];
  DATATYPE Jiz [xx];
  DATATYPE Jix0 = 0;
  DATATYPE Jiy0 = 0;

  // total current
  DATATYPE Jx [xx];
  DATATYPE Jy [xx];
  DATATYPE Jz [xx];

  // Electric field Parallel and Perpendicular components, for time-dependent stuff
  DATATYPE Epar, Eperp2;
  DATATYPE Eeff [xx];
  DATATYPE Emag [xx];

  // probe fields
  DATATYPE Exprobe [tsteps][nprobes];
  DATATYPE Eyprobe [tsteps][nprobes];
  DATATYPE Ezprobe [tsteps][nprobes];
  DATATYPE Hxprobe [tsteps][nprobes];
  DATATYPE Hyprobe [tsteps][nprobes];
  DATATYPE Hzprobe [tsteps][nprobes];
  int ix;

  // initialize

  memset(Ex,0,sizeof(Ex));
  memset(Ey,0,sizeof(Ey));
  memset(Ez,0,sizeof(Ez));
  memset(Hx,0,sizeof(Hx));
  memset(Hy,0,sizeof(Hy));
  memset(Hz,0,sizeof(Hz));
  memset(Jex,0,sizeof(Jex));
  memset(Jey,0,sizeof(Jey));
  memset(Jez,0,sizeof(Jez));
  memset(Jix,0,sizeof(Jix));
  memset(Jiy,0,sizeof(Jiy));
  memset(Jiz,0,sizeof(Jiz));
  memset(Jx,0,sizeof(Jx));
  memset(Jy,0,sizeof(Jy));
  memset(Jz,0,sizeof(Jz));
  memset(Eeff,0,sizeof(Eeff));
  memset(Emag,0,sizeof(Emag));
  
  //------------------------------------------------------------------------//
  //-------NUMBER DENSITIES ------------------------------------------------//

  // set up electron density and B field

  double ne [xx];        // 1D electron density
  DATATYPE wpe [xx];     // plasma frequency
  double nd [xx];        // neutral density
  DATATYPE mue;          // mobility
  DATATYPE nue [xx];     // collision frequency
  DATATYPE Ek [xx];      // breakdown field

  // do Ions as well!
  double ni [xx];
  DATATYPE wpi [xx];
  DATATYPE nui [xx];

  FILE * neFile;
  neFile = fopen("ne.dat","rb");
  fread(&ne,sizeof(double),xx,neFile);
  fclose(neFile);

  FILE * niFile;
  niFile = fopen("ni.dat","rb");
  fread(&ni,sizeof(double),xx,niFile);
  fclose(niFile);

  for (int i = 0; i < xx; i++) {
    //ne[i] = 0;
    //ni[i] = 0;
    wpe[i] = QE * sqrt(ne[i] / (ME * E0));
    wpi[i] = QE * sqrt(ni[i] / (MI * E0));
  }
  
  // temperature intermediates; private variables in openmp, except S and components
  
  DATATYPE S [xx];
  DATATYPE Sx,Sy,Sz;
  DATATYPE JdotE [xx];
  DATATYPE heat [xx];

  memset(S,0,sizeof(S));
  memset(JdotE,0,sizeof(JdotE));
  memset(heat,0,sizeof(heat));

  //------------------------------------------------------------------------//
  //-------OPTICAL EXCITATION RATES ----------------------------------------//

  // read and compute rate arrays to be used for nonlinear stuff

  FILE * rateFile;
  rateFile = fopen("rates.dat","rb");

  // code uses lookup table / interpolation, these are interpolation params
  int topE, botE;
  int nume;
  fread(&nume,sizeof(int),1,rateFile);
  double efield [nume];
  fread(&efield,sizeof(double),nume,rateFile);

  double ioniz [xx*nume];
  double attach [xx*nume];
  double mobility [xx*nume];
  double Ored [xx*nume];
  double Ogrn [xx*nume];
  double N21p [xx*nume];
  double N22p [xx*nume];
  double N2p1n [xx*nume];
  double N2pM [xx*nume];
  double O2p1n [xx*nume];

  fread(&ioniz,sizeof(double),xx*nume,rateFile);
  fread(&attach,sizeof(double),xx*nume,rateFile);
  fread(&mobility,sizeof(double),xx*nume,rateFile);
  fread(&Ored,sizeof(double),xx*nume,rateFile);
  fread(&Ogrn,sizeof(double),xx*nume,rateFile);
  fread(&N21p,sizeof(double),xx*nume,rateFile);
  fread(&N22p,sizeof(double),xx*nume,rateFile);
  fread(&N2p1n,sizeof(double),xx*nume,rateFile);
  fread(&N2pM,sizeof(double),xx*nume,rateFile);
  fread(&O2p1n,sizeof(double),xx*nume,rateFile);

  fclose(rateFile);

  // reorder them into a useful 2D array

  double viArray [xx][nume];
  double vaArray [xx][nume];
  double muArray [xx][nume];
  double OrArray [xx][nume];
  double OgArray [xx][nume];
  double N21pArray [xx][nume];
  double N22pArray [xx][nume];
  double N2p1nArray [xx][nume];
  double N2pMArray [xx][nume];
  double O2p1nArray [xx][nume];

  int index;

  for (int i = 0; i < xx; i++) {
    for (int m = 0; m < nume; m++) {
      index = xx*m + i;
      viArray[i][m] = ioniz[index];
      vaArray[i][m] = attach[index];
      muArray[i][m] = mobility[index];
      OrArray[i][m] = Ored[index];
      OgArray[i][m] = Ogrn[index];
      N21pArray[i][m] = N21p[index];
      N22pArray[i][m] = N22p[index];
      N2p1nArray[i][m] = N2p1n[index];
      N2pMArray[i][m] = N2pM[index];
      O2p1nArray[i][m] = O2p1n[index];
    }
  }

  // array of five values that are used to solve the ionization / detachment / attachment equations
  // will use this as a private variable in updates.
  // first four are matrix entries; fifth is determinant.
 
  DATATYPE Aion [5];
  DATATYPE nenew;


  //------------------------------------------------------------------------//
  //------- BREAKDOWN, MOBILITY, ETC----------------------------------------//

  FILE * ndFile;
  ndFile = fopen("nd.dat","rb");
  fread(&nd,sizeof(double),xx,ndFile);
  fclose(ndFile);

  for (int i = 0; i < xx; i++) {
    Ek[i] = 3.0e6 * nd[i] / nd[0];
    mue = muArray[i][0] / nd[i];
    nue[i] = (QE / ME) / mue;
    nui[i] = nue[i] / 100; //* ME / MI;    // wrong: fix later
  }
  
  //------------------------------------------------------------------------//
  //------- OPTICAL FIELDS -------------------------------------------------//

  DATATYPE EoEk;
  DATATYPE vi, va, vd;  // ionization
  DATATYPE nOm [xx]; // O- ions

  // optical arrays
  DATATYPE vN21P, vN22P, vN2P1N, vN2PM, vO2P1N, vOred, vOgrn;
  DATATYPE nN21P [xx];
  DATATYPE nN22P [xx];
  DATATYPE nN2P1N [xx];
  DATATYPE nN2PM [xx];
  DATATYPE nO2P1N [xx];
  DATATYPE nOred [xx];
  DATATYPE nOgrn [xx];
  DATATYPE tauN22P, tauN21P, tauN2P1N, tauN2PM, tauO2P1N, tauOred, tauOgrn;

  // initialize optical

  memset(nN21P,0,sizeof(nN21P));
  memset(nN22P,0,sizeof(nN22P));
  memset(nN2P1N,0,sizeof(nN2P1N));
  memset(nN2PM,0,sizeof(nN2PM));
  memset(nO2P1N,0,sizeof(nO2P1N));
  memset(nOred,0,sizeof(nOred));
  memset(nOgrn,0,sizeof(nOgrn));
  memset(nOm,0,sizeof(nOm));


  //------------------------------------------------------------------------//
  //------- A COUPLE OF DETAILS --------------------------------------------//

  // partialtime used to determine write times
  DATATYPE partialtime;

  // these are A and K matrices for Lee&Kalluri current solution
  DATATYPE Ae [3][3];
  DATATYPE Ke [3][3];
  DATATYPE Ai [3][3];
  DATATYPE Ki [3][3];
  DATATYPE Emidx;

  // intermediate calculations in Lee&Kalluri
  DATATYPE Ee1, Ee2, Se1, Ce1, Ce2, Ce3, Ce4;
  DATATYPE Ei1, Ei2, Si1, Ci1, Ci2, Ci3, Ci4;


// FEB 27 12:20 PM - I AM HERE

  //------------------------------------------------------------------------//
  //------- SET UP PML -----------------------------------------------------//

  // initialize A, B, and P arrays even if we don't use them.

  double pmlm = 4;
  int pmllen = 10;

  DATATYPE A_eyx [pmllen];
  DATATYPE A_ezx [pmllen];
  DATATYPE A_hyx [pmllen];
  DATATYPE A_hzx [pmllen];

  DATATYPE B_eyx [pmllen];
  DATATYPE B_ezx [pmllen];
  DATATYPE B_hyx [pmllen];
  DATATYPE B_hzx [pmllen];

  DATATYPE P_eyx [pmllen];
  DATATYPE P_ezx [pmllen];
  DATATYPE P_hyx [pmllen];
  DATATYPE P_hzx [pmllen];

  // initialize
  memset(P_eyx,0,sizeof(P_eyx));
  memset(P_ezx,0,sizeof(P_ezx));
  memset(P_hyx,0,sizeof(P_hyx));
  memset(P_hzx,0,sizeof(P_hzx));

  // end initialization

  DATATYPE simax = 1.5 * (DATATYPE)(pmlm + 1) / (150 * PI * dx[xx-2]);
  DATATYPE kamax = 1;
  DATATYPE almax = 0;
  
  DATATYPE sx [pmllen];
  DATATYPE sxm [pmllen];
  DATATYPE kx [pmllen];
  DATATYPE kxm [pmllen];
  DATATYPE ax [pmllen];
  DATATYPE axm [pmllen];
  DATATYPE mf;
  
  for (int m = 0; m < pmllen; m++) {
    // note shifts here compared to matlab, due to zero index
    mf = (DATATYPE)m;
    sx[m] = simax * pow(((mf+0.5)/pmllen),pmlm);
    sxm[m] = simax * pow(((mf+1)/pmllen),pmlm);
    kx[m] = 1 + (kamax-1) * pow(((mf+0.5)/pmllen),pmlm);
    kxm[m] = 1 + (kamax-1) * pow(((mf+1)/pmllen),pmlm);
    ax[m] = almax * pow(((pmllen-mf-0.5)/pmllen),pmlm);
    axm[m] = almax * pow(((pmllen-mf-1)/pmllen),pmlm);
  }
  
  // update / initialize A and B vectors
  for (int m = 0; m < pmllen; m++) {
    B_eyx[m] = exp(-((sx[m]/kx[m]) + ax[m])*dt/E0);
    A_eyx[m] = sx[m] / (sx[m]*kx[m] + pow(kx[m],2)*ax[m]) * (B_eyx[m] - 1);
    B_ezx[m] = exp(-((sx[m]/kx[m]) + ax[m])*dt/E0);
    A_ezx[m] = sx[m] / (sx[m]*kx[m] + pow(kx[m],2)*ax[m]) * (B_ezx[m] - 1);
    
    B_hyx[m] = exp(-((sxm[m]/kxm[m]) + axm[m])*dt/E0);
    A_hyx[m] = sxm[m] / (sxm[m]*kxm[m] + pow(kxm[m],2)*axm[m]) * (B_hyx[m] - 1);
    B_hzx[m] = exp(-((sxm[m]/kxm[m]) + axm[m])*dt/E0);
    A_hzx[m] = sxm[m] / (sxm[m]*kxm[m] + pow(kxm[m],2)*axm[m]) * (B_hzx[m] - 1);     
    
  }
  
  // offsets for PML
  int xshift = xx-pmllen-1;
  
  // measure memory usage, just for kicks
  struct rusage ru;
  double memusage;
  getrusage(RUSAGE_SELF, &ru);
  memusage = (double)ru.ru_maxrss / 1024;

  logfile << "Total memory usage = " << memusage << " MB\n";
  logfile.close();


  ////////////////////////////////////////////////////////////////////////////
  //                                                                      ////
  //------- MAIN TIME LOOP -----------------------------------------------////

  double tloopStart = omp_get_wtime();

  for (int t = 0; t < tsteps; t++) {

    // ---------------------------------------------------
    // Psi updates for H field
    for (int i = 0; i < pmllen; i++) {
      P_hyx[i] = B_hyx[i] * P_hyx[i]					\
	+ A_hyx[i] * ( Ez[i+xshift+1] - Ez[i+xshift] ) / dx[i+xshift];
    }
    for (int i = 0; i < pmllen; i++) {
      P_hzx[i] = B_hzx[i] * P_hzx[i]					\
	+ A_hzx[i] * ( Ey[i+xshift+1] - Ey[i+xshift] ) / dx[i+xshift];
    }
    

    // ----------------------------------------
    // Hx update: none, since there is no dy or dz!

    // ----------------------------------------
    // Hy update
#pragma omp parallel for
    for (int i = 0; i < xx-1; i++) {
      Hy[i] = c1h * Hy[i] + c2h * ( Ez[i+1] - Ez[i] ) / dx[i];
    }
    // pml correction
    for (int i = 0; i < pmllen; i++) {
      Hy[i+xshift] = Hy[i+xshift] + c2h * P_hyx[i];
    }

    // ----------------------------------------
    // Hz update
#pragma omp parallel for
    for (int i = 0; i < xx-1; i++) {
      Hz[i] = c1h * Hz[i] - c2h * ( Ey[i+1] - Ey[i] ) / dx[i];
    }
    // pml correction
    for (int i = 0; i < pmllen; i++) {
      Hz[i+xshift] = Hz[i+xshift] - c2h * P_hzx[i];
    }

    // ---------------------------------------
    // Psi updates for E field

      for (int i = 0; i < pmllen; i++) {
	  P_eyx[i] = B_eyx[i] * P_eyx[i] + A_eyx[i] * ( Hz[i+xshift] - Hz[i+xshift-1] ) / ((dx[i+xshift]+dx[i+xshift-1])/2);
      }
      for (int i = 0; i < pmllen; i++) {
	  P_ezx[i] = B_ezx[i] * P_ezx[i] + A_ezx[i] * ( Hy[i+xshift] - Hy[i+xshift-1] ) / ((dx[i+xshift]+dx[i+xshift-1])/2);
      }

    // ----------------------------------------
    // Ex update: there is no del x H for Ex, but there is Jx
#pragma omp parallel for
    for (int i = 1; i < xx-1; i++) {
      Ex[i] = c1e * Ex[i] - c2e * (Jx[i+1] + Jx[i]) / 2;
    }

    // ----------------------------------------
    // Ey update
#pragma omp parallel for
    for (int i = 1; i < xx-1; i++) {
	Ey[i] = c1e * Ey[i] - c2e * ( Hz[i] - Hz[i-1] ) / ((dx[i]+dx[i-1])/2)  - c2e * Jy[i];
    }
    // pml correction
    for (int i = 0; i < pmllen; i++) {
      Ey[i+xshift] = Ey[i+xshift] - c2e * P_eyx[i];
    }

    // Ey source
    Ey[0] = Ezsource[t];

    // ----------------------------------------
#pragma omp parallel for
    // Ez update
    for (int i = 1; i < xx-1; i++) {
	Ez[i] = c1e * Ez[i] + c2e * ( Hy[i] - Hy[i-1] ) / ((dx[i]+dx[i-1])/2) - c2e * Jz[i];
    }
    // pml corrections
    for (int i = 0; i < pmllen; i++) {
      Ez[i+xshift] = Ez[i+xshift] + c2e * P_ezx[i];
    }

    // Ez source
    Ez[0] = Ezsource[t];

    // output E and H probe points at desired points

    for (int i = 0; i < nprobes; i++) {
      ix = probex[i];
      Exprobe[t][i] = 0.5*(Ex[ix]+Ex[ix-1]);
      Eyprobe[t][i] = Ey[ix];
      Ezprobe[t][i] = Ez[ix];

      Hxprobe[t][i] = Hx[ix];
      Hyprobe[t][i] = 0.5*(Hy[ix]+Hy[ix-1]);
      Hzprobe[t][i] = 0.5*(Hz[ix]+Hz[ix-1]);
    }

    // ---------------------------------------
    // update Eeff, for use in ionization updates

    if (doionosphere) {
#pragma omp parallel for private(Epar,Eperp2,Sx,Sy,Sz,Exm,Hym,Hzm)
       for (int i = 0; i < xx-1; i++) {
 
 	  Exm = ( dx[i-1] * Ex[i] + dx[i] * Ex[i-1] ) / (dx[i] + dx[i-1]);
 	  Hym = ( dx[i-1] * Hy[i] + dx[i] * Hy[i-1] ) / (dx[i] + dx[i-1]);
 	  Hzm = ( dx[i-1] * Hz[i] + dx[i] * Hz[i-1] ) / (dx[i] + dx[i-1]);
 	  Epar = (wcex/wce0)*Exm + (wcey/wce0)*Ey[i] + (wcez/wce0)*Ez[i];
 	  Emag[i] = sqrt( Exm*Exm + Ey[i]*Ey[i] + Ez[i]*Ez[i] );
 	  Eperp2 = Emag[i]*Emag[i] - Epar*Epar;
 	  Eeff[i] = sqrt( Epar*Epar + Eperp2 * nue[i]*nue[i] / ( nue[i]*nue[i] + wce0*wce0 ) );
 
 	  Sx = Ey[i]*Hzm - Ez[i]*Hym;   // this is not quite right, because H is at half-time-steps...
 	  Sy = Ez[i]*Hx[i] - Exm*Hzm;
 	  Sz = Exm*Hym - Ey[i]*Hx[i];
 	  S[i] = sqrt( Sx*Sx + Sy*Sy + Sz*Sz );
       }
 
     } // if doionosphere
     

     // Ionization updates.
     // ------------------------------------------
     // update mue and nue, calculate vi, va, ne. Not going to do it inside PML! Makes things unstable!

    if (doionosphere & doioniz) {

#pragma omp parallel for private(EoEk,vN21P,vN22P,vN2P1N,vN2PM,vO2P1N,vOred,vOgrn,tauN21P,tauN22P,tauN2P1N,tauN2PM,tauO2P1N,tauOred,tauOgrn,Aion,nenew,topE,botE,vi,va,vd,mue)
      for (int i = 1; i < xx-pmllen+1; i++) {

	topE = 1;
	botE = 0;
	// find index into Efield rate array
	for (int m = 0; m < nume; m++) {
	  if (efield[m] > Eeff[i] / nd[i]) {
	    topE = m;
	    botE = topE - 1;
	    break;
	  }
	}

	// read off mue, vi, va, vN21P, vN22P from the arrays and interpolate
	  
	if (topE <= 1) {
	  // set everything to first value 
	  mue = muArray[i][0] / nd[i];
	  vi = 0;
	  va = 0;
	  vN21P = 0;
	  vN22P = 0;
	  vN2P1N = 0;
	  vN2PM = 0;
	  vO2P1N = 0;
	  vOgrn = 0;
	  vOred = 0;
	} else {
	  //  interpolate
	  mue = (1/nd[i]) * ( muArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((muArray[i][topE]-muArray[i][botE])/(efield[topE]-efield[botE])) );
	  vi = nd[i] * ( viArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((viArray[i][topE]-viArray[i][botE])/(efield[topE]-efield[botE])) );
	  va = nd[i] * ( vaArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((vaArray[i][topE]-vaArray[i][botE])/(efield[topE]-efield[botE])) );
	  vN21P = nd[i] * ( N21pArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((N21pArray[i][topE]-N21pArray[i][botE])/(efield[topE]-efield[botE])) );
	  vN22P = nd[i] * ( N22pArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((N22pArray[i][topE]-N22pArray[i][botE])/(efield[topE]-efield[botE])) );
	  vN2P1N = nd[i] * ( N2p1nArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((N2p1nArray[i][topE]-N2p1nArray[i][botE])/(efield[topE]-efield[botE])) );
	  vN2PM = nd[i] * ( N2pMArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((N2pMArray[i][topE]-N2pMArray[i][botE])/(efield[topE]-efield[botE])) );
	  vO2P1N = nd[i] * ( O2p1nArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((O2p1nArray[i][topE]-O2p1nArray[i][botE])/(efield[topE]-efield[botE])) );
	  vOred = nd[i] * ( OrArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((OrArray[i][topE]-OrArray[i][botE])/(efield[topE]-efield[botE])) );
	  vOgrn = nd[i] * ( OgArray[i][botE] + (Eeff[i]/nd[i] - efield[botE])*((OgArray[i][topE]-OgArray[i][botE])/(efield[topE]-efield[botE])) );
	}

	// old method for detachment coefficient
	EoEk = Eeff[i] / Ek[i];
	if (EoEk < 0.05) {
	  vd = 0;
	} else {
	  vd = 0.78 * nd[i] * 1.08e-18 * sqrt(EoEk) * exp(-0.078375 / EoEk);
	}
	
	// update electron density and plasma frequency, and excited state density
	// test to make sure vi is not too big - things will go unstable! this happens inside right PML.

	if ( (vi - va + vd)*dt < 0.5 ) {
	  Aion[4] = 1 - (vi - va - vd)*dt - vi*vd*dt*dt; // determinant of solution matrix
	  Aion[0] = (1 + vd*dt);
	  Aion[1] = vd*dt;
	  Aion[2] = va*dt;
	  Aion[3] = 1 - (vi - va)*dt;
	  nenew = ( Aion[0] * ne[i] + Aion[1] * nOm[i] ) / Aion[4];     // temporary value so we don't overwrite
	  nOm[i] = ( Aion[2] * ne[i] + Aion[3] * nOm[i] ) / Aion[4];
	  ne[i] = nenew;
	} else {
	  // throttle the ionization so it doesn't go crazy unstable
	  ne[i] = ne[i] * exp(0.5);
	}

	wpe[i] = QE * sqrt( ne[i] / (ME * E0) );
	nue[i]  = (QE / ME) / mue;

	// assume that ionization creates a positive ion, and attachment creates a negative ion?
	// for now, we will assume ion densities don't change.

	// optical: Need to input relative densities of O2, N2, O versus altitude to make these equations valid at higher altitudes.

	tauN22P = 1/(2e7 + 3e-19*(0.21*nd[i]) );  // quenching through O2, which is 21% of density
	tauN21P = 1/(1.7e5 + 1e-17*(0.78*nd[i]) ); // quenching through N2, which is 78% of density
	tauN2P1N = 1/(1.4e7 + 4e-16*(0.99*nd[i]) ); // quenched by N2 and O2
	tauN2PM = 1/(7e4 + 5e-16*(0.78*nd[i]) ); // quenched by N2
	tauO2P1N = 1/(8.3e5 + 2e-16*(0.78*nd[i]) ); // quenched by N2
	tauOred = 1/(0.0091 + 5e-17*(0.78*nd[i]) );
	tauOgrn = 1/(1.43 + 3e-19*(0.21*nd[i]) );

	// attemping backward Euler (implicit) for optics, since time constants are smaller than dt in some cases.

	nN22P[i] = tauN22P/(dt + tauN22P) * nN22P[i] + (tauN22P*dt)/(dt + tauN22P) * ( vN22P * ne[i] );
	nN21P[i] = tauN21P/(dt + tauN21P) * nN21P[i] + (tauN21P*dt)/(dt + tauN21P) * ( vN21P * ne[i] + 2e7 * nN22P[i] ); // cascading
	nN2P1N[i] = tauN2P1N/(dt + tauN2P1N) * nN2P1N[i] + (tauN2P1N*dt)/(dt + tauN2P1N) * ( vN2P1N * ne[i] );
	nN2PM[i] = tauN2PM/(dt + tauN2PM) * nN2PM[i] + (tauN2PM*dt)/(dt + tauN2PM) * ( vN2PM * ne[i] );
	nO2P1N[i] = tauO2P1N/(dt + tauO2P1N) * nO2P1N[i] + (tauO2P1N*dt)/(dt + tauO2P1N) * ( vO2P1N * ne[i] );

	nOgrn[i] = tauOgrn/(dt + tauOgrn) * nOgrn[i] + (tauOgrn*dt)/(dt + tauOgrn) * ( vOgrn * ne[i] );
	nOred[i] = tauOred/(dt + tauOred) * nOred[i] + (tauOred*dt)/(dt + tauOred) * ( vOred * ne[i] + 1.43 * nOgrn[i] );

      }

    } // if doioniz


    // ----------------------------------------
    // J update equations... parallelizing these actually slows it down, since there are so many private variables
    
    if (doionosphere) {
#pragma omp parallel for private(Ee1,Ee2,Se1,Ce1,Ce2,Ce3,Ce4,Ei1,Ei2,Si1,Ci1,Ci2,Ci3,Ci4,Ae,Ke,Ai,Ki,Emidx,Jex0,Jey0,Jix0,Jiy0)
      for (int i = 1; i < xx; i++) {

	  Ee1 = exp(-nue[i]*dt);
	  Ee2 = 1/(wce0*wce0 + nue[i]*nue[i]);
	  if (wce0 == 0) {
	    Se1 = 1;
	    Ce1 = 0;
	  } else {
	    Se1 = sin(wce0*dt)/wce0;
	    Ce1 = (1 - cos(wce0*dt))/(wce0*wce0);
	  }
	  Ce2 = (1 - Ee1)/nue[i] - Ee1*nue[i]*Ce1 - Ee1*Se1;
	  Ce3 = nue[i] * (1 - Ee1*cos(wce0*dt)) + Ee1*wce0*sin(wce0*dt);
	  Ce4 = 1 - Ee1*cos(wce0*dt) - Ee1*nue[i]*Se1;

	  // same for ions
	  Ei1 = exp(-nui[i]*dt);
	  Ei2 = 1/(wci0*wci0 + nui[i]*nui[i]);
	  if (wci0 == 0) {
	    Si1 = 1;
	    Ci1 = 0;
	  } else {
	    Si1 = sin(wci0*dt)/wci0;
	    Ci1 = (1 - cos(wci0*dt))/(wci0*wci0);
	  }
	  Ci2 = (1 - Ei1)/nui[i] - Ei1*nui[i]*Ci1 - Ei1*Si1;
	  Ci3 = nui[i] * (1 - Ei1*cos(wci0*dt)) + Ei1*wci0*sin(wci0*dt);
	  Ci4 = 1 - Ei1*cos(wci0*dt) - Ei1*nui[i]*Si1;

	  // A and K matrices
	  Ae[0][0] = Ee1 * ( Ce1*wcex*wcex + cos(wce0*dt) );
	  Ae[0][1] = Ee1 * ( Ce1*wcex*wcey - Se1*wcez );
	  Ae[0][2] = Ee1 * ( Ce1*wcex*wcez + Se1*wcey );
	  Ae[1][0] = Ee1 * ( Ce1*wcey*wcex + Se1*wcez );     
	  Ae[1][1] = Ee1 * ( Ce1*wcey*wcey + cos(wce0*dt) );
	  Ae[1][2] = Ee1 * ( Ce1*wcey*wcez - Se1*wcex );
	  Ae[2][0] = Ee1 * ( Ce1*wcez*wcex - Se1*wcey );
	  Ae[2][1] = Ee1 * ( Ce1*wcez*wcey + Se1*wcex );
	  Ae[2][2] = Ee1 * ( Ce1*wcez*wcez + cos(wce0*dt) );
            
	  Ke[0][0] = Ee2 * ( Ce2*wcex*wcex + Ce3 );
	  Ke[0][1] = Ee2 * ( Ce2*wcex*wcey - Ce4*wcez );
	  Ke[0][2] = Ee2 * ( Ce2*wcex*wcez + Ce4*wcey );
	  Ke[1][0] = Ee2 * ( Ce2*wcey*wcex + Ce4*wcez );
	  Ke[1][1] = Ee2 * ( Ce2*wcey*wcey + Ce3 );
	  Ke[1][2] = Ee2 * ( Ce2*wcey*wcez - Ce4*wcex );
	  Ke[2][0] = Ee2 * ( Ce2*wcez*wcex - Ce4*wcey );
	  Ke[2][1] = Ee2 * ( Ce2*wcez*wcey + Ce4*wcex );
	  Ke[2][2] = Ee2 * ( Ce2*wcez*wcez + Ce3 );

	  // A and K for ions
	  Ai[0][0] = Ei1 * ( Ci1*wcix*wcix + cos(wci0*dt) );
	  Ai[0][1] = Ei1 * ( Ci1*wcix*wciy - Si1*wciz );
	  Ai[0][2] = Ei1 * ( Ci1*wcix*wciz + Si1*wciy );
	  Ai[1][0] = Ei1 * ( Ci1*wciy*wcix + Si1*wciz );     
	  Ai[1][1] = Ei1 * ( Ci1*wciy*wciy + cos(wci0*dt) );
	  Ai[1][2] = Ei1 * ( Ci1*wciy*wciz - Si1*wcix );
	  Ai[2][0] = Ei1 * ( Ci1*wciz*wcix - Si1*wciy );
	  Ai[2][1] = Ei1 * ( Ci1*wciz*wciy + Si1*wcix );
	  Ai[2][2] = Ei1 * ( Ci1*wciz*wciz + cos(wci0*dt) );
            
	  Ki[0][0] = Ei2 * ( Ci2*wcix*wcix + Ci3 );
	  Ki[0][1] = Ei2 * ( Ci2*wcix*wciy - Ci4*wciz );
	  Ki[0][2] = Ei2 * ( Ci2*wcix*wciz + Ci4*wciy );
	  Ki[1][0] = Ei2 * ( Ci2*wciy*wcix + Ci4*wciz );
	  Ki[1][1] = Ei2 * ( Ci2*wciy*wciy + Ci3 );
	  Ki[1][2] = Ei2 * ( Ci2*wciy*wciz - Ci4*wcix );
	  Ki[2][0] = Ei2 * ( Ci2*wciz*wcix - Ci4*wciy );
	  Ki[2][1] = Ei2 * ( Ci2*wciz*wciy + Ci4*wcix );
	  Ki[2][2] = Ei2 * ( Ci2*wciz*wciz + Ci3 );

	  Emidx = ( dx[i-1] * Ex[i] + dx[i] * Ex[i-1] ) / (dx[i] + dx[i-1]);
	 
	  // okay, updates for J finally.

	    Jex0 = Ae[0][0] * Jex[i] + Ae[0][1] * Jey[i] + Ae[0][2] * Jez[i] \
	      + E0 * wpe[i]*wpe[i] * ( Ke[0][0] * Emidx + Ke[0][1] * Ey[i] + Ke[0][2] * Ez[i] );

	    Jey0 = Ae[1][0] * Jex[i] + Ae[1][1] * Jey[i] + Ae[1][2] * Jez[i] \
	      + E0 * wpe[i]*wpe[i] * ( Ke[1][0] * Emidx + Ke[1][1] * Ey[i] + Ke[1][2] * Ez[i] );

	    Jez[i] = Ae[2][0] * Jex[i] + Ae[2][1] * Jey[i] + Ae[2][2] * Jez[i] \
	      + E0 * wpe[i]*wpe[i] * ( Ke[2][0] * Emidx + Ke[2][1] * Ey[i] + Ke[2][2] * Ez[i] );

	    Jix0 = Ai[0][0] * Jix[i] + Ai[0][1] * Jiy[i] + Ai[0][2] * Jiz[i] \
	      + E0 * wpi[i]*wpi[i] * ( Ki[0][0] * Emidx + Ki[0][1] * Ey[i] + Ki[0][2] * Ez[i] );
	  
	    Jiy0 = Ai[1][0] * Jix[i] + Ai[1][1] * Jiy[i] + Ai[1][2] * Jiz[i] \
	      + E0 * wpi[i]*wpi[i] * ( Ki[1][0] * Emidx + Ki[1][1] * Ey[i] + Ki[1][2] * Ez[i] );
	
	    Jiz[i] = Ai[2][0] * Jix[i] + Ai[2][1] * Jiy[i] + Ai[2][2] * Jiz[i] \
	      + E0 * wpi[i]*wpi[i] * ( Ki[2][0] * Emidx + Ki[2][1] * Ey[i] + Ki[2][2] * Ez[i] );

	    // these temporary values are necessary so that the fields aren't updating twice!
	    Jex[i] = Jex0;
	    Jey[i] = Jey0;
	    Jix[i] = Jix0;
	    Jiy[i] = Jiy0;
      }

      // add currents together

      for (int i = 1; i < xx-1; i++) {
	  Jx[i] = Jex[i] + Jix[i];
	  Jy[i] = Jey[i] + Jiy[i];
	  Jz[i] = Jez[i] + Jiz[i];
      }

      // do heating: on axis
      for (int i = 1; i < xx-1; i++) {
	  JdotE[i] = Jx[i] * (dx[i-1] * Ex[i] + dx[i] * Ex[i-1])/(dx[i]+dx[i-1]) + Jy[i] * Ey[i] + Jz[i] * Ez[i];
	  heat[i] = heat[i] + dt * JdotE[i];
      }

    } // if doionosphere

    //////////////////////////////////////////////////////////////
    // ----------------------------------------

    // Figure out run time

    if (t == 50) {
      double runtimemin = (omp_get_wtime() - tloopStart) / 60.0;
      double totaltime = runtimemin * tsteps / 50.0;

      logfile.open("log.txt",std::fstream::app);
      logfile << "t = 50, time taken = " << runtimemin << " minutes.\n";
      logfile << "You can expect the total simulation to take " << totaltime << " minutes.\n";
      logfile.close();
    }
    
    if (t % (int)(tsteps/numfiles) == 0) {
      partialtime = 100 * t/tsteps;   
      logfile.open("log.txt",std::fstream::app); 
      logfile << "t = " << t << ": " << partialtime << " % Done...  ";

      DATATYPE Emax = 0;
      for (int i = 0; i < xx; i++) {
	if (fabs(Ez[i]) > Emax) {
	  Emax = fabs(Ez[i]);
	}
      }
      logfile << "Maximum abs(Ez) is " << Emax << "\n";
      logfile.close();

      // write to files

      FILE * filePtr;
 
      filePtr = fopen("output_E.dat","ab");
      fwrite(Ex,sizeof(DATATYPE),xx,filePtr);
      fwrite(Ey,sizeof(DATATYPE),xx,filePtr);
      fwrite(Ez,sizeof(DATATYPE),xx,filePtr);
      fclose(filePtr);

      filePtr = fopen("output_J.dat","ab");
      fwrite(Jex,sizeof(DATATYPE),xx,filePtr);
      fwrite(Jey,sizeof(DATATYPE),xx,filePtr);
      fwrite(Jez,sizeof(DATATYPE),xx,filePtr);
      fwrite(Jix,sizeof(DATATYPE),xx,filePtr);      
      fwrite(Jiy,sizeof(DATATYPE),xx,filePtr);      
      fwrite(Jiz,sizeof(DATATYPE),xx,filePtr);
      fclose(filePtr);

      filePtr = fopen("output_H.dat","ab");
      fwrite(Hx,sizeof(DATATYPE),xx,filePtr);
      fwrite(Hy,sizeof(DATATYPE),xx,filePtr);
      fwrite(Hz,sizeof(DATATYPE),xx,filePtr);
      fclose(filePtr);

      filePtr = fopen("output_K.dat","ab");
      fwrite(Eeff,sizeof(DATATYPE),xx,filePtr);
      fwrite(Ek,sizeof(DATATYPE),xx,filePtr);
      fwrite(heat,sizeof(DATATYPE),xx,filePtr);
      fwrite(S,sizeof(DATATYPE),xx,filePtr);
      fclose(filePtr);
      
      filePtr = fopen("output_D.dat","ab");
      fwrite(ne,sizeof(DATATYPE),xx,filePtr);
      fwrite(nOm,sizeof(DATATYPE),xx,filePtr);
      fwrite(nue,sizeof(DATATYPE),xx,filePtr);
      fclose(filePtr);

      filePtr = fopen("output_O.dat","ab");
      fwrite(nN21P,sizeof(DATATYPE),xx,filePtr);
      fwrite(nN22P,sizeof(DATATYPE),xx,filePtr);
      fwrite(nN2P1N,sizeof(DATATYPE),xx,filePtr);
      fwrite(nN2PM,sizeof(DATATYPE),xx,filePtr);
      fwrite(nO2P1N,sizeof(DATATYPE),xx,filePtr);
      fwrite(nOred,sizeof(DATATYPE),xx,filePtr);
      fwrite(nOgrn,sizeof(DATATYPE),xx,filePtr);
      fclose(filePtr);

    }
    
    // end giant time loop
  }

  /////////////////////////////////////////////////////////////////////////////

  FILE * ProbeFile;
  ProbeFile = fopen("Probe.dat","wb");
  fwrite(&nprobes,sizeof(int),1,ProbeFile);
  fwrite(&probex,sizeof(int),nprobes,ProbeFile);
  fwrite(&Exprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Eyprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Ezprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Hxprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Hyprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Hzprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fclose(ProbeFile);

  logfile.open("log.txt",std::fstream::app);
  logfile << "All done!\n";

  double totalruntime = (omp_get_wtime() - tStart) / 60.0;

  logfile << "Total run time = " << totalruntime << " minutes.\n";
  logfile.close();

  // end program - do not type below this!

}

//////////////////////////////////////////////////////////

