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

// above MI is mass of NO+ ion

//////////////////////

// version 2.0.1: November 2, 2013: now includes 2D neutral density
// 2.0.2: April 2014: can't remember the latest changes. Sorry.
// fork: adding SIBC to 3D code. Also modifying planet settings to match 2d code. Also changing source like 2D, to read from h-t file

///////////////////////////////////////////////////////////////////////////////////

int main()
{

  // time the whole ting
  double tStart = omp_get_wtime();

  //------------------------------------------------------------------------//
  //-------READ INPUTS FROM FILE -------------------------------------------//
  
  FILE * inputFile;
  inputFile = fopen("inputs.dat","rb");

  double RE;
  int dopml_top;
  int dopml_wall;
  int doionosphere;
  int doioniz;
  int doelve;
  int dodetach;
  int dotransmitter;
  int savefields [6];
  int groundmethod;
  double maxalt;
  double stepalt;
  double dr0;
  double dr1;
  double dr2;
  int nground;
  double range;
  double drange;
  double dt;
  int tsteps;
  double sig;
  double sigm;
  double camdist;
  double camalt;
  int elvesteps;
  int sourcedirection;
  int numfiles;
  int planet;
  int decfactor;
  int nprobes;
  int dogwave;
  double gwavemag;
  double gwavemaxalt;
  double gwavekr, gwavekh, gwavekp;

  fread(&RE,sizeof(double),1,inputFile);
  fread(&dopml_top,sizeof(int),1,inputFile);
  fread(&dopml_wall,sizeof(int),1,inputFile);
  fread(&doionosphere,sizeof(int),1,inputFile);
  fread(&doioniz,sizeof(int),1,inputFile);
  fread(&doelve,sizeof(int),1,inputFile);
  fread(&dodetach,sizeof(int),1,inputFile);
  fread(&dotransmitter,sizeof(int),1,inputFile);
  fread(&savefields,sizeof(int),6,inputFile);
  fread(&groundmethod,sizeof(int),1,inputFile);
  fread(&maxalt,sizeof(double),1,inputFile);
  fread(&stepalt,sizeof(double),1,inputFile);
  fread(&dr0,sizeof(double),1,inputFile); 
  fread(&dr1,sizeof(double),1,inputFile);
  fread(&dr2,sizeof(double),1,inputFile);
  fread(&nground,sizeof(int),1,inputFile);
  fread(&range,sizeof(double),1,inputFile);  
  fread(&drange,sizeof(double),1,inputFile);
  fread(&dt,sizeof(double),1,inputFile);
  fread(&tsteps,sizeof(int),1,inputFile);
  fread(&sig,sizeof(double),1,inputFile);
  fread(&sigm,sizeof(double),1,inputFile);
  fread(&camdist,sizeof(double),1,inputFile);
  fread(&camalt,sizeof(double),1,inputFile);
  fread(&elvesteps,sizeof(int),1,inputFile);
  fread(&sourcedirection,sizeof(int),1,inputFile);
  fread(&numfiles,sizeof(int),1,inputFile);
  fread(&planet,sizeof(int),1,inputFile);
  fread(&decfactor,sizeof(int),1,inputFile);
  fread(&nprobes,sizeof(int),1,inputFile);
  int prober[nprobes], probet[nprobes], probep[nprobes];
  fread(&prober,sizeof(int),nprobes,inputFile);
  fread(&probet,sizeof(int),nprobes,inputFile);
  fread(&probep,sizeof(int),nprobes,inputFile);
  fread(&dogwave,sizeof(int),1,inputFile);
  fread(&gwavemag,sizeof(double),1,inputFile);
  fread(&gwavemaxalt,sizeof(double),1,inputFile);
  fread(&gwavekr,sizeof(double),1,inputFile);
  fread(&gwavekh,sizeof(double),1,inputFile);
  fread(&gwavekp,sizeof(double),1,inputFile);
  fclose(inputFile);

  //------------------------------------------------------------------------//
  //-------READ SOURCE------------------------------------------------------//

  FILE * sourceFile;
  sourceFile = fopen("source.dat","rb");

  int nsalts;
  int nstimes;
  int halfchlength;
  fread(&nsalts,sizeof(int),1,sourceFile);
  fread(&nstimes,sizeof(int),1,sourceFile);
  fread(&halfchlength,sizeof(int),1,sourceFile);
  double Isource [nstimes][nsalts];
  fread(&Isource,sizeof(double),nsalts*nstimes,sourceFile);
  fclose(sourceFile);


  //------------------------------------------------------------------------//
  //-------GRID PARAMETERS -------------------------------------------------//

  // set up radial array

  // DATATYPE ngr = (DATATYPE)nground;
  int rr = stepalt/dr1 + (maxalt - stepalt)/dr2 + 1 + nground;
  int stepi = stepalt/dr1;
  DATATYPE r [rr];
  DATATYPE dr [rr-1];
  r[0] = RE - nground*dr0;  
  for (int i = 1; i < nground+1; i++) {
    r[i] = r[i-1] + dr0;
    dr[i-1] = r[i] - r[i-1];
  }
  for (int i = nground+1; i < rr; i++) {
    if (r[i-1] < RE + stepalt) { 
      r[i] = r[i-1] + dr1;
    } else { 
      r[i] = r[i-1] + dr2; 
    }
    dr[i-1] = r[i] - r[i-1];
  }

  // set up theta array

  DATATYPE thmax = range / RE;
  DATATYPE dth = drange / RE;
  int hh = floor(2 * thmax / dth) + 1;
  DATATYPE th [hh];
  th[0] = PI/2 - thmax;  
  for (int i = 1; i < hh; i++) {
    th[i] = th[i-1] + dth;
  }

  // set up phi array

  DATATYPE phmax = thmax;
  DATATYPE dph = dth;
  int pp = floor(2 * phmax / dph) + 1;
  DATATYPE ph [pp];
  ph[0] = -phmax;
  for (int i = 1; i < pp; i++) {
    ph[i] = ph[i-1] + dph;
  }

  // decimation for output fields - NOT USED IN 3D.
  int drr = rr;
  int dhh = hh;
  int dpp = pp;
  if (decfactor > 1) {
    drr = floor(rr/decfactor);
    dhh = floor(hh/decfactor);
    dpp = floor(pp/decfactor);
  }

  //------------------------------------------------------------------------//
  //-------GROUND PARAMETERS -----------------------------------------------//

  // input groundmethod can be 0 (PEC), 1 (SIBC), or 2 (Real parameters)

  // read ground sigma and epsilon. In 3D I only allow a scalar value (homogeneous ground)
  double gsig;
  double geps;
    
  FILE * groundFile;
  groundFile = fopen("ground.dat","rb");
  fread(&gsig,sizeof(double),1,groundFile);
  fread(&geps,sizeof(double),1,groundFile);
  fclose(groundFile);

  //------------------------------------------------------------------------//
  //-------MULTIPLY FACTORS FOR FDTD ---------------------------------------//

  DATATYPE c1h = (2*U0 - sigm*dt)/(2*U0 + sigm*dt);
  DATATYPE c2h = 2*dt/(2*U0 + sigm*dt);

  DATATYPE c1e [rr];
  DATATYPE c2e [rr];
  
  for (int i = 0; i < rr; i++) {
    c1e[i] = (2*E0 - sig*dt)/(2*E0 + sig*dt);
    c2e[i] = 2*dt/(2*E0 + sig*dt);
  }

  if (groundmethod == 2) {
    for (int i = 0; i < nground; i++) {
      c1e[i] = (2*E0*geps - gsig*dt)/(2*E0*geps + gsig*dt);
      c2e[i] = 2*dt/(2*E0*geps + gsig*dt);
    }
  }

  //-----------------------------------------------------------------------//
  //----   Surface Impedance Boundary Condition (SIBC) stuff --------------//
  //----   This is taken from Oh and Schutt-Aine, IEEE TRANSACTIONS ON ANTENNAS AND PROPAGATION, VOL. 43, NO. 7, JULY 1995 

  DATATYPE sibc_a;
  DATATYPE sibc_eta;
  DATATYPE sibc_C [6] = { 1.22646e-8, 2.56716e-6, 1.51777e-4, 4.42437e-3, 6.98268e-2, 0.42473 };
  DATATYPE sibc_omega [6] = { 4.06981e-6, 1.84651e-4, 3.24245e-3, 3.42849e-2, 0.23606, 0.83083 };
  DATATYPE sibc_K;
  DATATYPE sibc_pi1 [6];
  DATATYPE sibc_pi2 [6];
  DATATYPE sibc_pi3 [6];

  DATATYPE sumAi;
  DATATYPE sibc_Ai [6];
  DATATYPE sibc_Aiold [6];

  // will also need "old" tangential H fields
  DATATYPE Hpold [rr][hh][pp];
  DATATYPE Htold [rr][hh][pp];

  sibc_a = gsig/(E0*geps);
  sibc_eta = sqrt(U0/(E0*geps));
  for (int m = 0; m < 6; m++) {
    sibc_K = sibc_a * sibc_omega[m] * dt;
    sibc_pi3[m] = exp(-sibc_K);
    sibc_pi1[m] = sibc_eta * (sibc_C[m]/sibc_omega[m]) * (1 + (exp(-sibc_K) - 1)/sibc_K);
    sibc_pi2[m] = sibc_eta * (sibc_C[m]/sibc_omega[m]) * (1/sibc_K - exp(-sibc_K) * (1 + 1/sibc_K));
  }

  memset(sibc_Ai,0,sizeof(sibc_Ai));
  memset(sibc_Aiold,0,sizeof(sibc_Aiold));
  memset(Hpold,0,sizeof(Hpold));
  memset(Htold,0,sizeof(Htold));

  
  //------------------------------------------------------------------------//
  //-------MAGNETIC FIELD --------------------------------------------------//

  // read 0D magnetic field (varying B doesn't buy much in 3D)
  double Br;
  double Bt;
  double Bp;
  FILE * BFile;
  BFile = fopen("B0.dat","rb");
  fread(&Br,sizeof(double),1,BFile);
  fread(&Bt,sizeof(double),1,BFile);
  fread(&Bp,sizeof(double),1,BFile);
  fclose(BFile);

  // compute gyrofrequencies, store in 2D arrays

  DATATYPE wcer;
  DATATYPE wcet;
  DATATYPE wcep;
  DATATYPE wce0;
  DATATYPE wcir;
  DATATYPE wcit;
  DATATYPE wcip;
  DATATYPE wci0;

  wcer = -QE * Br / ME;
  wcet = -QE * Bt / ME;
  wcep = -QE * Bp / ME;
  wce0 = sqrt(wcer*wcer + wcet*wcet + wcep*wcep);
  wcir = QE * Br / MI;
  wcit = QE * Bt / MI;
  wcip = QE * Bp / MI;
  wci0 = sqrt(wcir*wcir + wcit*wcit + wcip*wcip);


  //------------------------------------------------------------------------//
  //-------OUTPUT TO LOG FILE ----------------------------------------------//

  ofstream logfile;
  logfile.open("log.txt");

  logfile << "\n";
  logfile << "------------------------------------------------------\n";
  logfile << " ------- 3D EMP Simulation ---------------------------\n";
  logfile << "DoPML = " << dopml_top << "(top) " << dopml_wall << "(wall)\n";
  logfile << "DoIonosphere = " << doionosphere << "\n";
  logfile << "DoIoniz = " << doioniz << "\n";
  logfile << "Magnetic field magnitude is " << wce0*ME/QE*1e9 << " nT\n";
  logfile << "Magnetic field vector is " << Br*1e9 << ", " << Bt*1e9 << ", " << Bp*1e9 << "\n";
  logfile << "Grid is " << rr << " (r) by " << hh << " (theta) by " << pp << " (phi) cells\n";
  logfile << "Time step is " << (dt*1e9) << " ns\n";
  logfile << "Simulation will run " << tsteps << " time steps.\n";
  logfile << "Will spit out fields " << numfiles << " times during the simulation.\n";
  double gwavelam = 2*PI/sqrt(gwavekh*gwavekh + gwavekp*gwavekp);
  logfile << "GWave = " << dogwave << "; GW wavelength = " << gwavelam << "\n";

  // write space parameters to sferic.dat file

  FILE * sfericFile;

  sfericFile = fopen("sferic.dat","wb");

  fwrite(&tsteps,sizeof(int),1,sfericFile);
  fwrite(&rr,sizeof(int),1,sfericFile);
  fwrite(&hh,sizeof(int),1,sfericFile);
  fwrite(&pp,sizeof(int),1,sfericFile);
  fwrite(&numfiles,sizeof(int),1,sfericFile);
  fwrite(&dt,sizeof(DATATYPE),1,sfericFile);
  fwrite(r,sizeof(DATATYPE),rr,sfericFile);
  fwrite(th,sizeof(DATATYPE),hh,sfericFile);
  fwrite(ph,sizeof(DATATYPE),pp,sfericFile);
  fwrite(&decfactor,sizeof(int),1,sfericFile);
  fclose(sfericFile);

  //------------------------------------------------------------------------//
  //-------ELVE INTEGRATION SETUP ------------------------------------------//

  //  int elvewritesteps = 1000;    // this is the number of time steps at which photons will be written to the elve cube.

  int numpixels;
  int cameratype;
  FILE * cameraFile;
  int elvemaxr;

  // elve max r is index into r of say, 100 km, above which we do not want to contribute photons to elve (they may come from boundary, or numerical issue, for example)
  for (int i = 0; i < rr; i++) {
    if (r[i] < (RE + 100e3)) {
      elvemaxr = i;
    }
  }
  std::cout << "elvemaxr = " << elvemaxr << "\n";

  cameraFile = fopen("camera.dat","rb");
  fread(&cameratype,sizeof(int),1,cameraFile);
  fread(&numpixels,sizeof(int),1,cameraFile);
  DATATYPE elveaz [numpixels];
  DATATYPE elveel [numpixels];
  fread(&elveaz,sizeof(double),numpixels,cameraFile);
  fread(&elveel,sizeof(double),numpixels,cameraFile);
  fclose(cameraFile);

  //int elvesteps = 1000;         // this is the number of time steps in the actual elve cube
  DATATYPE eti = (camdist - range)/C;
  DATATYPE etf = (camdist + range)/C + tsteps*dt;
  DATATYPE elvedt = (etf - eti)/(elvesteps-1);

  //  int elvesize [2] = {128, 64};

  DATATYPE elveN21P [numpixels][elvesteps];
  DATATYPE elveN22P [numpixels][elvesteps];
  DATATYPE elveN2P1N [numpixels][elvesteps];
  DATATYPE elveN2PM [numpixels][elvesteps];
  DATATYPE elveO2P1N [numpixels][elvesteps];

  logfile << "Size of elve array = " << sizeof(elveN21P) << "\n";

  memset(elveN21P,0,sizeof(elveN21P));
  memset(elveN22P,0,sizeof(elveN22P));
  memset(elveN2P1N,0,sizeof(elveN22P));
  memset(elveN2PM,0,sizeof(elveN22P));
  memset(elveO2P1N,0,sizeof(elveN22P));

  // need arrays thind, phind, tdelay, and raylength to use for integration

  int thind [numpixels][rr];
  int phind [numpixels][rr];
  DATATYPE tdelay [numpixels][rr];
  DATATYPE raylength [numpixels][rr];

  DATATYPE d, az, el, th2, phi, thetap;
  DATATYPE rc = RE + camalt;
  DATATYPE camth = camdist/RE;
  int elvetimeind;

  for (int m = 0; m < numpixels; m++) {

    az = elveaz[m];
    el = elveel[m];

    for (int k = stepi; k < rr; k++) { 
      
      d = rc * cos(PI/2 + el) + sqrt( r[k]*r[k] + rc*rc*(pow(cos(PI/2+el),2) - 1) );
      tdelay[m][k] = d/C/dt - eti/dt;
      th2 = asin( d * sin(PI/2+el) / r[k] );
      phi = asin( sin(az) * sin(th2) );
      phind[m][k] = floor(phi/dph + (pp-1)/2);
      if (az == 0) {
	thetap = th2;
      } else {
	thetap = asin( tan(phi) / tan(az) );
      }	
      thind[m][k] = floor((thetap-camth)/dth + (hh-1)/2);
      raylength[m][k] = dr2 / sin(el+th2);
      
    }
  }


  //------------------------------------------------------------------------//
  //------- FIELD INITIALIZATION --------------------------------------------//

  // electric field vector
  DATATYPE Er [rr][hh][pp];
  DATATYPE Et [rr][hh][pp];
  DATATYPE Ep [rr][hh][pp];
  // magnetic field vector
  DATATYPE Hr [rr][hh][pp];
  DATATYPE Ht [rr][hh][pp];
  DATATYPE Hp [rr][hh][pp];
  // spatial-averaged fields
  DATATYPE Erm, Etm, Epm, Hrm, Htm, Hpm;

  // probe fields
  DATATYPE Erprobe [tsteps][nprobes];
  DATATYPE Etprobe [tsteps][nprobes];
  DATATYPE Epprobe [tsteps][nprobes];
  DATATYPE Hrprobe [tsteps][nprobes];
  DATATYPE Htprobe [tsteps][nprobes];
  DATATYPE Hpprobe [tsteps][nprobes];
  int ir, it, ip;  // probe point indices

  // electron current
  DATATYPE Jer [rr][hh][pp];
  DATATYPE Jet [rr][hh][pp];
  DATATYPE Jep [rr][hh][pp];
  DATATYPE Jer0 = 0;
  DATATYPE Jet0 = 0;

  DATATYPE heat [rr][hh][pp];

  // electron Temperature
  DATATYPE Te [rr][hh][pp];
  DATATYPE Te0 [rr][hh][pp];

  // Electric field Parallel and Perpendicular components, for time-dependent stuff
  DATATYPE Epar, Eperp2, HpmR, HpmL, HtmR, HtmL;
  DATATYPE Eeff [rr][hh][pp];
  DATATYPE Emag [rr][hh][pp];

  // slices for output to file
  DATATYPE Eslicep [rr][hh];
  DATATYPE Eslicet [rr][pp];
  DATATYPE Eslicer [hh][pp];
  DATATYPE Jslicep [rr][hh];
  DATATYPE Jslicet [rr][pp];
  DATATYPE Jslicer [hh][pp];
  DATATYPE Hslicep [rr][hh];
  DATATYPE Hslicet [rr][pp];
  DATATYPE Hslicer [hh][pp];

  DATATYPE Eeffslice [rr][hh];
  DATATYPE Ekslice [rr][hh];
  DATATYPE heatslice [rr][hh];
  DATATYPE Sslice [rr][hh];
  DATATYPE Teslice [rr][hh];
  DATATYPE neslice [rr][hh];
  DATATYPE nOmslice [rr][hh];
  DATATYPE nuslice [rr][hh];
  DATATYPE nN21Pslice [rr][hh];
  DATATYPE nN22Pslice [rr][hh];
  DATATYPE nN2P1Nslice [rr][hh];
  DATATYPE nN2PMslice [rr][hh];
  DATATYPE nO2P1Nslice [rr][hh];

  // initialize 

  memset(Er,0,sizeof(Er));
  memset(Et,0,sizeof(Et));
  memset(Ep,0,sizeof(Ep));
  memset(Hr,0,sizeof(Hr));
  memset(Ht,0,sizeof(Ht));
  memset(Hp,0,sizeof(Hp));
  memset(Jer,0,sizeof(Jer));
  memset(Jet,0,sizeof(Jet));
  memset(Jep,0,sizeof(Jep));
  memset(Eeff,0,sizeof(Eeff));
  memset(Emag,0,sizeof(Emag));
  memset(heat,0,sizeof(heat));

  memset(Eslicep,0,sizeof(Eslicep));
  memset(Eslicet,0,sizeof(Eslicet));
  memset(Eslicer,0,sizeof(Eslicer));
  memset(Jslicep,0,sizeof(Jslicep));
  memset(Jslicet,0,sizeof(Jslicet));
  memset(Jslicer,0,sizeof(Jslicer));
  memset(Hslicep,0,sizeof(Hslicep));
  memset(Hslicet,0,sizeof(Hslicet));
  memset(Hslicer,0,sizeof(Hslicer));

  memset(Eeffslice,0,sizeof(Eeffslice));
  memset(Ekslice,0,sizeof(Ekslice));
  memset(heatslice,0,sizeof(heatslice));
  memset(Sslice,0,sizeof(Sslice));
  memset(Teslice,0,sizeof(Teslice));
  memset(neslice,0,sizeof(neslice));
  memset(nOmslice,0,sizeof(nOmslice));
  memset(nuslice,0,sizeof(nuslice));
  memset(nN21Pslice,0,sizeof(nN21Pslice));
  memset(nN22Pslice,0,sizeof(nN22Pslice));
  memset(nN2P1Nslice,0,sizeof(nN2P1Nslice));
  memset(nN2PMslice,0,sizeof(nN2PMslice));
  memset(nO2P1Nslice,0,sizeof(nO2P1Nslice));
 

  // some stuff for the current source

  int nextra = 5;
  int st = (hh+1)/2;
  int sp = (pp+1)/2;
  DATATYPE ifl;
  DATATYPE jfl;
  DATATYPE kfl;

  //------------------------------------------------------------------------//
  //-------NUMBER DENSITIES ------------------------------------------------//

  // set up electron density and B field

  double ne1 [rr];               // 1D electron density
  DATATYPE ne [rr][hh][pp];      // 2D electron density
  DATATYPE wpe [rr][hh][pp];     // plasma frequency
  double nd [rr];                // neutral density
  DATATYPE mue;                  // mobility
  DATATYPE nue [rr][hh][pp];     // collision frequency
  DATATYPE Ek [rr][hh][pp];      // breakdown field

  FILE * neFile;
  neFile = fopen("ne.dat","rb");
  fread(&ne1,sizeof(double),rr,neFile);
  fclose(neFile);

  // read ambient temperature
  double Tamb [rr];
  FILE * Tefile;
  Tefile = fopen("etemp.dat","rb");
  fread(&Tamb,sizeof(double),rr,Tefile);
  fclose(Tefile);

  // extrapolate 3D arrays
  for (int i = 0; i < rr; i++) {
    for (int j = 0; j < hh; j++) {
      for (int k = 0; k < pp; k++) {
	ne[i][j][k] = ne1[i];
	wpe[i][j][k] = QE * sqrt(ne[i][j][k] / (ME * E0));
	Te[i][j][k] = Tamb[i];
	Te0[i][j][k] = Tamb[i];
      }    
    }
  }

  // temperature intermediates; private variables in openmp, except S and components

  DATATYPE S [rr][hh][pp];
  DATATYPE Sr, St, Sp;
  DATATYPE gg,f_N2,f_O2,Lelast_N2,Lelast_O2,Lrot_N2,Lrot_O2,Lvib_N2,Lvib_O2,Le;
  DATATYPE JdotE [rr][hh][pp];

  memset(S,0,sizeof(S));
  memset(JdotE,0,sizeof(JdotE));


  //------------------------------------------------------------------------//
  //-------OPTICAL EXCITATION RATES ----------------------------------------//

  // // read and compute rate arrays to be used for nonlinear stuff

  FILE * rateFile;
  rateFile = fopen("rates.dat","rb");

  // code uses lookup table / interpolation, these are interpolation params

  int topE, botE;
  int nume;
  fread(&nume,sizeof(int),1,rateFile);
  double efield [nume];
  fread(&efield,sizeof(double),nume,rateFile);

  double ioniz [rr*nume];
  double attach [rr*nume];
  double mobility [rr*nume];
  double Ored [rr*nume];
  double Ogrn [rr*nume];
  double N21p [rr*nume];
  double N22p [rr*nume];
  double N2p1n [rr*nume];
  double N2pM [rr*nume];
  double O2p1n [rr*nume];

  fread(&ioniz,sizeof(double),rr*nume,rateFile);
  fread(&attach,sizeof(double),rr*nume,rateFile);
  fread(&mobility,sizeof(double),rr*nume,rateFile);
  fread(&Ored,sizeof(double),rr*nume,rateFile);
  fread(&Ogrn,sizeof(double),rr*nume,rateFile);
  fread(&N21p,sizeof(double),rr*nume,rateFile);
  fread(&N22p,sizeof(double),rr*nume,rateFile);
  fread(&N2p1n,sizeof(double),rr*nume,rateFile);
  fread(&N2pM,sizeof(double),rr*nume,rateFile);
  fread(&O2p1n,sizeof(double),rr*nume,rateFile);

  fclose(rateFile);

  // reorder them into a useful 2D array

  double viArray [rr][nume];
  double vaArray [rr][nume];
  double muArray [rr][nume];
  double OrArray [rr][nume];
  double OgArray [rr][nume];
  double N21pArray [rr][nume];
  double N22pArray [rr][nume];
  double N2p1nArray [rr][nume];
  double N2pMArray [rr][nume];
  double O2p1nArray [rr][nume];

  int index;

  for (int i = 0; i < rr; i++) {
    for (int m = 0; m < nume; m++) {
      index = rr*m + i;
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
  fread(&nd,sizeof(double),rr,ndFile);
  fclose(ndFile);

  // create 3D nd and ne arrays in order to input gravity wave
  DATATYPE nd3 [rr][hh][pp];
  DATATYPE gwfac, gwamp;
  DATATYPE scaleheight = 7e3;
  for (int i = 0; i < rr; i++) {
    gwamp = gwavemag / exp(gwavemaxalt/(2*scaleheight)) * exp((r[i]-RE)/(2*scaleheight));
    if (gwamp > gwavemag) gwamp = gwavemag;
    if (gwamp < -gwavemag) gwamp = -gwavemag;
    for (int j = 0; j < hh; j++) {
      for (int k = 0; k < pp; k++) {
	gwfac = gwamp * cos(gwavekh*j*dth*RE + gwavekp*k*dph*RE + gwavekr*(r[i]-RE));
	if (dogwave) {
	  nd3[i][j][k] = nd[i] * ( 1 + gwfac );
	  ne[i][j][k] = ne[i][j][k] * ( 1 + gwfac );
	} else {
	  nd3[i][j][k] = nd[i];
	}
      }
    }
  }

  for (int i = 0; i < rr; i++) {
    for (int j = 0; j < hh; j++) {
      for (int k = 0; k < pp; k++) {
	if (planet == 1) { // venus
	  Ek[i][j][k] = 7.07e7 * nd3[i][j][k] / nd[nground];
	  mue = 0.0018 * nd[nground] / nd3[i][j][k];
	} else if (planet == 0) { // earth
	  Ek[i][j][k] = 3.0e6 * nd3[i][j][k] / nd[nground];
	  mue = muArray[i][0] / nd3[i][j][k];
	} else if (planet == 2) { // saturn
	  Ek[i][j][k] = 7.07e7 * nd3[i][j][k] / nd[nground];
	  mue = 0.0018 * nd[nground] / nd3[i][j][k];
	} else {
	  // assume earth
	  Ek[i][j][k] = 3.0e6 * nd3[i][j][k] / nd[nground];
	  mue = muArray[i][0] / nd3[i][j][k];
	}
	nue[i][j][k] = (QE / ME) / mue;
      }
    }
  }

  // if do transmitter, collision frequency depends on Temperature and density

  if (dotransmitter) {
    for (int i = 0; i < rr; i++) {
      for (int j = 0; j < hh; j++) {
	for (int k = 0; k < pp; k++) {
	  nue[i][j][k] = 1.6 * ( 2.33e-17*(0.78*nd3[i][j][k])*(1.0 - 1.25e-4*Te[i][j][k])*Te[i][j][k] + 1.82e-16*(0.21*nd3[i][j][k])*(1.0 + 3.6e-2*sqrt(Te[i][j][k]))*sqrt(Te[i][j][k]) );
	}
      }
    }
  }

  //------------------------------------------------------------------------//
  //------- OPTICAL FIELDS -------------------------------------------------//

  DATATYPE EoEk;
  DATATYPE vi, va, vd;
  DATATYPE nOm [rr][hh][pp];

  // optical arrays
  DATATYPE vN21P, vN22P, vN2P1N, vN2PM, vO2P1N;
  DATATYPE nN21P [rr][hh][pp];
  DATATYPE nN22P [rr][hh][pp];
  DATATYPE nN2P1N [rr][hh][pp];
  DATATYPE nN2PM [rr][hh][pp];
  DATATYPE nO2P1N [rr][hh][pp];
  DATATYPE tauN22P, tauN21P, tauN2P1N, tauN2PM, tauO2P1N;

  // initialize optical

  memset(nN21P,0,sizeof(nN21P));
  memset(nN22P,0,sizeof(nN22P));
  memset(nN2P1N,0,sizeof(nN2P1N));
  memset(nN2PM,0,sizeof(nN2PM));
  memset(nO2P1N,0,sizeof(nO2P1N));
  memset(nOm,0,sizeof(nOm));


  //------------------------------------------------------------------------//
  //------- A COUPLE OF DETAILS --------------------------------------------//

  // partialtime used to determine write times
  DATATYPE partialtime;

  // these are A and K matrices for Lee&Kalluri current solution
  DATATYPE Ae [3][3];
  DATATYPE Ke [3][3];
  DATATYPE Emid [3];

  // intermediate calculations in Lee&Kalluri
  DATATYPE Ee1, Ee2, Se1, Ce1, Ce2, Ce3, Ce4;

  //------------------------------------------------------------------------//
  //------- SET UP PML -----------------------------------------------------//

  // initialize A, B, and P arrays even if we don't use them

  double pmlm = 4;
  int pmllen = 10;

  DATATYPE A_ert [2*pmllen];
  DATATYPE A_erp [2*pmllen];
  DATATYPE A_etr [2*pmllen];
  DATATYPE A_etp [2*pmllen];
  DATATYPE A_epr [2*pmllen];
  DATATYPE A_ept [2*pmllen];
  DATATYPE A_hrt [2*pmllen];
  DATATYPE A_hrp [2*pmllen];
  DATATYPE A_htr [2*pmllen];
  DATATYPE A_htp [2*pmllen];
  DATATYPE A_hpr [2*pmllen];
  DATATYPE A_hpt [2*pmllen];

  DATATYPE B_ert [2*pmllen];
  DATATYPE B_erp [2*pmllen];
  DATATYPE B_etr [2*pmllen];
  DATATYPE B_etp [2*pmllen];
  DATATYPE B_epr [2*pmllen];
  DATATYPE B_ept [2*pmllen];
  DATATYPE B_hrt [2*pmllen];
  DATATYPE B_hrp [2*pmllen];
  DATATYPE B_htr [2*pmllen];
  DATATYPE B_htp [2*pmllen];
  DATATYPE B_hpr [2*pmllen];
  DATATYPE B_hpt [2*pmllen];

  DATATYPE P_ert [rr][2*pmllen][pp];
  DATATYPE P_erp [rr][hh][2*pmllen];
  DATATYPE P_etr [2*pmllen][hh][pp];
  DATATYPE P_etp [rr][hh][2*pmllen];
  DATATYPE P_epr [2*pmllen][hh][pp];
  DATATYPE P_ept [rr][2*pmllen][pp];
  
  DATATYPE P_hrt [rr][2*pmllen][pp];
  DATATYPE P_hrp [rr][hh][2*pmllen];
  DATATYPE P_htr [2*pmllen][hh][pp];
  DATATYPE P_htp [rr][hh][2*pmllen];
  DATATYPE P_hpr [2*pmllen][hh][pp];
  DATATYPE P_hpt [rr][2*pmllen][pp];

  // initialize
  memset(P_erp,0,sizeof(P_erp));
  memset(P_etp,0,sizeof(P_etp));
  memset(P_hrp,0,sizeof(P_hrp));
  memset(P_htp,0,sizeof(P_htp));
  memset(P_ert,0,sizeof(P_ert));
  memset(P_ept,0,sizeof(P_ept));
  memset(P_hrt,0,sizeof(P_hrt));
  memset(P_hpt,0,sizeof(P_hpt));
  memset(P_etr,0,sizeof(P_etr));
  memset(P_epr,0,sizeof(P_epr));
  memset(P_htr,0,sizeof(P_htr));
  memset(P_hpr,0,sizeof(P_hpr));

  // end initialization

  if (dopml_top || dopml_wall) {
    DATATYPE simaxr = 1.5 * (DATATYPE)(pmlm + 1) / (150 * PI * dr[rr-2]);
    DATATYPE simaxt = 1.5 * (DATATYPE)(pmlm + 1) / (150 * PI * drange);
    DATATYPE kamax = 1;
    DATATYPE almax = 0;

    DATATYPE st [pmllen];
    DATATYPE stm [pmllen];
    DATATYPE sr [pmllen];
    DATATYPE srm [pmllen];

    DATATYPE kt [pmllen];
    DATATYPE ktm [pmllen];
    DATATYPE kr [pmllen];
    DATATYPE krm [pmllen];

    DATATYPE at [pmllen];
    DATATYPE atm [pmllen];
    DATATYPE ar [pmllen];
    DATATYPE arm [pmllen];
    DATATYPE mf;

    for (int m = 0; m < pmllen; m++) {
      // note shifts here compared to matlab, due to zero index
      mf = (DATATYPE)m;
      st[m] = simaxt * pow(((mf+0.5)/pmllen),pmlm);
      stm[m] = simaxt * pow(((mf+1)/pmllen),pmlm);
      sr[m] = simaxr * pow(((mf+0.5)/pmllen),pmlm);
      srm[m] = simaxr * pow(((mf+1)/pmllen),pmlm);

      kt[m] = 1 + (kamax-1) * pow(((mf+0.5)/pmllen),pmlm);
      ktm[m] = 1 + (kamax-1) * pow(((mf+1)/pmllen),pmlm);
      kr[m] = 1 + (kamax-1) * pow(((mf+0.5)/pmllen),pmlm);
      krm[m] = 1 + (kamax-1) * pow(((mf+1)/pmllen),pmlm);

      at[m] = almax * pow(((pmllen-mf-0.5)/pmllen),pmlm);
      atm[m] = almax * pow(((pmllen-mf-1)/pmllen),pmlm);
      ar[m] = almax * pow(((pmllen-mf-0.5)/pmllen),pmlm);
      arm[m] = almax * pow(((pmllen-mf-1)/pmllen),pmlm);
    }

    // update / initialize A and B vectors
    for (int m = pmllen; m < 2*pmllen; m++) {
      int n = m - pmllen;
      B_ert[m] = exp(-((st[n]/kt[n]) + at[n])*dt/E0);
      A_ert[m] = st[n] / (st[n]*kt[n] + pow(kt[n],2)*at[n]) * (B_ert[m] - 1);
      B_etr[m] = exp(-((sr[n]/kr[n]) + ar[n])*dt/E0);
      A_etr[m] = sr[n] / (sr[n]*kr[n] + pow(kr[n],2)*ar[n]) * (B_etr[m] - 1);
      B_ept[m] = exp(-((st[n]/kt[n]) + at[n])*dt/E0);
      A_ept[m] = st[n] / (st[n]*kt[n] + pow(kt[n],2)*at[n]) * (B_ept[m] - 1);
      B_epr[m] = exp(-((sr[n]/kr[n]) + ar[n])*dt/E0);
      A_epr[m] = sr[n] / (sr[n]*kr[n] + pow(kr[n],2)*ar[n]) * (B_epr[m] - 1);
       
      B_hrt[m] = exp(-((stm[n]/ktm[n]) + atm[n])*dt/E0);
      A_hrt[m] = stm[n] / (stm[n]*ktm[n] + pow(ktm[n],2)*atm[n]) * (B_hrt[m] - 1);
      B_htr[m] = exp(-((srm[n]/krm[n]) + arm[n])*dt/E0);
      A_htr[m] = srm[n] / (srm[n]*krm[n] + pow(krm[n],2)*arm[n]) * (B_htr[m] - 1);
      B_hpt[m] = exp(-((stm[n]/ktm[n]) + atm[n])*dt/E0);
      A_hpt[m] = stm[n] / (stm[n]*ktm[n] + pow(ktm[n],2)*atm[n]) * (B_hpt[m] - 1);
      B_hpr[m] = exp(-((srm[n]/krm[n]) + arm[n])*dt/E0);
      A_hpr[m] = srm[n] / (srm[n]*krm[n] + pow(krm[n],2)*arm[n]) * (B_hpr[m] - 1);  
      
      // using theta conductivities for phi direction. 

      B_erp[m] = exp(-((st[n]/kt[n]) + at[n])*dt/E0);
      A_erp[m] = st[n] / (st[n]*kt[n] + pow(kt[n],2)*at[n]) * (B_erp[m] - 1);
      B_etp[m] = exp(-((st[n]/kt[n]) + at[n])*dt/E0);
      A_etp[m] = st[n] / (st[n]*kt[n] + pow(kt[n],2)*at[n]) * (B_etp[m] - 1);
      B_hrp[m] = exp(-((stm[n]/ktm[n]) + atm[n])*dt/E0);
      A_hrp[m] = stm[n] / (stm[n]*ktm[n] + pow(ktm[n],2)*atm[n]) * (B_hrp[m] - 1);
      B_htp[m] = exp(-((stm[n]/ktm[n]) + atm[n])*dt/E0);
      A_htp[m] = stm[n] / (stm[n]*ktm[n] + pow(ktm[n],2)*atm[n]) * (B_htp[m] - 1);
    }
    for (int m = 0; m < pmllen; m++) {
      B_ert[m] = B_ert[2*pmllen-1-m];
      A_ert[m] = A_ert[2*pmllen-1-m];
      B_etr[m] = B_etr[2*pmllen-1-m];
      A_etr[m] = A_etr[2*pmllen-1-m];
      B_ept[m] = B_ept[2*pmllen-1-m];
      A_ept[m] = A_ept[2*pmllen-1-m];
      B_epr[m] = B_epr[2*pmllen-1-m];
      A_epr[m] = A_epr[2*pmllen-1-m];

      B_hrt[m] = B_hrt[2*pmllen-1-m];
      A_hrt[m] = A_hrt[2*pmllen-1-m];
      B_htr[m] = B_htr[2*pmllen-1-m];
      A_htr[m] = A_htr[2*pmllen-1-m];
      B_hpt[m] = B_hpt[2*pmllen-1-m];
      A_hpt[m] = A_hpt[2*pmllen-1-m];
      B_hpr[m] = B_hpr[2*pmllen-1-m];
      A_hpr[m] = A_hpr[2*pmllen-1-m];

      B_erp[m] = B_erp[2*pmllen-1-m];
      A_erp[m] = A_erp[2*pmllen-1-m];
      B_etp[m] = B_etp[2*pmllen-1-m];
      A_etp[m] = A_etp[2*pmllen-1-m];
      B_hrp[m] = B_hrp[2*pmllen-1-m];
      A_hrp[m] = A_hrp[2*pmllen-1-m];
      B_htp[m] = B_htp[2*pmllen-1-m];
      A_htp[m] = A_htp[2*pmllen-1-m];
    }

    // PML: everything up to here seems gtg


  } // if dopml

  // offsets for PML
  int rshift = rr-2*pmllen-1;
  int hshift = hh-2*pmllen-1;
  int pshift = pp-2*pmllen-1;

  //-----------------------------------------------------------------------//
  // measure memory usage, just for kicks
  struct rusage ru;
  double memusage;
  getrusage(RUSAGE_SELF, &ru);
  memusage = (double)ru.ru_maxrss / 1024 / 1024;

  logfile << "Total memory usage = " << memusage << " GB\n";
  logfile.close();


  ////////////////////////////////////////////////////////////////////////////
  //                                                                      ////
  //------- MAIN TIME LOOP -----------------------------------------------////

  double tloopStart = omp_get_wtime();

  for (int t = 0; t < tsteps; t++) {

    // ----------------------------------------
    // Psi updates for H field
    if (dopml_wall) {

      // Psi's applied to Hr
      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < pmllen; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    P_hrt[i][j][k] = B_hrt[j] * P_hrt[i][j][k] \
	      + A_hrt[j] * ( sin(th[j+1]) * Ep[i][j+1][k] - sin(th[j]) * Ep[i][j][k] ) / dth;
	  }
	}
	for (int j = hh-pmllen-1; j < hh-1; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    P_hrt[i][j-hshift][k] = B_hrt[j-hshift] * P_hrt[i][j-hshift][k] \
	      + A_hrt[j-hshift] * ( sin(th[j+1]) * Ep[i][j+1][k] - sin(th[j]) * Ep[i][j][k] ) / dth;
	  }
	}
      }

      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < hh-1; j++) {
	  for (int k = 0; k < pmllen; k++) {
	    P_hrp[i][j][k] = B_hrp[k] * P_hrp[i][j][k] \
	      + A_hrp[k] * ( Et[i][j][k+1] - Et[i][j][k] ) / dph;
	  }
	  for (int k = pp-pmllen-1; k < pp-1; k++) {
	    P_hrp[i][j][k-pshift] = B_hrp[k-pshift] * P_hrp[i][j][k-pshift] \
	      + A_hrp[k-pshift] * ( Et[i][j][k+1] - Et[i][j][k] ) / dph;
	  }
	}
      }
    }

    // Psi's applied to Ht
    if (dopml_top) {
      for (int i = rr-pmllen-1; i < rr-1; i++) {
	for (int j = 0; j < hh; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    P_htr[i-rshift][j][k] = B_htr[i-rshift] * P_htr[i-rshift][j][k] \
	      + A_htr[i-rshift] * ( r[i+1] * Ep[i+1][j][k] - r[i] * Ep[i][j][k] ) / dr[i]; 
	  }
	}
      }
    }

    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < hh; j++) {
	  for (int k = 0; k < pmllen; k++) {
	    P_htp[i][j][k] = B_htp[k] * P_htp[i][j][k] \
	      + A_htp[k] * ( Er[i][j][k+1] - Er[i][j][k] ) / dph;
	  }
	  for (int k = pp-pmllen-1; k < pp-1; k++) {
	    P_htp[i][j][k-pshift] = B_htp[k-pshift] * P_htp[i][j][k-pshift] \
	      + A_htp[k-pshift] * ( Er[i][j][k+1] - Er[i][j][k] ) / dph;
	  }
	}
      }
    }

    // Psi's applied to Hp
    if (dopml_top) {
      for (int i = rr-pmllen-1; i < rr-1; i++) {
	for (int j = 0; j < hh-1; j++) {
	  for (int k = 0; k < pp; k++) {
	    P_hpr[i-rshift][j][k] = B_hpr[i-rshift] * P_hpr[i-rshift][j][k] \
	      + A_hpr[i-rshift] * ( r[i+1] * Et[i+1][j][k] - r[i] * Et[i][j][k] ) / dr[i];
	  }
	}
      }
    }

    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  for (int k = 0; k < pp; k++) {
	    P_hpt[i][j][k] = B_hpt[j] * P_hpt[i][j][k] \
	      + A_hpt[j] * ( Er[i][j+1][k] - Er[i][j][k] ) / dth;
	  }
	}
	for (int j = hh-pmllen-1; j < hh-1; j++) {
	  for (int k = 0; k < pp; k++) {
	    P_hpt[i][j-hshift][k] = B_hpt[j-hshift] * P_hpt[i][j-hshift][k] \
	      + A_hpt[j-hshift] * (Er[i][j+1][k] - Er[i][j][k] ) / dth;
	  }
	}
      }
    }


    // ----------------------------------------
    // Hr update

#pragma omp parallel for
    for (int i = 0; i < rr; i++ ) {
      for (int j = 0; j < hh-1; j++) {
	for (int k = 0; k < pp-1; k++) {
	  Hr[i][j][k] = c1h * Hr[i][j][k] - c2h / (r[i] * sin(th[j]+dth/2)) * \
	    ( ( sin(th[j+1]) * Ep[i][j+1][k] - sin(th[j]) * Ep[i][j][k] ) / dth - ( Et[i][j][k+1] - Et[i][j][k] ) / dph );
	}
      }
    }
    // spit it out to file
    for (int i = 0; i < rr; i++) {
      for (int j = 0; j < hh-1; j++) {
	Hslicep[i][j] = Hr[i][j][sp];
      }
    }
    for (int i = 0; i < rr; i++) {
      for (int k = 0; k < pp-1; k++) {
	Hslicet[i][k] = Hr[i][st][k];
      }
    }
    for (int j = 0; j < hh-1; j++) {
      for (int k = 0; k < pp-1; k++) {
	Hslicer[j][k] = Hr[stepi+10][j][k];
      }
    }

    // pml corrections
    if (dopml_wall) {
      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < pmllen; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    Hr[i][j][k] = Hr[i][j][k] - c2h / (r[i] * sin(th[j]+dth/2) ) * P_hrt[i][j][k];
	  }
	}
	for (int j = hh-pmllen-1; j < hh-1; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    Hr[i][j][k] = Hr[i][j][k] - c2h / (r[i] * sin(th[j]+dth/2) ) * P_hrt[i][j-hshift][k];
	  }
	}
      }

      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < hh-1; j++) {
	  for (int k = 0; k < pmllen; k++) {
	    Hr[i][j][k] = Hr[i][j][k] + c2h / (r[i] * sin(th[j]+dth/2) ) * P_hrp[i][j][k];
	  }
	  for (int k = pp-pmllen-1; k < pp-1; k++) {
	    Hr[i][j][k] = Hr[i][j][k] + c2h / (r[i] * sin(th[j]+dth/2) ) * P_hrp[i][j][k-pshift];
	  }
	}
      }

    } // if dopml
    

      // ----------------------------------------
      // Ht update
#pragma omp parallel for
    for (int i = 0; i < rr-1; i++) {
      for (int j = 0; j < hh; j++) {
	for (int k = 0; k < pp-1; k++) {
	  Ht[i][j][k] = c1h * Ht[i][j][k] - c2h / (r[i]+dr[i]/2) * \
	    ( ( Er[i][j][k+1] - Er[i][j][k] ) / (sin(th[j]) * dph ) - ( r[i+1] * Ep[i+1][j][k] - r[i] * Ep[i][j][k] ) / dr[i] );
	}
      }
    }

    // pml corrections
    if (dopml_top) {
      for (int i = rr-pmllen-1; i < rr-1; i++) {
	for (int j = 0; j < hh; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    Ht[i][j][k] = Ht[i][j][k] + c2h / (r[i]+dr[i]/2) * P_htr[i-rshift][j][k];
	  }
	}
      }
    }

    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < hh; j++) {
	  for (int k = 0; k < pmllen; k++) {
	    Ht[i][j][k] = Ht[i][j][k] - c2h / ( (r[i]+dr[i]/2) * sin(th[j]) ) * P_htp[i][j][k];
	  }
	  for (int k = pp-pmllen-1; k < pp-1; k++) {
	    Ht[i][j][k] = Ht[i][j][k] - c2h / ( (r[i]+dr[i]/2) * sin(th[j]) ) * P_htp[i][j][k-pshift];
	  }
	}
      }

    } // if dopml
    

      // ----------------------------------------
      // Hp update
#pragma omp parallel for
    for (int i = 0; i < rr-1; i++) {
      for (int j = 0; j < hh-1; j++) {
	for (int k = 0; k < pp; k++) {
	  Hp[i][j][k] = c1h * Hp[i][j][k] - c2h / (r[i]+dr[i]/2) * \
	    ( ( r[i+1] * Et[i+1][j][k] - r[i] * Et[i][j][k] ) / dr[i] - ( Er[i][j+1][k] - Er[i][j][k] ) / dth );
	}
      }
    }

    // pml corrections
    if (dopml_top) {
      for (int i = rr-pmllen-1; i < rr-1; i++) {
	for (int j = 0; j < hh-1; j++) {
	  for (int k = 0; k < pp; k++) {
	    Hp[i][j][k] = Hp[i][j][k] - c2h / (r[i]+dr[i]/2) * P_hpr[i-rshift][j][k];
	  }	
	}
      }
    }

    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  for (int k = 0; k < pp; k++) {
	    Hp[i][j][k] = Hp[i][j][k] + c2h / (r[i]+dr[i]/2) * P_hpt[i][j][k];
	  }
	}
	for (int j = hh-pmllen-1; j < hh-1; j++) {
	  for (int k = 0; k < pp; k++) {
	    Hp[i][j][k] = Hp[i][j][k] + c2h / (r[i]+dr[i]/2) * P_hpt[i][j-hshift][k];
	  }
	}
      }

    } // if dopml

      // ---------------------------------------
      // Psi updates for E field
    
    if (dopml_wall) {

      // Psi's applied to Er
      for (int i = 0; i < rr-1; i++) {
	for (int j = 1; j < pmllen+1; j++) {
	  for (int k = 1; k < pp-1; k++) {
	    P_ert[i][j-1][k] = B_ert[j-1] * P_ert[i][j-1][k] \
	      + A_ert[j-1] * ( sin(th[j]+dth/2) * Hp[i][j][k] - sin(th[j]-dth/2) * Hp[i][j-1][k] ) / dth;
	  }
	}
	for (int j = hh-pmllen-1; j < hh-1; j++) {
	  for (int k = 1; k < pp-1; k++) {
	    P_ert[i][j-hshift][k] = B_ert[j-hshift] * P_ert[i][j-hshift][k] \
	      + A_ert[j-hshift] * ( sin(th[j]+dth/2) * Hp[i][j][k] - sin(th[j]-dth/2) * Hp[i][j-1][k] ) / dth;
	  }
	}
      }

      for (int i = 0; i < rr; i++) {
	for (int j = 1; j < hh-1; j++) {
	  for (int k = 1; k < pmllen+1; k++) {
	    P_erp[i][j][k-1] = B_erp[k-1] * P_erp[i][j][k-1] \
	      + A_erp[k-1] * ( Ht[i][j][k] - Ht[i][j][k-1] ) / dph;
	  }
	  for (int k = pp-pmllen-1; k < pp-1; k++) {
	    P_erp[i][j][k-pshift] = B_erp[k-pshift] * P_erp[i][j][k-pshift] \
	      + A_erp[k-pshift] * ( Ht[i][j][k] - Ht[i][j][k-1] ) / dph;
	  }
	}
      }
    }

    // Psi's applied to Et
    if (dopml_top) {
      for (int i = rr-pmllen-1; i < rr-1; i++) {
	for (int j = 0; j < hh-1; j++) {
	  for (int k = 1; k < pp-1; k++) {
	    P_etr[i-rshift][j][k] = B_etr[i-rshift] * P_etr[i-rshift][j][k] \
	      + A_etr[i-rshift] * ( (r[i]+dr[i]/2) * Hp[i][j][k] - (r[i]-dr[i-1]/2) * Hp[i-1][j][k] ) / ((dr[i]+dr[i-1])/2); 
	  }
	}
      }
    }

    if (dopml_wall) {
      for (int i = 1; i < rr-1; i++) {
	for (int j = 0; j < hh-1; j++) {
	  for (int k = 1; k < pmllen+1; k++) {
	    P_etp[i][j][k-1] = B_etp[k-1] * P_etp[i][j][k-1] \
	      + A_etp[k-1] * ( Hr[i][j][k] - Hr[i][j][k-1] ) / dph;
	  }
	  for (int k = pp-pmllen-1; k < pp-1; k++) {
	    P_etp[i][j][k-pshift] = B_etp[k-pshift] * P_etp[i][j][k-pshift] \
	      + A_etp[k-pshift] * ( Hr[i][j][k] - Hr[i][j][k-1] ) / dph;
	  }
	}
      }
    }

    // Psi's applied to Ep
    if (dopml_top) {
      for (int i = rr-pmllen-1; i < rr-1; i++) {
	for (int j = 1; j < hh-1; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    P_epr[i-rshift][j][k] = B_epr[i-rshift] * P_epr[i-rshift][j][k] \
	      + A_epr[i-rshift] * ( (r[i]+dr[i]/2) * Ht[i][j][k] - (r[i]-dr[i-1]/2) * Ht[i-1][j][k] ) / ((dr[i]+dr[i-1])/2);
	  }
	}
      }
    }

    if (dopml_wall) {
      for (int i = 1; i < rr-1; i++) {
	for (int j = 1; j < pmllen+1; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    P_ept[i][j-1][k] = B_ept[j-1] * P_ept[i][j-1][k] \
	      + A_ept[j-1] * ( Hr[i][j][k] - Hr[i][j-1][k] ) / dth;
	  }
	}
	for (int j = hh-pmllen-1; j < hh-1; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    P_ept[i][j-hshift][k] = B_ept[j-hshift] * P_ept[i][j-hshift][k] \
	      + A_ept[j-hshift] * (Hr[i][j][k] - Hr[i][j-1][k] ) / dth;
	  }
	}
      }
    } // if dopml
    

      // ----------------------------------------
      // Er update
#pragma omp parallel for
    for (int i = 0; i < rr-1; i++) {
      for (int j = 1; j < hh-1; j++) {
	for (int k = 1; k < pp-1; k++) {
	  Er[i][j][k] = c1e[i] * Er[i][j][k] + c2e[i] / ( (r[i]+dr[i]/2) * sin(th[j]) ) * \
	    ( ( sin(th[j]+dth/2) * Hp[i][j][k] - sin(th[j]-dth/2) * Hp[i][j-1][k] ) / dth - ( Ht[i][j][k] - Ht[i][j][k-1] ) / dph ) \
	    - c2e[i] * ( Jer[i+1][j][k] + Jer[i][j][k] ) / 2;
	}
      }
    }
    // spit it out to file
    for (int i = 0; i < rr; i++) {
      for (int j = 0; j < hh; j++) {
	Eslicep[i][j] = Er[i][j][sp];
      }
    }
    for (int i = 0; i < rr; i++) {
      for (int k = 0; k < pp; k++) {
	Eslicet[i][k] = Er[i][st][k];
      }
    }
    for (int j = 0; j < hh; j++) {
      for (int k = 0; k < pp; k++) {
	Eslicer[j][k] = Er[stepi+10][j][k];
      }
    }



    // source =======================================

    if (t < nstimes) {

      // vertical discharge
      if (sourcedirection == 0) {
	for (int i = 0; i < nsalts; i++) {
	  for (int j = st-nextra; j < st+nextra; j++) {
	    jfl = (DATATYPE)(j-st);
	    for (int k = sp-nextra; k < sp+nextra; k++) {
	      kfl = (DATATYPE)(k-sp);
	      Er[i][j][k] = Er[i][j][k] - c2e[i] * Isource[t][i] * exp(-jfl*jfl/9.0) * exp(-kfl*kfl/9.0) / (PI * 9.0*drange*drange);
	    }
	  }
	}
	// horizontal discharge. Altitude index is nsalts!
      } else if (sourcedirection == 1) {
	for (int i = nsalts-nextra; i < nsalts+nextra; i++) {
	  ifl = (DATATYPE)(i-nsalts);
	  for (int j = st-halfchlength; j < st+halfchlength; j++) {
	    for (int k = sp-nextra; k < sp+nextra; k++) {
	      kfl = (DATATYPE)(k-sp);
	      Er[i][j][k] = Er[i][j][k] - c2e[i] * Isource[t][j] * exp(-ifl*ifl/9.0) * exp(-kfl*kfl/9.0) / (PI * 9.0*drange*drange);

	    }
	  }
	}
      }

    }

    // ==============================================

    // pml corrections
    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 1; j < pmllen+1; j++) {
	  for (int k = 1; k < pp-1; k++) {
	    Er[i][j][k] = Er[i][j][k] + c2e[i] / ( (r[i]+dr[i]/2) * sin(th[j]) ) * P_ert[i][j-1][k];
	  }
	}
	for (int j = hh-pmllen-1; j < hh-1; j++) {
	  for (int k = 1; k < pp-1; k++) {
	    Er[i][j][k] = Er[i][j][k] + c2e[i] / ( (r[i]+dr[i]/2) * sin(th[j]) ) * P_ert[i][j-hshift][k];
	  }
	}
      }

      for (int i = 0; i < rr-1; i++) {
	for (int j = 1; j < hh-1; j++) {
	  for (int k = 1; k < pmllen+1; k++) {
	    Er[i][j][k] = Er[i][j][k] - c2e[i] / ( (r[i]+dr[i]/2) * sin(th[j]) ) * P_erp[i][j][k-1];
	  }
	  for (int k = pp-pmllen-1; k < pp-1; k++) {
	    Er[i][j][k] = Er[i][j][k] - c2e[i] / ( (r[i]+dr[i]/2) * sin(th[j]) ) * P_erp[i][j][k-pshift];
	  }
	}
      }

    } // if dopml
   

      // ----------------------------------------
      // Et update
#pragma omp parallel for
    for (int i = 1; i < rr-1; i++) {
      for (int j = 0; j < hh-1; j++) {
	for (int k = 1; k < pp-1; k++) {
	  Et[i][j][k] = c1e[i] * Et[i][j][k] + c2e[i] / r[i] *	\
	    ( ( Hr[i][j][k] - Hr[i][j][k-1] ) / ( sin(th[j]+dth/2) * dph ) - \
	      ( (r[i]+dr[i]/2) * Hp[i][j][k] - (r[i]-dr[i-1]/2) * Hp[i-1][j][k] ) / ((dr[i]+dr[i-1])/2) ) \
	    - c2e[i] * ( Jet[i][j+1][k] + Jet[i][j][k] ) / 2;
	  if (groundmethod == 0) Et[nground-1][j][k] = 0;
	}
      }
    }

    // SIBC boundary for Et
    if (groundmethod == 1) {
      for (int j = 0; j < hh-1; j++) {
	for (int k = 0; k < pp-1; k++) {
	  sumAi = 0;
	  for (int m = 0; m < 6; m++) {
	    sibc_Ai[m] = -sibc_pi1[m] * Hp[nground][j][k] - sibc_pi2[m] * Hpold[nground][j][k] + sibc_pi3[m] * sibc_Ai[m];
	    sumAi = sumAi + sibc_Ai[m];
	  }
	  Et[nground][j][k] = -sibc_eta * Hp[nground][j][k] - sumAi;
	}
      }
    }

    // pml corrections
    if (dopml_top) {
      for (int i = rr-pmllen-1; i < rr-1; i++) {
	for (int j = 0; j < hh-1; j++) {
	  for (int k = 1; k < pp-1; k++) {
	    Et[i][j][k] = Et[i][j][k] - c2e[i] / r[i] * P_etr[i-rshift][j][k];
	  }
	}
      }
    }

    if (dopml_wall) {
      for (int i = 1; i < rr-1; i++) {
	for (int j = 0; j < hh-1; j++) {
	  for (int k = 1; k < pmllen+1; k++) {
	    Et[i][j][k] = Et[i][j][k] + c2e[i] / ( r[i] * sin(th[j]+dth/2) ) * P_etp[i][j][k-1];
	  }
	  for (int k = pp-pmllen-1; k < pp-1; k++) {
	    Et[i][j][k] = Et[i][j][k] + c2e[i] / ( r[i] * sin(th[j]+dth/2) ) * P_etp[i][j][k-pshift];
	  }
	}
      }
    } // if dopml

      // ----------------------------------------
      // Ep update
#pragma omp parallel for
    for (int i = 1; i < rr-1; i++) {
      for (int j = 1; j < hh-1; j++) {
	for (int k = 0; k < pp-1; k++) {
	  Ep[i][j][k] = c1e[i] * Ep[i][j][k] + c2e[i] / r[i] *		\
	    ( ( (r[i]+dr[i]/2) * Ht[i][j][k] - (r[i]-dr[i-1]/2) * Ht[i-1][j][k] ) / ((dr[i]+dr[i-1])/2) - \
	      ( Hr[i][j][k] - Hr[i][j-1][k] ) / dth ) - c2e[i] * ( Jep[i][j][k+1] + Jep[i][j][k] ) / 2;
	  if (groundmethod == 0) Ep[nground-1][j][k] = 0;
	}
      }
    }

    // pml corrections
    if (dopml_top) {
      for (int i = rr-pmllen-1; i < rr-1; i++) {
	for (int j = 1; j < hh-1; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    Ep[i][j][k] = Ep[i][j][k] + c2e[i] / r[i] * P_epr[i-rshift][j][k];
	  }
	}
      }
    }

    if (dopml_wall) {
      for (int i = 1; i < rr-1; i++) {
	for (int j = 1; j < pmllen+1; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    Ep[i][j][k] = Ep[i][j][k] - c2e[i] / r[i] * P_ept[i][j-1][k];
	  }
	}
	for (int j = hh-pmllen-1; j < pp-1; j++) {
	  for (int k = 0; k < pp-1; k++) {
	    Ep[i][j][k] = Ep[i][j][k] - c2e[i] / r[i] * P_ept[i][j-hshift][k];
	  }
	}
      }
    } // if dopml

    // output E and H probe points at desired points

    for (int i = 0; i < nprobes; i++) {
      ir = prober[i];
      it = probet[i];
      ip = probep[i];
      Erprobe[t][i] = 0.5*(Er[ir][it][ip]+Er[ir-1][it][ip]);
      Etprobe[t][i] = 0.5*(Et[ir][it][ip]+Et[ir][it-1][ip]);
      Epprobe[t][i] = 0.5*(Ep[ir][it][ip]+Ep[ir][it][ip-1]);

      Hrprobe[t][i] = 0.25*(Hr[ir][it][ip]+Hr[ir][it-1][ip]+Hr[ir][it][ip-1]+Hr[ir][it-1][ip-1]);
      Htprobe[t][i] = 0.25*(Ht[ir][it][ip]+Ht[ir-1][it][ip]+Ht[ir][it][ip-1]+Ht[ir-1][it][ip-1]);
      Hpprobe[t][i] = 0.25*(Hp[ir][it][ip]+Hp[ir][it-1][ip]+Hp[ir-1][it][ip]+Hp[ir-1][it-1][ip]);

    }
    
    // ---------------------------------------
    // update E_eff, for use in ionization updates

    if (doionosphere) {
#pragma omp parallel for private(HtmL,HtmR,HpmL,HpmR,Epar,Eperp2,Sr,St,Sp,Erm,Etm,Epm,Hrm,Htm,Hpm)
      for (int i = nground+1; i < rr-1; i++) {
	for (int j = 1; j < hh-1; j++) {
	  for (int k = 1; k < pp-1; k++) {
	    
	    Erm = ( dr[i-1] * Er[i][j][k] + dr[i] * Er[i-1][j][k] ) / (dr[i] + dr[i-1]);
	    Etm = (Et[i][j][k] + Et[i][j-1][k])/2;
	    Epm = (Ep[i][j][k] + Ep[i][j][k-1])/2;
	    Hrm = (Hr[i][j][k] + Hr[i][j-1][k] + Hr[i][j][k-1] + Hr[i][j-1][k-1])/4; // need to think about field placements for these! what is at k = 0?
	    HtmL = ( dr[i-1] * Ht[i][j][k] + dr[i] * Ht[i-1][j][k] ) / (dr[i] + dr[i-1]);
	    HtmR = ( dr[i-1] * Ht[i][j][k-1] + dr[i] * Ht[i-1][j][k-1] ) / (dr[i] + dr[i-1]);
	    HpmL = ( dr[i-1] * Hp[i][j][k] + dr[i] * Hp[i-1][j][k] ) / (dr[i] + dr[i-1]);
	    HpmR = ( dr[i-1] * Hp[i][j-1][k] + dr[i] * Hp[i-1][j-1][k] ) / (dr[i] + dr[i-1]);
	    Htm = (HtmL + HtmR)/2;
	    Hpm = (HpmL + HpmR)/2;
	    Epar = (wcer/wce0)*Erm + (wcet/wce0)*Etm + (wcep/wce0)*Epm;
	    Emag[i][j][k] = sqrt( Erm*Erm + Etm*Etm + Epm*Epm );
	    Eperp2 = Emag[i][j][k]*Emag[i][j][k] - Epar*Epar;
	    Eeff[i][j][k] = sqrt( Epar*Epar + Eperp2 * nue[i][j][k]*nue[i][j][k] / ( nue[i][j][k]*nue[i][j][k] + wce0*wce0 ) );
	    
	    Sr = Etm*Hpm - Epm*Htm;
	    St = Epm*Hrm - Erm*Hpm;
	    Sp = Erm*Htm - Etm*Hrm;
	    S[i][j][k] = sqrt( Sr*Sr + St*St + Sp*Sp );
	    
	  }
	}
      }
      
      // test to find numerical instability!
      if (0) {
	for (int i = nground+1; i < rr-1; i++) {
	  for (int j = 0; j < hh-1; j++) {
	    for (int k = 0; k < pp-1; k++) {
	      if (Eeff[i][j][k] > 100*Ek[i][j][k]) {
		logfile.open("log.txt",std::fstream::app);	      
		logfile << "Unstable! Ending main loop after interation " << t << "...\n";
		logfile.close();
		goto endofbigloop;
	      }
	    }
	  }
	}
      }
      
      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < hh; j++) {
	  Eeffslice[i][j] = Eeff[i][j][sp];
	  Ekslice[i][j] = Ek[i][j][sp];
	  Sslice[i][j] = S[i][j][sp];
	}
      }
      
    } // if doionosphere
    
    // Ionization updates.
    // ------------------------------------------
    // update mue and nu, calculate vi, va, ne

    if (doionosphere & doioniz & !dotransmitter) {
      //if (0) {

#pragma omp parallel for private(EoEk,vN22P,vN21P,vN2P1N,vN2PM,vO2P1N,tauN22P,tauN21P,tauN2P1N,tauN2PM,tauO2P1N,Aion,nenew,botE,topE,vi,va,vd,mue)
      for (int i = nground+1; i < rr; i++) {
	for (int j = 0; j < hh; j++) {
	  for (int k = 0; k < pp; k++) {
	    
	    // earth
	    topE = 1;
	    botE = 0;
	    // find index into Efield rate array
	    for (int m = 0; m < nume; m++) {
	      if (efield[m] > Eeff[i][j][k] / nd3[i][j][k]) {
		topE = m;
		botE = topE - 1;
		break;
	      }
	    }
	      
	    // read off mue, vi, va, vN21P, vN22P from the arrays and interpolate
	  
	    if (topE <= 1) {
	      // set everything to first value 
	      mue = muArray[i][0] / nd3[i][j][k];
	      vi = 0;
	      va = 0;
	      vN21P = 0;
	      vN22P = 0;
	      vN2P1N = 0;
	      vN2PM = 0;
	      vO2P1N = 0;
	    } else {
	      //  interpolate
	      mue = (1/nd3[i][j][k]) * ( muArray[i][botE] + (Eeff[i][j][k]/nd3[i][j][k] - efield[botE])*((muArray[i][topE]-muArray[i][botE])/(efield[topE]-efield[botE])) );
	      vi = nd3[i][j][k] * ( viArray[i][botE] + (Eeff[i][j][k]/nd3[i][j][k] - efield[botE])*((viArray[i][topE]-viArray[i][botE])/(efield[topE]-efield[botE])) );
	      va = nd3[i][j][k] * ( vaArray[i][botE] + (Eeff[i][j][k]/nd3[i][j][k] - efield[botE])*((vaArray[i][topE]-vaArray[i][botE])/(efield[topE]-efield[botE])) );
	      vN21P = nd3[i][j][k] * ( N21pArray[i][botE] + (Eeff[i][j][k]/nd3[i][j][k] - efield[botE])*((N21pArray[i][topE]-N21pArray[i][botE])/(efield[topE]-efield[botE])) );
	      vN22P = nd3[i][j][k] * ( N22pArray[i][botE] + (Eeff[i][j][k]/nd3[i][j][k] - efield[botE])*((N22pArray[i][topE]-N22pArray[i][botE])/(efield[topE]-efield[botE])) );
	      vN2P1N = nd3[i][j][k] * ( N2p1nArray[i][botE] + (Eeff[i][j][k]/nd3[i][j][k] - efield[botE])*((N2p1nArray[i][topE]-N2p1nArray[i][botE])/(efield[topE]-efield[botE])) );
	      vN2PM = nd3[i][j][k] * ( N2pMArray[i][botE] + (Eeff[i][j][k]/nd3[i][j][k] - efield[botE])*((N2pMArray[i][topE]-N2pMArray[i][botE])/(efield[topE]-efield[botE])) );
	      vO2P1N = nd3[i][j][k] * ( O2p1nArray[i][botE] + (Eeff[i][j][k]/nd3[i][j][k] - efield[botE])*((O2p1nArray[i][topE]-O2p1nArray[i][botE])/(efield[topE]-efield[botE])) );
	    }
	      
	    // old method for detachment coefficient
	    EoEk = Eeff[i][j][k] / Ek[i][j][k];
	    if (EoEk < 0.05) {
	      vd = 0;
	    } else {
	      vd = 0.78 * nd3[i][j][k] * 1.08e-18 * sqrt(EoEk) * exp(-0.078375 / EoEk);
	    }
	      
	    // update electron density and plasma frequency, and excited state density
	    // test to make sure vi is not too big - things will go unstable! this happens inside right PML.
	      
	    if ( (vi - va + vd)*dt < 0.5 ) {
	      if ( dodetach ) {
		Aion[4] = 1 - (vi - va - vd)*dt - vi*vd*dt*dt; // determinant of solution matrix
		Aion[0] = (1 + vd*dt);
		Aion[1] = vd*dt;
		Aion[2] = va*dt;
		Aion[3] = 1 - (vi - va)*dt;
		nenew = ( Aion[0] * ne[i][j][k] + Aion[1] * nOm[i][j][k] ) / Aion[4];     // temporary value so we don't overwrite
		nOm[i][j][k] = ( Aion[2] * ne[i][j][k] + Aion[3] * nOm[i][j][k] ) / Aion[4];
		ne[i][j][k] = nenew;
	      } else {
		ne[i][j][k] = ne[i][j][k] / (1 - (vi - va)*dt);
		nOm[i][j][k] = nOm[i][j][k] + va*dt*ne[i][j][k];
	      }
	    } else {
	      // throttle the ionization so it doesn't go crazy unstable
	      ne[i][j][k] = ne[i][j][k] * exp(0.5);
	    }
	      
	    wpe[i][j][k] = QE * sqrt( ne[i][j][k] / (ME * E0) );
	    nue[i][j][k]  = (QE / ME) / mue;

	    // assume that ionization creates a positive ion, and attachment creates a negative ion?
	    // for now, we will assume ion densities don't change.

	    // optical: Need to input relative densities of O2, N2, O versus altitude to make these equations valid at higher altitudes.

	    tauN22P = 1/(2e7 + 3e-16*(0.21*nd3[i][j][k]) );  // quenching through O2, which is 21% of density
	    tauN21P = 1/(1.7e5 + 1e-17*(0.78*nd3[i][j][k]) );
	    tauN2P1N = 1/(1.4e7 + 4e-16*(0.99*nd3[i][j][k]) );
	    tauN2PM = 1/(7e4 + 5e-16*(0.78*nd3[i][j][k]) );
	    tauO2P1N = 1/(8.3e5 + 2e-16*(0.78*nd3[i][j][k]) );

	    // attemping backward Euler (implicit) for optics, since time constants are smaller than dt in some cases.

	    nN22P[i][j][k] = tauN22P/(dt + tauN22P) * nN22P[i][j][k] + (tauN22P*dt)/(dt + tauN22P) * ( vN22P * ne[i][j][k] );
	    nN21P[i][j][k] = tauN21P/(dt + tauN21P) * nN21P[i][j][k] + (tauN21P*dt)/(dt + tauN21P) * ( vN21P * ne[i][j][k] + 2e7 * nN22P[i][j][k] );
	    nN2P1N[i][j][k] = tauN2P1N/(dt + tauN2P1N) * nN2P1N[i][j][k] + (tauN2P1N*dt)/(dt + tauN2P1N) * ( vN2P1N * ne[i][j][k] );
	    nN2PM[i][j][k] = tauN2PM/(dt + tauN2PM) * nN2PM[i][j][k] + (tauN2PM*dt)/(dt + tauN2PM) * ( vN2PM * ne[i][j][k] );
	    nO2P1N[i][j][k] = tauO2P1N/(dt + tauO2P1N) * nO2P1N[i][j][k] + (tauO2P1N*dt)/(dt + tauO2P1N) * ( vO2P1N * ne[i][j][k] );

	  }
	}
      }

      // slice of nu, ne, and optics
      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < hh; j++) {
	  nuslice[i][j] = nue[i][j][sp];
	  neslice[i][j] = ne[i][j][sp];
	  nOmslice[i][j] = nOm[i][j][sp];
	  nN21Pslice[i][j] = nN21P[i][j][sp];
	  nN22Pslice[i][j] = nN22P[i][j][sp];
	  nN2P1Nslice[i][j] = nN2P1N[i][j][sp];
	  nN2PMslice[i][j] = nN2PM[i][j][sp];
	  nO2P1Nslice[i][j] = nO2P1N[i][j][sp];
	}
      }
      
    } // if doioniz

    // ------------------------------------------
    // optical integration to camera view

    if (doelve) {
#pragma omp parallel for private(elvetimeind) // parallel might be a problem here!
      for (int m = 0; m < numpixels; m++) {
	for (int k = stepi; k < elvemaxr; k++) {
	  if ((thind[m][k] > pmllen) && (thind[m][k] <= hh-pmllen) && (phind[m][k] > pmllen) && (phind[m][k] <= pp-pmllen)) {
	    
	    elvetimeind = round((tdelay[m][k] + t)*dt/elvedt);
	    if ((elvetimeind >= 0) && (elvetimeind < elvesteps)) {
	      // the constant below (e.g. 1.7e-5) is the Einstein coefficient (1.7e5) * 1e-6 (for Rayleigh definition) * 1e-4 (convert from m^2 to cm^2)
	      elveN21P[m][elvetimeind] = elveN21P[m][elvetimeind] + 1.7e-5 * nN21P[k][thind[m][k]][phind[m][k]] * raylength[m][k] * dt;
	      elveN22P[m][elvetimeind] = elveN22P[m][elvetimeind] + 2.0e-3 * nN22P[k][thind[m][k]][phind[m][k]] * raylength[m][k] * dt;
	      elveN2P1N[m][elvetimeind] = elveN2P1N[m][elvetimeind] + 1.4e-3 * nN2P1N[k][thind[m][k]][phind[m][k]] * raylength[m][k] * dt;
	      elveN2PM[m][elvetimeind] = elveN2PM[m][elvetimeind] + 7.1e-6 * nN2PM[k][thind[m][k]][phind[m][k]] * raylength[m][k] * dt;
	      elveO2P1N[m][elvetimeind] = elveO2P1N[m][elvetimeind] + 8.3e-5 * nO2P1N[k][thind[m][k]][phind[m][k]] * raylength[m][k] * dt;
	      // since we multiply by the effective dt, what we are measuring here is the number of photons, in R-s, collected at each pixel in each elvetimeind.
	    }
	    
	  }
	}
      }
      
    } // if doelve
    
      // ----------------------------------------
    
      // J update equations: I'm going to ignore the walls and let J = 0. The fields should be near zero anyways since it's all PML.
    
    if (doionosphere) {
#pragma omp parallel for private(Ee1,Ee2,Se1,Ce1,Ce2,Ce3,Ce4,Ae,Ke,Emid,Jer0,Jet0)
      for (int i = nground+1; i < rr-1; i++) {
	for (int j = 1; j < hh-1; j++) {
	  for (int k = 1; k < pp-1; k++) {

	    Ee1 = exp(-nue[i][j][k]*dt);
	    Ee2 = 1/(wce0*wce0 + nue[i][j][k]*nue[i][j][k]);
	    if (wce0 == 0) {
	      Se1 = 1;
	      Ce1 = 0;
	    } else {
	      Se1 = sin(wce0*dt)/wce0;
	      Ce1 = (1 - cos(wce0*dt))/(wce0*wce0);
	    }
	    Ce2 = (1 - Ee1)/nue[i][j][k] - Ee1*nue[i][j][k]*Ce1 - Ee1*Se1;
	    Ce3 = nue[i][j][k] * (1 - Ee1*cos(wce0*dt)) + Ee1*wce0*sin(wce0*dt);
	    Ce4 = 1 - Ee1*cos(wce0*dt) - Ee1*nue[i][j][k]*Se1;
	    
	    // A and K matrices
	    Ae[0][0] = Ee1 * ( Ce1*wcer*wcer + cos(wce0*dt) );
	    Ae[0][1] = Ee1 * ( Ce1*wcer*wcet - Se1*wcep );
	    Ae[0][2] = Ee1 * ( Ce1*wcer*wcep + Se1*wcet );
	    Ae[1][0] = Ee1 * ( Ce1*wcet*wcer + Se1*wcep );     
	    Ae[1][1] = Ee1 * ( Ce1*wcet*wcet + cos(wce0*dt) );
	    Ae[1][2] = Ee1 * ( Ce1*wcet*wcep - Se1*wcer );
	    Ae[2][0] = Ee1 * ( Ce1*wcep*wcer - Se1*wcet );
	    Ae[2][1] = Ee1 * ( Ce1*wcep*wcet + Se1*wcer );
	    Ae[2][2] = Ee1 * ( Ce1*wcep*wcep + cos(wce0*dt) );
            
	    Ke[0][0] = Ee2 * ( Ce2*wcer*wcer + Ce3 );
	    Ke[0][1] = Ee2 * ( Ce2*wcer*wcet - Ce4*wcep );
	    Ke[0][2] = Ee2 * ( Ce2*wcer*wcep + Ce4*wcet );
	    Ke[1][0] = Ee2 * ( Ce2*wcet*wcer + Ce4*wcep );
	    Ke[1][1] = Ee2 * ( Ce2*wcet*wcet + Ce3 );
	    Ke[1][2] = Ee2 * ( Ce2*wcet*wcep - Ce4*wcer );
	    Ke[2][0] = Ee2 * ( Ce2*wcep*wcer - Ce4*wcet );
	    Ke[2][1] = Ee2 * ( Ce2*wcep*wcet + Ce4*wcer );
	    Ke[2][2] = Ee2 * ( Ce2*wcep*wcep + Ce3 );

	    // okay, updates for J finally.

	    Emid[0] = ( dr[i-1] * Er[i][j][k] + dr[i] * Er[i-1][j][k] ) / (dr[i] + dr[i-1]);
	    Emid[1] = (Et[i][j][k] + Et[i][j-1][k])/2;
	    Emid[2] = (Ep[i][j][k] + Ep[i][j][k-1])/2;
	    
	    Jer0 = Ae[0][0] * Jer[i][j][k] + Ae[0][1] * Jet[i][j][k] + Ae[0][2] * Jep[i][j][k] \
	      + E0 * wpe[i][j][k]*wpe[i][j][k] * ( Ke[0][0] * Emid[0] +  Ke[0][1] * Emid[1] + Ke[0][2] * Emid[2] );
	    
	    Jet0 = Ae[1][0] * Jer[i][j][k] + Ae[1][1] * Jet[i][j][k] + Ae[1][2] * Jep[i][j][k] \
	      + E0 * wpe[i][j][k]*wpe[i][j][k] * ( Ke[1][0] * Emid[0] +  Ke[1][1] * Emid[1] + Ke[1][2] * Emid[2] );
	    
	    Jep[i][j][k] = Ae[2][0] * Jer[i][j][k] + Ae[2][1] * Jet[i][j][k] + Ae[2][2] * Jep[i][j][k] \
	      + E0 * wpe[i][j][k]*wpe[i][j][k] * ( Ke[2][0] * Emid[0] +  Ke[2][1] * Emid[1] + Ke[2][2] * Emid[2] );
	    
	    Jer[i][j][k] = Jer0;
	    Jet[i][j][k] = Jet0;
	    
	  }
	}
      }  
      
      // do heating
      for (int i = nground+2; i < rr; i++) {
	for (int j = 1; j < hh; j++) {
	  for (int k = 1; k < pp; k++) {
	    JdotE[i][j][k] = Jer[i][j][k] * (dr[i-1] * Er[i][j][k] + dr[i] * Er[i-1][j][k])/(dr[i]+dr[i-1]) + Jet[i][j][k] * (Et[i][j][k] + Et[i][j-1][k])/2 + Jep[i][j][k] * (Ep[i][j][k] + Ep[i][j][k-1])/2;
	    heat[i][j][k] = heat[i][j][k] + dt * JdotE[i][j][k];
	  }
	}
      }
      
      // take a slice of currents
      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < hh; j++) {
	  Jslicep[i][j] = Jer[i][j][sp];
	}
      }
      for (int i = 0; i < rr; i++) {
	for (int k = 0; k < pp; k++) {
	  Jslicet[i][k] = Jer[i][st][k];
	}
      }
      for (int j = 0; j < hh; j++) {
	for (int k = 0; k < pp; k++) {
	  Jslicer[j][k] = Jer[stepi+10][j][k];
	}
      }
      
    } // do ionosphere


    // VLF Transmitter induced heating, using J.E heating and cooling from Rodriquez (1994) (updated). 
    // ---------------------------------------------------------------------

    if (dotransmitter & doioniz) {
    
      // first index for heating: start at 40 km
      int txfirsti = round(60e3/dr1 + nground);

#pragma omp parallel for private(f_N2,f_O2,gg,Lelast_N2,Lelast_O2,Lrot_N2,Lrot_O2,Lvib_N2,Lvib_O2,Le)
      for (int i = txfirsti; i < rr; i++) {
	for (int j = 1; j < hh; j++) {
	  for (int k = 1; k < pp; k++) { 

	    // cooling
	    f_N2 = 1.06e4 + 7.51e3*tanh(0.0011*(Te[i][j][k]-1800.0));
	    f_O2 = 3300.0 - 839.0*sin(0.000191*(Te[i][j][k]-2700.0));
	    gg = 3300.0 + 1.233*(Te[i][j][k]-1000) - (2.056e-4)*(Te[i][j][k]-1000.0)*(Te[i][j][k]-4000.0);
	    Lelast_N2 = (1.89e-44)*ne[i][j][k]*(0.78*nd3[i][j][k])*(1.0 - 1.21e-4*Te[i][j][k])*Te[i][j][k]*(Te[i][j][k]-Te0[i][j][k]);
	    Lelast_O2 = (1.29e-43)*ne[i][j][k]*(0.21*nd3[i][j][k])*(1.0 + 3.6e-2*sqrt(Te[i][j][k]))*sqrt(Te[i][j][k])*(Te[i][j][k]-Te0[i][j][k]);
	    Lrot_N2 = (4.65e-39)*ne[i][j][k]*(0.78*nd3[i][j][k])*(Te[i][j][k]-Te0[i][j][k])/sqrt(Te[i][j][k]);
	    Lrot_O2 = (1.11e-38)*ne[i][j][k]*(0.21*nd3[i][j][k])*(Te[i][j][k]-Te0[i][j][k])/sqrt(Te[i][j][k]);
	    Lvib_N2 = (4.79e-37)*ne[i][j][k]*(0.78*nd3[i][j][k])*exp(f_N2*(Te[i][j][k]-2000.0)/(2000.0*Te[i][j][k])) * (1-exp(-gg*(Te[i][j][k]-Te0[i][j][k])/(Te[i][j][k]*Te0[i][j][k])));
	    Lvib_O2 = (8.32e-38)*ne[i][j][k]*(0.21*nd3[i][j][k])*exp(f_O2*(Te[i][j][k]-700.0)/(700.0*Te[i][j][k])) * (1-exp(-2700.0*(Te[i][j][k]-Te0[i][j][k])/(Te[i][j][k]*Te0[i][j][k])));
	    Le = Lelast_N2 + Lrot_N2 + Lvib_N2+ Lelast_O2 + Lrot_O2 + Lvib_O2;
	    
	    Te[i][j][k] = Te[i][j][k] + dt * (2.0/(3.0*ne[i][j][k]*KB)) * (JdotE[i][j][k] - Le);
	    
	    // update collision frequency
	    nue[i][j][k] = 1.6 * ( 2.33e-17*(0.78*nd3[i][j][k])*(1.0 - 1.25e-4*Te[i][j][k])*Te[i][j][k] + 1.82e-16*(0.21*nd3[i][j][k])*(1.0 + 3.6e-2*sqrt(Te[i][j][k]))*sqrt(Te[i][j][k]) );
	    
	  }
	}
      }
      
      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < hh; j++) {
	  Teslice[i][j] = Te[i][j][sp];
	  nuslice[i][j] = nue[i][j][sp];
	}
      }
      
    } // if do transmitter
    


    ////////////////////////////////////////////////////////////////////////
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
      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < hh; j++) {
	  for(int k = 0; k < pp; k++) {
	    if (fabs(Er[i][j][k]) > Emax) {
	      Emax = fabs(Er[i][j][k]);
	    }
	  }
	}
      }
      logfile << "Maximum abs(Er) is " << Emax << "\n";
      logfile.close();

      // write to files

      FILE * filePtr;
 
      // E, J, and H are sliced in three directions, using the component selected.
      if (savefields[0]) {
	filePtr = fopen("output_E.dat","ab");
	fwrite(Eslicep,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(Eslicet,sizeof(DATATYPE),rr*pp,filePtr);
	fwrite(Eslicer,sizeof(DATATYPE),hh*pp,filePtr);
	fclose(filePtr);
      }

      if (savefields[1]) {
	filePtr = fopen("output_J.dat","ab");
	fwrite(Jslicep,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(Jslicet,sizeof(DATATYPE),rr*pp,filePtr);
	fwrite(Jslicer,sizeof(DATATYPE),hh*pp,filePtr);
	fclose(filePtr);
      }

      if (savefields[2]) {
	filePtr = fopen("output_H.dat","ab");
	fwrite(Hslicep,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(Hslicet,sizeof(DATATYPE),rr*pp,filePtr);
	fwrite(Hslicer,sizeof(DATATYPE),hh*pp,filePtr);
	fclose(filePtr);
      }

      // remaining fields are always sliced in r-theta plane, at phi = 0.

      if (savefields[3]) {
	filePtr = fopen("output_K.dat","ab");
	fwrite(Eeffslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(Ekslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(heatslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(Sslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(Teslice,sizeof(DATATYPE),rr*hh,filePtr);
	fclose(filePtr);
      }
      
      if (savefields[4]) {
	filePtr = fopen("output_D.dat","ab");
	fwrite(neslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(nOmslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(nuslice,sizeof(DATATYPE),rr*hh,filePtr);
	fclose(filePtr);
      }

      if (savefields[5]) {
	filePtr = fopen("output_O.dat","ab");
	fwrite(nN21Pslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(nN22Pslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(nN2P1Nslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(nN2PMslice,sizeof(DATATYPE),rr*hh,filePtr);
	fwrite(nO2P1Nslice,sizeof(DATATYPE),rr*hh,filePtr);
	fclose(filePtr);
      }

    }

    // end of big time loop
  }
    
 endofbigloop:

  /////////////////////////////////////////////////////////////////////////////

  // write the elve

  if (doelve) {
    logfile.open("log.txt",std::fstream::app);
    logfile << "Now writing elve data to elve.dat\n";
    logfile.close();

    FILE * elveFile;
    
    elveFile = fopen("elve.dat","wb");
    
    fwrite(&elveN21P,sizeof(DATATYPE),numpixels*elvesteps,elveFile);
    fwrite(&elveN22P,sizeof(DATATYPE),numpixels*elvesteps,elveFile);
    fwrite(&elveN2P1N,sizeof(DATATYPE),numpixels*elvesteps,elveFile);
    fwrite(&elveN2PM,sizeof(DATATYPE),numpixels*elvesteps,elveFile);
    fwrite(&elveO2P1N,sizeof(DATATYPE),numpixels*elvesteps,elveFile);
    fclose(elveFile);
    
  }
  
  FILE * ProbeFile;
  ProbeFile = fopen("Probe.dat","wb");
  fwrite(&nprobes,sizeof(int),1,ProbeFile);
  fwrite(&prober,sizeof(int),nprobes,ProbeFile);
  fwrite(&probet,sizeof(int),nprobes,ProbeFile);
  fwrite(&probep,sizeof(int),nprobes,ProbeFile);
  fwrite(&Erprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Etprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Epprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Hrprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Htprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fwrite(&Hpprobe,sizeof(DATATYPE),tsteps*nprobes,ProbeFile);
  fclose(ProbeFile);

  logfile.open("log.txt",std::fstream::app);
  logfile << "All done!\n";

  double totalruntime = (omp_get_wtime() - tStart) / 60.0;

  logfile << "Total run time = " << totalruntime << "minutes.\n";
  logfile.close();

  // end program - do not type below this!

}

//////////////////////////////////////////////////////////
