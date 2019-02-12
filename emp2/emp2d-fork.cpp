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

// version 2.0.1: November 2, 2013: now includes 2D neutral density
// version 2.0.2: February 27, 2014: includes Surface Impedance Boundary Condition (SIBC)
// version 2.0.3: April 2014: now gives option to output six different savefields files individually, with a vector of six flags input.
// version 2.0.4: May 2014: implified venus option to make it external to code. Just load ne, nd, and rates profiles as you do for earth. difference will lie in generation of those files in matlab setup.
// version 3.0: modified to read source from file as alt-vs-time 2D array. 
// version 3.1: added DFT computations, to get transmitter response from impulsive source
// fork: adding ability to read 2D ionosphere from file.

///////////////////////////////////////////////////////////////////////////////////

int main()
{

  // time the whole thing
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
  int numfiles;
  int planet;
  int decfactor;
  int nprobes;
  int dogwave;
  double gwavemag;
  double gwavemaxalt;
  double gwavekh;
  int nonlinearstart;  // index to start nonlinear calculations
  int doDFT; // output amplitudes along ground to extract amplitude and phase later
  int numDFTfreqs;
  int read2Dionosphere;

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
  fread(&numfiles,sizeof(int),1,inputFile);
  fread(&planet,sizeof(int),1,inputFile);
  fread(&decfactor,sizeof(int),1,inputFile);
  fread(&nprobes,sizeof(int),1,inputFile);
  int prober[nprobes], probet[nprobes];
  fread(&prober,sizeof(int),nprobes,inputFile);
  fread(&probet,sizeof(int),nprobes,inputFile);
  fread(&dogwave,sizeof(int),1,inputFile);
  fread(&gwavemag,sizeof(double),1,inputFile);
  fread(&gwavemaxalt,sizeof(double),1,inputFile); 
  fread(&gwavekh,sizeof(double),1,inputFile); 
  fread(&nonlinearstart,sizeof(int),1,inputFile);
  fread(&doDFT,sizeof(int),1,inputFile);
  fread(&numDFTfreqs,sizeof(int),1,inputFile);
  double DFTfreqs[numDFTfreqs];
  fread(&DFTfreqs,sizeof(double),numDFTfreqs,inputFile);
  fread(&read2Dionosphere,sizeof(int),1,inputFile);

  fclose(inputFile);

  ofstream logfile;

  logfile.open("log.txt"); 
  logfile << "Successfully read inputs file!\n";
  logfile.close();

  //------------------------------------------------------------------------//
  //-------READ SOURCE------------------------------------------------------//

  FILE * sourceFile;
  sourceFile = fopen("source.dat","rb");

  int nsalts;
  int nstimes;
  int halfchlength; // ignored by 2D code
  fread(&nsalts,sizeof(int),1,sourceFile);
  fread(&nstimes,sizeof(int),1,sourceFile);
  fread(&halfchlength,sizeof(int),1,sourceFile);
  double Isource [nstimes][nsalts];
  fread(&Isource,sizeof(double),nsalts*nstimes,sourceFile);
  fclose(sourceFile);

  logfile.open("log.txt", std::ofstream::app);
  logfile << "Successfully read source file!\n";
  logfile.close();
    
  //------------------------------------------------------------------------//
  //-------GRID PARAMETERS -------------------------------------------------//

  // set up radial array
  
  // DATATYPE ngr = (DATATYPE)nground;
  int rr = stepalt/dr1 + (maxalt - stepalt)/dr2 + 1 + nground;
  int stepi = stepalt/dr1;
  DATATYPE r [rr];
  DATATYPE dr [rr-1];
  
  r[0] = RE - nground*dr0;
  if (nground > 0) {
    for (int i = 1; i < nground+1; i++) {
      r[i] = r[i-1] + dr0;
      dr[i-1] = r[i] - r[i-1];
    }
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
  if (thmax > PI) {
    thmax = PI;
  }
  DATATYPE dth = drange / RE;
  int hh = round(thmax / dth) + 1;
  DATATYPE th [hh];
  th[0] = 0;  
  for (int i = 1; i < hh; i++) {
    th[i] = th[i-1] + dth;
  }

  // decimation for output fields
  int drr = rr;
  int dhh = hh;
  if (decfactor > 1) {
    drr = floor(rr/decfactor);
    dhh = floor(hh/decfactor);
  }

  

  //------------------------------------------------------------------------//
  //-------GROUND PARAMETERS -----------------------------------------------//

  // input groundmethod can be 0 (PEC), 1 (SIBC), or 2 (Real parameters)

  // read ground sigma and epsilon
  double gsig [hh];
  double geps [hh];
    
  FILE * groundFile;
  groundFile = fopen("ground.dat","rb");
  fread(&gsig,sizeof(double),hh,groundFile);
  fread(&geps,sizeof(double),hh,groundFile);
  fclose(groundFile);

  logfile.open("log.txt");
  logfile << "Successfully read ground file!\n";
  logfile.close();
    
  //------------------------------------------------------------------------//
  //-------MULTIPLY FACTORS FOR FDTD ---------------------------------------//
  
    
  DATATYPE c1h = (2*U0 - sigm*dt)/(2*U0 + sigm*dt);
  DATATYPE c2h = 2*dt/(2*U0 + sigm*dt);
  DATATYPE c1e [rr][hh];
  DATATYPE c2e [rr][hh];
  
    
  for (int i = 0; i < rr; i++) {
    for (int j = 0; j < hh; j++) {
      c1e[i][j] = (2*E0 - sig*dt)/(2*E0 + sig*dt);
      c2e[i][j] = 2*dt/(2*E0 + sig*dt);
    }
  }
  
  if (groundmethod == 2) {
    for (int i = 0; i < nground; i++) {
      for (int j = 0; j < hh; j++) {
	c1e[i][j] = (2*E0*geps[j] - gsig[j]*dt)/(2*E0*geps[j] + gsig[j]*dt);
	c2e[i][j] = 2*dt/(2*E0*geps[j] + gsig[j]*dt);
      }
    }
  }

  //-----------------------------------------------------------------------//
  //----   Surface Impedance Boundary Condition (SIBC) stuff --------------//
  //----   This is taken from Oh and Schutt-Aine, IEEE TRANSACTIONS ON ANTENNAS AND PROPAGATION, VOL. 43, NO. 7, JULY 1995 

  DATATYPE sibc_a [hh];
  DATATYPE sibc_eta [hh];
  DATATYPE sibc_C [6] = { 1.22646e-8, 2.56716e-6, 1.51777e-4, 4.42437e-3, 6.98268e-2, 0.42473 };
  DATATYPE sibc_omega [6] = { 4.06981e-6, 1.84651e-4, 3.24245e-3, 3.42849e-2, 0.23606, 0.83083 };
  DATATYPE sibc_K;
  DATATYPE sibc_pi1 [hh][6];
  DATATYPE sibc_pi2 [hh][6];
  DATATYPE sibc_pi3 [hh][6];

  DATATYPE sumAi;
  DATATYPE sibc_Ai [hh][6];
  DATATYPE sibc_Aiold [hh][6];

  // will also need "old" tangential H fields
  DATATYPE Hpold [rr][hh];
  DATATYPE Htold [rr][hh];
    

  for (int j = 0; j < hh; j++) {
    sibc_a[j] = gsig[j]/(E0*geps[j]);
    sibc_eta[j] = sqrt(U0/(E0*geps[j]));
    for (int m = 0; m < 6; m++) {
      sibc_K = sibc_a[j] * sibc_omega[m] * dt;
      sibc_pi3[j][m] = exp(-sibc_K);
      sibc_pi1[j][m] = sibc_eta[j] * (sibc_C[m]/sibc_omega[m]) * (1 + (exp(-sibc_K) - 1)/sibc_K);
      sibc_pi2[j][m] = sibc_eta[j] * (sibc_C[m]/sibc_omega[m]) * (1/sibc_K - exp(-sibc_K) * (1 + 1/sibc_K));
    }
  }

  memset(sibc_Ai,0,sizeof(sibc_Ai));
  memset(sibc_Aiold,0,sizeof(sibc_Aiold));
  memset(Hpold,0,sizeof(Hpold));
  memset(Htold,0,sizeof(Htold));


  //------------------------------------------------------------------------//
  //-------MAGNETIC FIELD --------------------------------------------------//

  // read magnetic field: varying in theta direction
  double Br [hh];
  double Bt [hh];
  double Bp [hh];
  FILE * BFile;
  BFile = fopen("B0.dat","rb");
  fread(&Br,sizeof(double),hh,BFile);
  fread(&Bt,sizeof(double),hh,BFile);
  fread(&Bp,sizeof(double),hh,BFile);
  fclose(BFile);

  logfile.open("log.txt");
  logfile << "Successfully read B0 file!\n";
  logfile.close();
    
  // compute gyrofrequencies, store in 2D arrays

  DATATYPE wcer [hh];
  DATATYPE wcet [hh];
  DATATYPE wcep [hh];
  DATATYPE wce0 [hh];
  DATATYPE wcir [hh];
  DATATYPE wcit [hh];
  DATATYPE wcip [hh];
  DATATYPE wci0 [hh];

  for (int j = 0; j < hh; j++) {
    wcer[j] = -QE * Br[j] / ME;
    wcet[j] = -QE * Bt[j] / ME;
    wcep[j] = -QE * Bp[j] / ME;
    wce0[j] = sqrt(wcer[j]*wcer[j] + wcet[j]*wcet[j] + wcep[j]*wcep[j]);
    wcir[j] = QE * Br[j] / MI;
    wcit[j] = QE * Bt[j] / MI;
    wcip[j] = QE * Bp[j] / MI;
    wci0[j] = sqrt(wcir[j]*wcir[j] + wcit[j]*wcit[j] + wcip[j]*wcip[j]);
  }


  //------------------------------------------------------------------------//
  //-------OUTPUT TO LOG FILE ----------------------------------------------//

  logfile.open("log.txt");

  logfile << "\n";
  logfile << "------------------------------------------------------------------\n";
  logfile << " ------------------- 2D EMP Simulation ---------------------------\n";
  logfile << "DoPML = " << dopml_top << " (top), " << dopml_wall << " (wall)\n";
  logfile << "DoIonosphere = " << doionosphere << "\n";
  logfile << "DoIoniz = " << doioniz << "\n";
  logfile << "Grid is " << rr << " (in r) by " << hh << " (in theta) cells\n";
  logfile << "Output grid is " << drr << " (in r) by " << dhh << " (in theta) cells\n";
  logfile << "Ground conductivity at source is " << gsig[0] << "\n";
  logfile << "Ground permittivity at source is " << geps[0] << "\n";
  logfile << "Br on axis is " << Br[0] << "; Bt on axis is " << Bt[0] << "\n";
  logfile << "Time step is " << (dt*1e9) << " ns\n";
  logfile << "Camera located " << camdist/1e3 << " km away, at " << camalt/1e3 << " km altitude\n"; 
  logfile << "drange is " << drange << "\n";
  logfile << "Simulation will run " << tsteps << " time steps.\n";
  logfile << "Nonlinear calculations start at i = " << nonlinearstart << "\n";
  logfile << "------------------------------------------------------------------\n\n";

  // write space parameters to sferic.dat file

  FILE * sfericFile;

  sfericFile = fopen("sferic.dat","wb");
  fwrite(&tsteps,sizeof(int),1,sfericFile);
  fwrite(&rr,sizeof(int),1,sfericFile);
  fwrite(&hh,sizeof(int),1,sfericFile);
  fwrite(&numfiles,sizeof(int),1,sfericFile);
  fwrite(&dt,sizeof(DATATYPE),1,sfericFile);
  fwrite(r,sizeof(DATATYPE),rr,sfericFile);
  fwrite(th,sizeof(DATATYPE),hh,sfericFile);
  fwrite(&decfactor,sizeof(int),1,sfericFile);
  fclose(sfericFile);

  logfile.open("log.txt");
  logfile << "Successfully wrote sferic.dat file!\n";
  logfile.close();

  //------------------------------------------------------------------------//
  //-------ELVE INTEGRATION SETUP ------------------------------------------//

  //int elvewritesteps = 1000;    // this is the number of time steps at which photons will be written to the elve cube.

  int numpixels;
  int cameratype;
  FILE * cameraFile;

  cameraFile = fopen("camera.dat","rb");
  fread(&cameratype,sizeof(int),1,cameraFile);
  fread(&numpixels,sizeof(int),1,cameraFile);
  DATATYPE elveaz [numpixels];
  DATATYPE elveel [numpixels];
  fread(&elveaz,sizeof(double),numpixels,cameraFile);
  fread(&elveel,sizeof(double),numpixels,cameraFile);
  fclose(cameraFile);

  logfile.open("log.txt");
  logfile << "Successfully read camera file!\n";
  logfile.close();

  //int elvesteps = 1000;         // this is the number of time steps in the actual elve cube
  DATATYPE eti = (camdist - range)/C;
  DATATYPE etf = (camdist + range)/C + tsteps*dt;
  DATATYPE elvedt = (etf - eti)/(elvesteps-1);

  //int elvesize [2] = {128, 64};

  DATATYPE elveN21P [numpixels][elvesteps];
  DATATYPE elveN22P [numpixels][elvesteps];
  DATATYPE elveN2P1N [numpixels][elvesteps];
  DATATYPE elveN2PM [numpixels][elvesteps];
  DATATYPE elveO2P1N [numpixels][elvesteps];
  DATATYPE elveOred [numpixels][elvesteps];
  DATATYPE elveOgrn [numpixels][elvesteps];

  logfile << "Size of elve array = " << sizeof(elveN21P) << "\n";

  memset(elveN21P,0,sizeof(elveN21P));
  memset(elveN22P,0,sizeof(elveN22P));
  memset(elveN2P1N,0,sizeof(elveN22P));
  memset(elveN2PM,0,sizeof(elveN22P));
  memset(elveO2P1N,0,sizeof(elveN22P));
  memset(elveOred,0,sizeof(elveOred));
  memset(elveOgrn,0,sizeof(elveOgrn));

  // need arrays thind, tdelay, and raylength to use for integration

  int thind [numpixels][rr];
  DATATYPE tdelay [numpixels][rr];
  DATATYPE raylength [numpixels][rr];

  DATATYPE d, az, el, th2, theta;
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
      theta = acos( cos(th2)*cos(camth) + sin(th2)*sin(camth)*cos(az) );
      thind[m][k] = floor(theta/dth);
      raylength[m][k] = dr2 / sin(el+th2);
    }
  }
  

  //------------------------------------------------------------------------//
  //------- FIELD INITIALIZATION --------------------------------------------//

  // electric field vector
  DATATYPE Er [rr][hh];
  DATATYPE Et [rr][hh];
  DATATYPE Ep [rr][hh];
  // magnetic field vector
  DATATYPE Hr [rr][hh];
  DATATYPE Ht [rr][hh];
  DATATYPE Hp [rr][hh];
  // spatial-averaged electric field
  DATATYPE Erm, Etm, Hrm, Htm, Hpm;

  // DFT fields
  DATATYPE DFTfieldsEr [hh][2*numDFTfreqs];
  DATATYPE DFTfieldsEt [hh][2*numDFTfreqs];
  DATATYPE DFTfieldsEp [hh][2*numDFTfreqs];
  DATATYPE DFTfieldsHr [hh][2*numDFTfreqs];
  DATATYPE DFTfieldsHt [hh][2*numDFTfreqs];
  DATATYPE DFTfieldsHp [hh][2*numDFTfreqs];

  // probe fields
  DATATYPE Erprobe [tsteps][nprobes];
  DATATYPE Etprobe [tsteps][nprobes];
  DATATYPE Epprobe [tsteps][nprobes];
  DATATYPE Hrprobe [tsteps][nprobes];
  DATATYPE Htprobe [tsteps][nprobes];
  DATATYPE Hpprobe [tsteps][nprobes];
  int ir, it;  // probe point indices

  // electron current
  DATATYPE Jer [rr][hh];
  DATATYPE Jet [rr][hh];
  DATATYPE Jep [rr][hh];
  DATATYPE Jer0 = 0;
  DATATYPE Jet0 = 0;

  // ion current
  DATATYPE Jir [rr][hh];
  DATATYPE Jit [rr][hh];
  DATATYPE Jip [rr][hh];
  DATATYPE Jir0 = 0;
  DATATYPE Jit0 = 0;

  // total current
  DATATYPE Jr [rr][hh];
  DATATYPE Jt [rr][hh];
  DATATYPE Jp [rr][hh];
  DATATYPE heat [rr][hh];   // heat dissipated through J.E

  // electron temperature
  DATATYPE Te [rr][hh];
  DATATYPE Te0 [rr][hh];

  // Electric field Parallel and Perpendicular components, for time-dependent stuff
  DATATYPE Epar, Eperp2, HpmR, HpmL;
  DATATYPE Eeff [rr][hh];
  DATATYPE Emag [rr][hh];

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
  memset(Jir,0,sizeof(Jir));
  memset(Jit,0,sizeof(Jit));
  memset(Jip,0,sizeof(Jip));
  memset(Jr,0,sizeof(Jr));
  memset(Jt,0,sizeof(Jt));
  memset(Jp,0,sizeof(Jp));
  memset(Eeff,0,sizeof(Eeff));
  memset(Emag,0,sizeof(Emag));
  memset(heat,0,sizeof(heat));

  memset(DFTfieldsEr,0,sizeof(DFTfieldsEr));
  memset(DFTfieldsEt,0,sizeof(DFTfieldsEt));
  memset(DFTfieldsEp,0,sizeof(DFTfieldsEp));
  memset(DFTfieldsHr,0,sizeof(DFTfieldsHr));
  memset(DFTfieldsHt,0,sizeof(DFTfieldsHt));
  memset(DFTfieldsHp,0,sizeof(DFTfieldsHp));
  
  // decimated fields, to reduce output file size

  DATATYPE Erdec[drr][dhh];
  DATATYPE Etdec[drr][dhh];
  DATATYPE Epdec[drr][dhh];
  DATATYPE Hrdec[drr][dhh];
  DATATYPE Htdec[drr][dhh];
  DATATYPE Hpdec[drr][dhh];
  DATATYPE Jerdec[drr][dhh];
  DATATYPE Jetdec[drr][dhh];
  DATATYPE Jepdec[drr][dhh];
  DATATYPE Jirdec[drr][dhh];
  DATATYPE Jitdec[drr][dhh];
  DATATYPE Jipdec[drr][dhh];

  DATATYPE Eeffdec[drr][dhh];
  DATATYPE Ekdec[drr][dhh];
  DATATYPE heatdec[drr][dhh];
  DATATYPE Sdec[drr][dhh];
  DATATYPE Tedec[drr][dhh];
  DATATYPE nedec[drr][dhh];
  DATATYPE nOmdec[drr][dhh];
  DATATYPE nuedec[drr][dhh];
  DATATYPE nN21Pdec[drr][dhh];
  DATATYPE nN22Pdec[drr][dhh];
  DATATYPE nN2P1Ndec[drr][dhh];
  DATATYPE nN2PMdec[drr][dhh];
  DATATYPE nO2P1Ndec[drr][dhh];
  DATATYPE nOreddec[drr][dhh];
  DATATYPE nOgrndec[drr][dhh];

  // some stuff for the current source

  int nextra = 5;
  DATATYPE jfl;

  //------------------------------------------------------------------------//
  //-------NUMBER DENSITIES ------------------------------------------------//

  // set up electron density and B field

  double ne1 [rr];       // 1D electron density
  DATATYPE ne [rr][hh];  // 2D electron density
  DATATYPE wpe [rr][hh]; // plasma frequency
  double nd [rr];        // neutral density
  DATATYPE mue;          // mobility
  DATATYPE nue [rr][hh]; // collision frequency
  DATATYPE Ek [rr][hh];  // breakdown field

  // do Ions as well!
  double ni1 [rr];
  DATATYPE ni [rr][hh];
  DATATYPE wpi [rr][hh];
  DATATYPE nui [rr][hh];

  FILE * neFile;
  neFile = fopen("ne.dat","rb");
  if (read2Dionosphere) {
    fread(&ne,sizeof(double),hh*rr,neFile);
  } else {
    fread(&ne1,sizeof(double),rr,neFile);
  }
  fclose(neFile);

  FILE * niFile;
  niFile = fopen("ni.dat","rb");
  fread(&ni1,sizeof(double),rr,niFile);
  fclose(niFile);

  // read ambient temperature
  double Tamb [rr];
  FILE * Tefile;
  Tefile = fopen("etemp.dat","rb");
  fread(&Tamb,sizeof(double),rr,Tefile);
  fclose(Tefile);

  // extrapolate 2D arrays
  for (int i = 0; i < rr; i++) {
    for (int j = 0; j < hh; j++) {
      if (!read2Dionosphere) {
	ne[i][j] = ne1[i];
      }
      ni[i][j] = ni1[i];
      wpe[i][j] = QE * sqrt(ne[i][j] / (ME * E0));
      wpi[i][j] = QE * sqrt(ni[i][j] / (MI * E0));
      Te0[i][j] = Tamb[i];
      Te[i][j] = Tamb[i];
    }
  }

  // temperature intermediates; private variables in openmp, except S and components

  DATATYPE S [rr][hh];
  DATATYPE Sr,St,Sp;
  DATATYPE gg,f_N2,f_O2,Lelast_N2,Lelast_O2,Lrot_N2,Lrot_O2,Lvib_N2,Lvib_O2,Le;
  DATATYPE JdotE [rr][hh];

  memset(S,0,sizeof(S));
  memset(JdotE,0,sizeof(JdotE));

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

  // create 2D nd array in order to input gravity wave
  DATATYPE nd2 [rr][hh];
  DATATYPE gwfac;
  for (int i = 0; i < rr; i++) {
    for (int j = 0; j < hh; j++) {
      gwfac = gwavemag / exp(gwavemaxalt/12e3) * exp((r[i]-RE)/12e3) * sin(gwavekh*th[j]*RE);
      if (gwfac > gwavemag) gwfac = gwavemag;
      if (gwfac < -gwavemag) gwfac = -gwavemag;
      if (dogwave) {
	nd2[i][j] = nd[i] * ( 1 + gwfac );
      } else {
	nd2[i][j] = nd[i];
      }
    }
  }

  // now use 2D nd to create 2D mue, Ek
  for (int i = 0; i < rr; i++) {
    for (int j = 0; j < hh; j++) {
      if (planet == 1) {   // venus
	Ek[i][j] = 7.07e7 * nd2[i][j] / nd[nground];
	mue = 0.0018 * nd[nground] / nd2[i][j];
      } else if (planet == 0) {  // earth
	Ek[i][j] = 3.0e6 * nd2[i][j] / nd[nground];
	mue = muArray[i][0] / nd2[i][j];
      } else if (planet == 2) {  // saturn
	Ek[i][j] = 7.07e7 * nd2[i][j] / nd[nground];
	mue = 0.0018 * nd[nground] / nd2[i][j];
      } else {
	// assume earth
	Ek[i][j] = 3.0e6 * nd2[i][j] / nd[nground];
	mue = muArray[i][0] / nd2[i][j];
      }
      nue[i][j] = (QE / ME) / mue;
      nui[i][j] = nue[i][j] / 100;
    }
  }

  // if we are using a 2D ionosphere, we also need a 2D collision frequency (perturbation studies)
  
  if (read2Dionosphere) {
  FILE * nuFile;
  nuFile = fopen("nu.dat","rb");
    fread(&nue,sizeof(double),hh*rr,nuFile);
    for (int i = 0; i < rr; i++) {
      for (int j = 0; j < hh; j++) {
	nui[i][j] = nue[i][j] / 100;
      }
    }
  }


  // if do transmitter, collision frequency depends on Temperature and density

  if (dotransmitter) {
    for (int i = 0; i < rr; i++) {
      for (int j = 0; j < hh; j++) {
	nue[i][j] = 1.6 * ( 2.33e-17*(0.78*nd2[i][j])*(1.0 - 1.25e-4*Te[i][j])*Te[i][j] + 1.82e-16*(0.21*nd2[i][j])*(1.0 + 3.6e-2*sqrt(Te[i][j]))*sqrt(Te[i][j]) );
	// override!
	nue[i][j] = 1.815775e11 * exp(-0.15*(r[i]-RE)/1000);
	nui[i][j] = nue[i][j] / 100;
      }
    }
  }

  //------------------------------------------------------------------------//
  //------- OPTICAL FIELDS -------------------------------------------------//

  DATATYPE EoEk;
  DATATYPE vi, va, vd;  // ionization
  DATATYPE nOm [rr][hh]; // O- ions

  // optical arrays
  DATATYPE vN21P, vN22P, vN2P1N, vN2PM, vO2P1N, vOred, vOgrn;
  DATATYPE nN21P [rr][hh];
  DATATYPE nN22P [rr][hh];
  DATATYPE nN2P1N [rr][hh];
  DATATYPE nN2PM [rr][hh];
  DATATYPE nO2P1N [rr][hh];
  DATATYPE nOred [rr][hh];
  DATATYPE nOgrn [rr][hh];
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
  DATATYPE Emidr;

  // intermediate calculations in Lee&Kalluri
  DATATYPE Ee1, Ee2, Se1, Ce1, Ce2, Ce3, Ce4;
  DATATYPE Ei1, Ei2, Si1, Ci1, Ci2, Ci3, Ci4;


  //------------------------------------------------------------------------//
  //------- SET UP PML -----------------------------------------------------//

  // initialize A, B, and P arrays even if we don't use them.

  double pmlm = 4;
  int pmllen = 10;

  DATATYPE A_ert [pmllen];
  DATATYPE A_etr [pmllen];
  DATATYPE A_epr [pmllen];
  DATATYPE A_ept [pmllen];
  DATATYPE A_hrt [pmllen];
  DATATYPE A_htr [pmllen];
  DATATYPE A_hpr [pmllen];
  DATATYPE A_hpt [pmllen];

  DATATYPE B_ert [pmllen];
  DATATYPE B_etr [pmllen];
  DATATYPE B_epr [pmllen];
  DATATYPE B_ept [pmllen];
  DATATYPE B_hrt [pmllen];
  DATATYPE B_htr [pmllen];
  DATATYPE B_hpr [pmllen];
  DATATYPE B_hpt [pmllen];

  DATATYPE P_etr [pmllen][hh];
  DATATYPE P_ert [rr][pmllen];
  DATATYPE P_epr [pmllen][hh];
  DATATYPE P_ept [rr][pmllen];
  DATATYPE P_hrt [rr][pmllen];
  DATATYPE P_htr [pmllen][hh];
  DATATYPE P_hpr [pmllen][hh];
  DATATYPE P_hpt [rr][pmllen];

  // initialize
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
    for (int m = 0; m < pmllen; m++) {
      B_ert[m] = exp(-((st[m]/kt[m]) + at[m])*dt/E0);
      A_ert[m] = st[m] / (st[m]*kt[m] + pow(kt[m],2)*at[m]) * (B_ert[m] - 1);
      B_etr[m] = exp(-((sr[m]/kr[m]) + ar[m])*dt/E0);
      A_etr[m] = sr[m] / (sr[m]*kr[m] + pow(kr[m],2)*ar[m]) * (B_etr[m] - 1);
      B_ept[m] = exp(-((st[m]/kt[m]) + at[m])*dt/E0);
      A_ept[m] = st[m] / (st[m]*kt[m] + pow(kt[m],2)*at[m]) * (B_ept[m] - 1);
      B_epr[m] = exp(-((sr[m]/kr[m]) + ar[m])*dt/E0);
      A_epr[m] = sr[m] / (sr[m]*kr[m] + pow(kr[m],2)*ar[m]) * (B_epr[m] - 1);
       
      B_hrt[m] = exp(-((stm[m]/ktm[m]) + atm[m])*dt/E0);
      A_hrt[m] = stm[m] / (stm[m]*ktm[m] + pow(ktm[m],2)*atm[m]) * (B_hrt[m] - 1);
      B_htr[m] = exp(-((srm[m]/krm[m]) + arm[m])*dt/E0);
      A_htr[m] = srm[m] / (srm[m]*krm[m] + pow(krm[m],2)*arm[m]) * (B_htr[m] - 1);
      B_hpt[m] = exp(-((stm[m]/ktm[m]) + atm[m])*dt/E0);
      A_hpt[m] = stm[m] / (stm[m]*ktm[m] + pow(ktm[m],2)*atm[m]) * (B_hpt[m] - 1);
      B_hpr[m] = exp(-((srm[m]/krm[m]) + arm[m])*dt/E0);
      A_hpr[m] = srm[m] / (srm[m]*krm[m] + pow(krm[m],2)*arm[m]) * (B_hpr[m] - 1);     

    }

  } // if dopml

  // offsets for PML
  int rshift = rr-pmllen-1;
  int hshift = hh-pmllen-1;

  //-----------------------------------------------------------------------//
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
    if (dopml_top) {
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < hh-1; j++) {
	  P_htr[i][j] = B_htr[i] * P_htr[i][j] \
	    + A_htr[i] * ( r[i+rshift+1] * Ep[i+rshift+1][j] - r[i+rshift] * Ep[i+rshift][j] ) / dr[i+rshift];
	}
      }
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < hh-1; j++) {
	  P_hpr[i][j] = B_hpr[i] * P_hpr[i][j] \
	    + A_hpr[i] * ( r[i+rshift+1] * Et[i+rshift+1][j] - r[i+rshift] * Et[i+rshift][j] ) / dr[i+rshift];
	}
      }
    }
    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  P_hrt[i][j] = B_hrt[j] * P_hrt[i][j] \
	    + A_hrt[j] * ( Ep[i][j+hshift+1] - Ep[i][j+hshift] ) / dth / r[i];
	}
      }
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  P_hpt[i][j] = B_hpt[j] * P_hpt[i][j] \
	    + A_hpt[j] * ( Er[i][j+hshift+1] - Er[i][j+hshift] ) / dth;
	}
      }
      // end Psi updates
    }

    // ----------------------------------------
    // Hr update
#pragma omp parallel for
    for (int i = 0; i < rr; i++ ) {
      for (int j = 0; j < hh-1; j++) {
	Hr[i][j] = c1h * Hr[i][j] \
	  - c2h / (r[i] * sin(th[j]+dth/2)) * ( sin(th[j+1]) * Ep[i][j+1] - sin(th[j]) * Ep[i][j] ) / dth;
      }
    }
    // pml correction
    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  Hr[i][j+hshift] = Hr[i][j+hshift] - c2h * P_hrt[i][j];
	}
      }
    }

    // ----------------------------------------
    // Ht update
#pragma omp parallel for
    for (int i = 0; i < rr-1; i++) {
      for (int j = 0; j < hh; j++) {
	Ht[i][j] = c1h * Htold[i][j] \
	  + c2h / (r[i]+dr[i]/2) * ( r[i+1] * Ep[i+1][j] - r[i] * Ep[i][j] ) / dr[i];
      }
    }
    // pml correction
    if (dopml_top) {
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < hh-1; j++) {
	  Ht[i+rshift][j] = Ht[i+rshift][j] + c2h / (r[i+rshift]+dr[i+rshift]/2) * P_htr[i][j];
	}
      }
    }

    // ----------------------------------------
    // Hp update
#pragma omp parallel for
    for (int i = 0; i < rr-1; i++) {
      for (int j = 0; j < hh-1; j++) {
	Hp[i][j] = c1h * Hpold[i][j] \
	  - c2h / (r[i]+dr[i]/2) * ( ( r[i+1] * Et[i+1][j] - r[i] * Et[i][j] ) / dr[i] - ( Er[i][j+1] - Er[i][j] ) / dth );
      }
    }
    // pml correction
    if (dopml_top) {
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < hh-1; j++) {
	  Hp[i+rshift][j] = Hp[i+rshift][j] - c2h / (r[i+rshift]+dr[i+rshift]/2) * P_hpr[i][j];
	}
      }
    }
    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  Hp[i][j+hshift] = Hp[i][j+hshift] + c2h / (r[i]+dr[i]/2) * P_hpt[i][j];
	}
      }
    }

    // ---------------------------------------
    // Psi updates for E field

    if (dopml_top) {
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < hh-1; j++) {
	  P_etr[i][j] = B_etr[i] * P_etr[i][j] \
	    + A_etr[i] * ( (r[i+rshift]+dr[i+rshift]/2) * Hp[i+rshift][j] \
			   - (r[i+rshift]-dr[i+rshift-1]/2) * Hp[i+rshift-1][j] ) / ((dr[i+rshift]+dr[i+rshift-1])/2);
	}
      }
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < hh-1; j++) {
	  P_epr[i][j] = B_epr[i] * P_epr[i][j] \
	    + A_epr[i] * ( (r[i+rshift]+dr[i+rshift]/2) * Ht[i+rshift][j] \
			   - (r[i+rshift]-dr[i+rshift-1]/2) * Ht[i+rshift-1][j] ) / ((dr[i+rshift]+dr[i+rshift-1])/2);
	}
      }
    }
    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  P_ert[i][j] = B_ert[j] * P_ert[i][j] \
	    + A_ert[j] * ( Hp[i][j+hshift] - Hp[i][j+hshift-1] ) / dth;
	}
      }
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  P_ept[i][j] = B_ept[j] * P_ept[i][j] \
	    + A_ept[j] * ( Hr[i][j+hshift] - Hr[i][j+hshift-1] ) / dth;
	}
      }
      // end Psi updates
    }

    // ----------------------------------------
    // Er update
#pragma omp parallel for
    for (int i = 0; i < rr-1; i++) {
      for (int j = 1; j < hh-1; j++) {
	Er[i][j] = c1e[i][j] * Er[i][j] \
	  + c2e[i][j] / ( (r[i]+dr[i]/2) * sin(th[j]) ) * ( sin(th[j]+dth/2) * Hp[i][j] - sin(th[j]-dth/2) * Hp[i][j-1] ) / dth \
	  - c2e[i][j] * ( Jr[i+1][j] + Jr[i][j] ) / 2;
      }
    }
    // on-axis correction
    for (int i = 0; i < rr-1; i++) {
      Er[i][0] = c1e[i][0] * Er[i][0] \
	+ c2e[i][0] * sin(dth/2) / ( (r[i]+dr[i]/2) * (1 - cos(dth/2)) ) * Hp[i][0] \
	- c2e[i][0] * ( Jr[i+1][0] + Jr[i][0] ) / 2;
    }
    // other axis, if thmax = pi
    if (thmax == PI) {
      for (int i = 0; i < rr-1; i++) {
	Er[i][hh] = c1e[i][hh] * Er[i][hh]					\
	  - c2e[i][hh] * sin(dth/2) / ( (r[i]+dr[i]/2) * (1 - cos(dth/2)) ) * Hp[i][hh-1] \
	  - c2e[i][hh] * ( Jr[i+1][hh] + Jr[i][hh] ) / 2;
      }
    }

    // source ===================================

    if (t < nstimes) {
      for (int i = 0; i < nsalts; i++) {
	for (int j = 0; j < nextra; j++) {
	  jfl = (DATATYPE)j;
	  Er[i][j] = Er[i][j] - c2e[i][j] * Isource[t][i] * exp(-jfl*jfl/9.0) / (PI * 9.0*drange*drange);
	}
      }
    }

    // ==========================================

    // pml corrections
    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  Er[i][j+hshift] = Er[i][j+hshift] + c2e[i][j+hshift] / (r[i]+dr[i]/2) * P_ert[i][j];
	}
      }
    }

    // ----------------------------------------
    // Et update
#pragma omp parallel for
    for (int i = 1; i < rr-1; i++) {
      for (int j = 0; j < hh-1; j++) {
	Et[i][j] = c1e[i][j] * Et[i][j] \
	  - c2e[i][j] / r[i] * ( (r[i]+dr[i]/2) * Hp[i][j] - (r[i]-dr[i-1]/2) * Hp[i-1][j] ) / ((dr[i]+dr[i-1])/2) \
	  - c2e[i][j] * ( Jt[i][j+1] + Jt[i][j] ) / 2;
	// perfect conducting ground to test
	if (groundmethod == 0) Et[nground][j] = 0;
      }
    }

    // SIBC boundary for Et
    if (groundmethod == 1) {
      for (int j = 0; j < hh-1; j++) {
	sumAi = 0;
	for (int m = 0; m < 6; m++) {
	  sibc_Ai[j][m] = -sibc_pi1[j][m] * Hp[nground][j] - sibc_pi2[j][m] * Hpold[nground][j] + sibc_pi3[j][m] * sibc_Aiold[j][m];
	  sumAi = sumAi + sibc_Ai[j][m];
	}
	Et[nground][j] = -sibc_eta[j] * Hp[nground][j] - sumAi;
      }
    }
    
    // pml correction
    if (dopml_top) {
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < hh-1; j++) {
	  Et[i+rshift][j] = Et[i+rshift][j] - c2e[i+rshift][j] / r[i+rshift] * P_etr[i][j];
	}
      }
    }

    // ----------------------------------------
#pragma omp parallel for
    // Ep update
    for (int i = 1; i < rr-1; i++) {
      for (int j = 1; j < hh-1; j++) {
	Ep[i][j] = c1e[i][j] * Ep[i][j] \
	  + c2e[i][j] / r[i] * ( (r[i]+dr[i]/2) * Ht[i][j] - (r[i]-dr[i-1]/2) * Ht[i-1][j] ) / ((dr[i]+dr[i-1])/2) \
	  - c2e[i][j] / r[i] * ( Hr[i][j] - Hr[i][j-1] ) / dth - c2e[i][j] * Jp[i][j];
	// perfect conducting ground
	if (groundmethod == 0) Ep[nground][j] = 0;
      }
    }
    // pml corrections
    if (dopml_top) {
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < hh-1; j++) {
	  Ep[i+rshift][j] = Ep[i+rshift][j] + c2e[i+rshift][j] / r[i+rshift] * P_epr[i][j];
	}
      }
    }
    if (dopml_wall) {
      for (int i = 0; i < rr-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  Ep[i][j+hshift] = Ep[i][j+hshift] - c2e[i][j+hshift] / r[i] * P_ept[i][j];
	}
      }
    }

    // output E and H probe points at desired points

    for (int i = 0; i < nprobes; i++) {
      ir = prober[i];
      it = probet[i];
      Erprobe[t][i] = 0.5*(Er[ir][it]+Er[ir-1][it]);
      Etprobe[t][i] = 0.5*(Et[ir][it]+Et[ir][it-1]);
      Epprobe[t][i] = Ep[ir][it];

      Hrprobe[t][i] = 0.5*(Hr[ir][it]+Hr[ir][it-1]);
      Htprobe[t][i] = 0.5*(Ht[ir][it]+Ht[ir-1][it]);
      Hpprobe[t][i] = 0.25*(Hp[ir][it]+Hp[ir][it-1]+Hp[ir-1][it]+Hp[ir-1][it-1]);
    }

    // update "old" fields
    for (int i = 0; i < rr; i++) {
      for (int j = 0; j < hh; j++) {
	Htold[i][j] = Ht[i][j];
	Hpold[i][j] = Hp[i][j];
      }
    }

    for (int j = 0; j < hh; j++) {
      for (int m = 0; m < 6; m++) {
	sibc_Aiold[j][m] = sibc_Ai[j][m];
      }
    }

    // ---------------------------------------
    // update Eeff, for use in ionization updates

    if (doionosphere) {
#pragma omp parallel for private(HpmL,HpmR,Epar,Eperp2,Sr,St,Sp,Erm,Etm,Hrm,Htm,Hpm)
      for (int i = nground+1; i < rr-1; i++) {

	// on axis
	Erm = ( dr[i-1] * Er[i][0] + dr[i] * Er[i-1][0] ) / (dr[i] + dr[i-1]);
	// based on fields on axis, all S (poynting) components should be zero. They are already initialized, so no need to update.
	Epar = Erm * wcer[0] / wce0[0];
	Emag[i][0] = Erm;
	Eperp2 = Emag[i][0]*Emag[i][0] - Epar*Epar;
	Eeff[i][0] = sqrt( Epar*Epar + Eperp2 * nue[i][0]*nue[i][0] / ( nue[i][0]*nue[i][0] + wce0[0]*wce0[0] ) );

	// everywhere else
	for (int j = 1; j < hh-1; j++) {
	  Erm = ( dr[i-1] * Er[i][j] + dr[i] * Er[i-1][j] ) / (dr[i] + dr[i-1]);
	  Etm = (Et[i][j] + Et[i][j-1])/2;
	  Hrm = (Hr[i][j+1] + Hr[i][j])/2;
	  Htm = ( dr[i-1] * Ht[i][j] + dr[i] * Ht[i-1][j] ) / (dr[i] + dr[i-1]);
	  HpmL = ( dr[i-1] * Hp[i][j-1] + dr[i] * Hp[i-1][j-1] ) / (dr[i] + dr[i-1]);
	  HpmR = ( dr[i-1] * Hp[i][j] + dr[i] * Hp[i-1][j] ) / (dr[i] + dr[i-1]);
	  Hpm = (HpmL + HpmR)/2;
	  Epar = (wcer[j]/wce0[j])*Erm + (wcet[j]/wce0[j])*Etm + (wcep[j]/wce0[j])*Ep[i][j];
	  Emag[i][j] = sqrt( Erm*Erm + Etm*Etm + Ep[i][j]*Ep[i][j] );
	  Eperp2 = Emag[i][j]*Emag[i][j] - Epar*Epar;
	  Eeff[i][j] = sqrt( Epar*Epar + Eperp2 * nue[i][j]*nue[i][j] / ( nue[i][j]*nue[i][j] + wce0[j]*wce0[j] ) );


	  Sr = Etm*Hpm - Ep[i][j]*Htm;   // this is not quite right, because H is at half-time-steps...
	  St = Ep[i][j]*Hrm - Erm*Hpm;
	  Sp = Erm*Htm - Etm*Hrm;
	  S[i][j] = sqrt( Sr*Sr + St*St + Sp*Sp );

	}
      }

      // test to find numerical instability!
      if (0) {
	for (int i = nground+1; i < rr-1; i++) {
	  for (int j = 0; j < hh-1; j++) {
	    if (Eeff[i][j] > 100*Ek[i][j]) {
	      logfile.open("log.txt",std::fstream::app);
	      logfile << "Unstable! Ending main loop after iteration " << t << "...\n";
	      logfile.close();
	      goto endofbigloop;
	    }
	  }
	}
      }

    } // if doionosphere
    

    // Ionization updates.
    // ------------------------------------------
    // update mue and nue, calculate vi, va, ne. Not going to do it inside PML! Makes things unstable!

    if (doionosphere & doioniz & !dotransmitter) {

#pragma omp parallel for private(EoEk,vN21P,vN22P,vN2P1N,vN2PM,vO2P1N,vOred,vOgrn,tauN21P,tauN22P,tauN2P1N,tauN2PM,tauO2P1N,tauOred,tauOgrn,Aion,nenew,topE,botE,vi,va,vd,mue)
      for (int i = nonlinearstart; i < rr-pmllen+1; i++) {
	for (int j = 0; j < hh-pmllen+1; j++) {

	  topE = 1;
	  botE = 0;
	  // find index into Efield rate array
	  for (int m = 0; m < nume; m++) {
	    if (efield[m] > Eeff[i][j] / nd2[i][j]) {
	      topE = m;
	      botE = topE - 1;
	      break;
	    }
	  }

	  // read off mue, vi, va, vN21P, vN22P from the arrays and interpolate
	  
	  if (topE <= 1) {
	    // set everything to first value 
	    mue = muArray[i][0] / nd2[i][j];
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
	    mue = (1/nd2[i][j]) * ( muArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((muArray[i][topE]-muArray[i][botE])/(efield[topE]-efield[botE])) );
	    vi = nd2[i][j] * ( viArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((viArray[i][topE]-viArray[i][botE])/(efield[topE]-efield[botE])) );
	    va = nd2[i][j] * ( vaArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((vaArray[i][topE]-vaArray[i][botE])/(efield[topE]-efield[botE])) );
	    vN21P = nd2[i][j] * ( N21pArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((N21pArray[i][topE]-N21pArray[i][botE])/(efield[topE]-efield[botE])) );
	    vN22P = nd2[i][j] * ( N22pArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((N22pArray[i][topE]-N22pArray[i][botE])/(efield[topE]-efield[botE])) );
	    vN2P1N = nd2[i][j] * ( N2p1nArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((N2p1nArray[i][topE]-N2p1nArray[i][botE])/(efield[topE]-efield[botE])) );
	    vN2PM = nd2[i][j] * ( N2pMArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((N2pMArray[i][topE]-N2pMArray[i][botE])/(efield[topE]-efield[botE])) );
	    vO2P1N = nd2[i][j] * ( O2p1nArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((O2p1nArray[i][topE]-O2p1nArray[i][botE])/(efield[topE]-efield[botE])) );
	    vOred = nd2[i][j] * ( OrArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((OrArray[i][topE]-OrArray[i][botE])/(efield[topE]-efield[botE])) );
	    vOgrn = nd2[i][j] * ( OgArray[i][botE] + (Eeff[i][j]/nd2[i][j] - efield[botE])*((OgArray[i][topE]-OgArray[i][botE])/(efield[topE]-efield[botE])) );
	  }

	  // old method for detachment coefficient
	  EoEk = Eeff[i][j] / Ek[i][j];
	  if (EoEk < 0.05) {
	    vd = 0;
	  } else {
	    vd = 0.78 * nd2[i][j] * 1.08e-18 * sqrt(EoEk) * exp(-0.078375 / EoEk);
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
	      nenew = ( Aion[0] * ne[i][j] + Aion[1] * nOm[i][j] ) / Aion[4];     // temporary value so we don't overwrite
	      nOm[i][j] = ( Aion[2] * ne[i][j] + Aion[3] * nOm[i][j] ) / Aion[4];
	      ne[i][j] = nenew;
	    } else {
	      ne[i][j] = ne[i][j] / (1 - (vi - va)*dt);
	      nOm[i][j] = nOm[i][j] + va*dt*ne[i][j];
	    }
	  } else {
	    // throttle the ionization so it doesn't go crazy unstable
	    ne[i][j] = ne[i][j] * exp(0.5);
	  }

	  wpe[i][j] = QE * sqrt( ne[i][j] / (ME * E0) );
	  nue[i][j]  = (QE / ME) / mue;

	  // assume that ionization creates a positive ion, and attachment creates a negative ion?
	  // for now, we will assume ion densities don't change.

	  // optical: Need to input relative densities of O2, N2, O versus altitude to make these equations valid at higher altitudes.

	  tauN22P = 1/(2e7 + 3e-19*(0.21*nd2[i][j]) );  // quenching through O2, which is 21% of density
	  tauN21P = 1/(1.7e5 + 1e-17*(0.78*nd2[i][j]) ); // quenching through N2, which is 78% of density
	  tauN2P1N = 1/(1.4e7 + 4e-16*(0.99*nd2[i][j]) ); // quenched by N2 and O2
	  tauN2PM = 1/(7e4 + 5e-16*(0.78*nd2[i][j]) ); // quenched by N2
	  tauO2P1N = 1/(8.3e5 + 2e-16*(0.78*nd2[i][j]) ); // quenched by N2
	  tauOred = 1/(0.0091 + 5e-17*(0.78*nd2[i][j]) );
	  tauOgrn = 1/(1.43 + 3e-19*(0.21*nd2[i][j]) );

	  // attemping backward Euler (implicit) for optics, since time constants are smaller than dt in some cases.

	  nN22P[i][j] = tauN22P/(dt + tauN22P) * nN22P[i][j] + (tauN22P*dt)/(dt + tauN22P) * ( vN22P * ne[i][j] );
	  nN21P[i][j] = tauN21P/(dt + tauN21P) * nN21P[i][j] + (tauN21P*dt)/(dt + tauN21P) * ( vN21P * ne[i][j] + 2e7 * nN22P[i][j] ); // cascading
	  nN2P1N[i][j] = tauN2P1N/(dt + tauN2P1N) * nN2P1N[i][j] + (tauN2P1N*dt)/(dt + tauN2P1N) * ( vN2P1N * ne[i][j] );
	  nN2PM[i][j] = tauN2PM/(dt + tauN2PM) * nN2PM[i][j] + (tauN2PM*dt)/(dt + tauN2PM) * ( vN2PM * ne[i][j] );
	  nO2P1N[i][j] = tauO2P1N/(dt + tauO2P1N) * nO2P1N[i][j] + (tauO2P1N*dt)/(dt + tauO2P1N) * ( vO2P1N * ne[i][j] );

	  nOgrn[i][j] = tauOgrn/(dt + tauOgrn) * nOgrn[i][j] + (tauOgrn*dt)/(dt + tauOgrn) * ( vOgrn * ne[i][j] );
	  nOred[i][j] = tauOred/(dt + tauOred) * nOred[i][j] + (tauOred*dt)/(dt + tauOred) * ( vOred * ne[i][j] + 1.43 * nOgrn[i][j] );

	}
      }

    } // if doioniz

    // ------------------------------------------
    // write optical outputs to elve cube
    
    if (doelve) {
#pragma omp parallel for private(elvetimeind)
      for (int m = 0; m < numpixels; m++) {
	for (int k = stepi; k < rr-pmllen; k++) {
	  if ( (thind[m][k] > pmllen) && (thind[m][k] <= (hh-pmllen) ) ) {
	    
	    elvetimeind = round((tdelay[m][k] + t)*dt/elvedt);
	    if (elvetimeind <= elvesteps) {
	      // the constant below (e.g. 1.7e-5) is the Einstein coefficient (1.7e5) * 1e-6 (for Rayleigh definition) * 1e-4 (convert from m^2 to cm^2)
	      elveN21P[m][elvetimeind] = elveN21P[m][elvetimeind] + 1.7e-5 * nN21P[k][thind[m][k]] * raylength[m][k] * dt;
	      elveN22P[m][elvetimeind] = elveN22P[m][elvetimeind] + 2.0e-3 * nN22P[k][thind[m][k]] * raylength[m][k] * dt;
	      elveN2P1N[m][elvetimeind] = elveN2P1N[m][elvetimeind] + 1.4e-3 * nN2P1N[k][thind[m][k]] * raylength[m][k] * dt;
	      elveN2PM[m][elvetimeind] = elveN2PM[m][elvetimeind] + 7.1e-6 * nN2PM[k][thind[m][k]] * raylength[m][k] * dt;
	      elveO2P1N[m][elvetimeind] = elveO2P1N[m][elvetimeind] + 8.3e-5 * nO2P1N[k][thind[m][k]] * raylength[m][k] * dt;
	      elveOred[m][elvetimeind] = elveOred[m][elvetimeind] + 0.0091e-10 * nOred[k][thind[m][k]] * raylength[m][k] * dt;
	      elveOgrn[m][elvetimeind] = elveOgrn[m][elvetimeind] + 1.43e-10 * nOgrn[k][thind[m][k]] * raylength[m][k] * dt;
	      // since we multiply by the effective dt, what we are measuring here is the number of photons, in R-s, collected at each pixel in each elvetimeind.
	    }
	    
	  }
	}
      }

    } // if doelve

    // ------------------------------------------
    // assuming ions densities don't really change, or that the changes are not sufficient to affect currents etc. 
    // I may need to correct this assumption - perhaps Dni = -Dne ? 

    // ----------------------------------------
    // J update equations... parallelizing these actually slows it down, since there are so many private variables
    
    if (doionosphere) {
#pragma omp parallel for private(Ee1,Ee2,Se1,Ce1,Ce2,Ce3,Ce4,Ei1,Ei2,Si1,Ci1,Ci2,Ci3,Ci4,Ae,Ke,Ai,Ki,Emidr,Jer0,Jet0,Jir0,Jit0)
      for (int i = nground+1; i < rr; i++) {
	for (int j = 0; j < hh-1; j++) {

	  Ee1 = exp(-nue[i][j]*dt);
	  Ee2 = 1/(wce0[j]*wce0[j] + nue[i][j]*nue[i][j]);
	  if (wce0[j] == 0) {
	    Se1 = 1;
	    Ce1 = 0.5;  // convergence value
	  } else {
	    Se1 = sin(wce0[j]*dt)/wce0[j];
	    Ce1 = (1 - cos(wce0[j]*dt))/(wce0[j]*wce0[j]);
	  }
	  Ce2 = (1 - Ee1)/nue[i][j] - Ee1*nue[i][j]*Ce1 - Ee1*Se1;
	  Ce3 = nue[i][j] * (1 - Ee1*cos(wce0[j]*dt)) + Ee1*wce0[j]*sin(wce0[j]*dt);
	  Ce4 = 1 - Ee1*cos(wce0[j]*dt) - Ee1*nue[i][j]*Se1;

	  // same for ions
	  Ei1 = exp(-nui[i][j]*dt);
	  Ei2 = 1/(wci0[j]*wci0[j] + nui[i][j]*nui[i][j]);
	  if (wci0[j] == 0) {
	    Si1 = 1;
	    Ci1 = 0.5;
	  } else {
	    Si1 = sin(wci0[j]*dt)/wci0[j];
	    Ci1 = (1 - cos(wci0[j]*dt))/(wci0[j]*wci0[j]);
	  }
	  Ci2 = (1 - Ei1)/nui[i][j] - Ei1*nui[i][j]*Ci1 - Ei1*Si1;
	  Ci3 = nui[i][j] * (1 - Ei1*cos(wci0[j]*dt)) + Ei1*wci0[j]*sin(wci0[j]*dt);
	  Ci4 = 1 - Ei1*cos(wci0[j]*dt) - Ei1*nui[i][j]*Si1;

	  // A and K matrices
	  Ae[0][0] = Ee1 * ( Ce1*wcer[j]*wcer[j] + cos(wce0[j]*dt) );
	  Ae[0][1] = Ee1 * ( Ce1*wcer[j]*wcet[j] - Se1*wcep[j] );
	  Ae[0][2] = Ee1 * ( Ce1*wcer[j]*wcep[j] + Se1*wcet[j] );
	  Ae[1][0] = Ee1 * ( Ce1*wcet[j]*wcer[j] + Se1*wcep[j] );     
	  Ae[1][1] = Ee1 * ( Ce1*wcet[j]*wcet[j] + cos(wce0[j]*dt) );
	  Ae[1][2] = Ee1 * ( Ce1*wcet[j]*wcep[j] - Se1*wcer[j] );
	  Ae[2][0] = Ee1 * ( Ce1*wcep[j]*wcer[j] - Se1*wcet[j] );
	  Ae[2][1] = Ee1 * ( Ce1*wcep[j]*wcet[j] + Se1*wcer[j] );
	  Ae[2][2] = Ee1 * ( Ce1*wcep[j]*wcep[j] + cos(wce0[j]*dt) );
            
	  Ke[0][0] = Ee2 * ( Ce2*wcer[j]*wcer[j] + Ce3 );
	  Ke[0][1] = Ee2 * ( Ce2*wcer[j]*wcet[j] - Ce4*wcep[j] );
	  Ke[0][2] = Ee2 * ( Ce2*wcer[j]*wcep[j] + Ce4*wcet[j] );
	  Ke[1][0] = Ee2 * ( Ce2*wcet[j]*wcer[j] + Ce4*wcep[j] );
	  Ke[1][1] = Ee2 * ( Ce2*wcet[j]*wcet[j] + Ce3 );
	  Ke[1][2] = Ee2 * ( Ce2*wcet[j]*wcep[j] - Ce4*wcer[j] );
	  Ke[2][0] = Ee2 * ( Ce2*wcep[j]*wcer[j] - Ce4*wcet[j] );
	  Ke[2][1] = Ee2 * ( Ce2*wcep[j]*wcet[j] + Ce4*wcer[j] );
	  Ke[2][2] = Ee2 * ( Ce2*wcep[j]*wcep[j] + Ce3 );

	  // A and K for ions
	  Ai[0][0] = Ei1 * ( Ci1*wcir[j]*wcir[j] + cos(wci0[j]*dt) );
	  Ai[0][1] = Ei1 * ( Ci1*wcir[j]*wcit[j] - Si1*wcip[j] );
	  Ai[0][2] = Ei1 * ( Ci1*wcir[j]*wcip[j] + Si1*wcit[j] );
	  Ai[1][0] = Ei1 * ( Ci1*wcit[j]*wcir[j] + Si1*wcip[j] );     
	  Ai[1][1] = Ei1 * ( Ci1*wcit[j]*wcit[j] + cos(wci0[j]*dt) );
	  Ai[1][2] = Ei1 * ( Ci1*wcit[j]*wcip[j] - Si1*wcir[j] );
	  Ai[2][0] = Ei1 * ( Ci1*wcip[j]*wcir[j] - Si1*wcit[j] );
	  Ai[2][1] = Ei1 * ( Ci1*wcip[j]*wcit[j] + Si1*wcir[j] );
	  Ai[2][2] = Ei1 * ( Ci1*wcip[j]*wcip[j] + cos(wci0[j]*dt) );
            
	  Ki[0][0] = Ei2 * ( Ci2*wcir[j]*wcir[j] + Ci3 );
	  Ki[0][1] = Ei2 * ( Ci2*wcir[j]*wcit[j] - Ci4*wcip[j] );
	  Ki[0][2] = Ei2 * ( Ci2*wcir[j]*wcip[j] + Ci4*wcit[j] );
	  Ki[1][0] = Ei2 * ( Ci2*wcit[j]*wcir[j] + Ci4*wcip[j] );
	  Ki[1][1] = Ei2 * ( Ci2*wcit[j]*wcit[j] + Ci3 );
	  Ki[1][2] = Ei2 * ( Ci2*wcit[j]*wcip[j] - Ci4*wcir[j] );
	  Ki[2][0] = Ei2 * ( Ci2*wcip[j]*wcir[j] - Ci4*wcit[j] );
	  Ki[2][1] = Ei2 * ( Ci2*wcip[j]*wcit[j] + Ci4*wcir[j] );
	  Ki[2][2] = Ei2 * ( Ci2*wcip[j]*wcip[j] + Ci3 );

	  if (i == rr-1) {
	    Emidr = Er[i-1][j];
	  } else {
	    Emidr = ( dr[i-1] * Er[i][j] + dr[i] * Er[i-1][j] ) / (dr[i] + dr[i-1]);
	  }

	  // okay, updates for J finally.
	  // on axis: only Er is non-zero
	  if (j == 0) {

	    Jer[i][j] = Ae[0][0] * Jer[i][j] + E0 * wpe[i][0]*wpe[i][0] * Ke[0][0] * Emidr;
	    Jir[i][j] = Ai[0][0] * Jir[i][j] + E0 * wpi[i][0]*wpi[i][0] * Ki[0][0] * Emidr;

	  } else {

	    Jer0 = Ae[0][0] * Jer[i][j] + Ae[0][1] * Jet[i][j] + Ae[0][2] * Jep[i][j] \
	      + E0 * wpe[i][j]*wpe[i][j] * ( Ke[0][0] * Emidr + Ke[0][1] * ( Et[i][j] + Et[i][j-1] )/2 + Ke[0][2] * Ep[i][j] );

	    Jet0 = Ae[1][0] * Jer[i][j] + Ae[1][1] * Jet[i][j] + Ae[1][2] * Jep[i][j] \
	      + E0 * wpe[i][j]*wpe[i][j] * ( Ke[1][0] * Emidr + Ke[1][1] * ( Et[i][j] + Et[i][j-1] )/2 + Ke[1][2] * Ep[i][j] );

	    Jep[i][j] = Ae[2][0] * Jer[i][j] + Ae[2][1] * Jet[i][j] + Ae[2][2] * Jep[i][j] \
	      + E0 * wpe[i][j]*wpe[i][j] * ( Ke[2][0] * Emidr + Ke[2][1] * ( Et[i][j] + Et[i][j-1] )/2 + Ke[2][2] * Ep[i][j] );

	    Jir0 = Ai[0][0] * Jir[i][j] + Ai[0][1] * Jit[i][j] + Ai[0][2] * Jip[i][j] \
	      + E0 * wpi[i][j]*wpi[i][j] * ( Ki[0][0] * Emidr + Ki[0][1] * ( Et[i][j] + Et[i][j-1] )/2 + Ki[0][2] * Ep[i][j] );
	  
	    Jit0 = Ai[1][0] * Jir[i][j] + Ai[1][1] * Jit[i][j] + Ai[1][2] * Jip[i][j] \
	      + E0 * wpi[i][j]*wpi[i][j] * ( Ki[1][0] * Emidr + Ki[1][1] * ( Et[i][j] + Et[i][j-1] )/2 + Ki[1][2] * Ep[i][j] );
	
	    Jip[i][j] = Ai[2][0] * Jir[i][j] + Ai[2][1] * Jit[i][j] + Ai[2][2] * Jip[i][j] \
	      + E0 * wpi[i][j]*wpi[i][j] * ( Ki[2][0] * Emidr + Ki[2][1] * ( Et[i][j] + Et[i][j-1] )/2 + Ki[2][2] * Ep[i][j] );

	    // these temporary values are necessary so that the fields aren't updating twice!
	    Jer[i][j] = Jer0;
	    Jet[i][j] = Jet0;
	    Jir[i][j] = Jir0;
	    Jit[i][j] = Jit0;
	  }

	}
      }

      // add currents together

      for (int i = nground+1; i < rr-1; i++) {
	for (int j = 0; j < hh-1; j++) {
	  Jr[i][j] = Jer[i][j] + Jir[i][j];
	  Jt[i][j] = Jet[i][j] + Jit[i][j];
	  Jp[i][j] = Jep[i][j] + Jip[i][j];
	}
      }

      // do heating: on axis
      for (int i = 1; i < rr-1; i++) {
	// on axis (j = 0)
	JdotE[i][0] = (Jr[i][0] * (Er[i][0] + Er[i-1][0]) / 2);
	heat[i][0] = heat[i][0] + dt * JdotE[i][0];
	for (int j = 1; j < hh-1; j++) {
	  JdotE[i][j] = Jr[i][j] * (dr[i-1] * Er[i][j] + dr[i] * Er[i-1][j])/(dr[i]+dr[i-1]) + Jt[i][j] * (Et[i][j] + Et[i][j-1])/2 + Jp[i][j] * Ep[i][j];
	  heat[i][j] = heat[i][j] + dt * JdotE[i][j];
	}
      }

    } // if doionosphere


    // VLF Transmitter induced heating, using J.E heating and cooling from Rodriquez (1994) (updated). 
    // ---------------------------------------------------------------------

    // override!
    if (0) {
      //if (dotransmitter) {

      // first index for heating: start at 40 km. Below causes instability...
      // int nonlinearstart = round(60e3/dr1 + nground);

#pragma omp parallel for private(f_N2,f_O2,gg,Lelast_N2,Lelast_O2,Lrot_N2,Lrot_O2,Lvib_N2,Lvib_O2,Le)
      for (int i = nonlinearstart; i < rr-pmllen+1; i++) {
	for (int j = 0; j < hh-pmllen+1; j++) {

	  // cooling
	  f_N2 = 1.06e4 + 7.51e3*tanh(0.0011*(Te[i][j]-1800.0));
	  f_O2 = 3300.0 - 839.0*sin(0.000191*(Te[i][j]-2700.0));
	  gg = 3300.0 + 1.233*(Te[i][j]-1000) - (2.056e-4)*(Te[i][j]-1000.0)*(Te[i][j]-4000.0);
	  Lelast_N2 = (1.89e-44)*ne[i][j]*(0.78*nd2[i][j])*(1.0 - 1.21e-4*Te[i][j])*Te[i][j]*(Te[i][j]-Te0[i][j]);
	  Lelast_O2 = (1.29e-43)*ne[i][j]*(0.21*nd2[i][j])*(1.0 + 3.6e-2*sqrt(Te[i][j]))*sqrt(Te[i][j])*(Te[i][j]-Te0[i][j]);
	  Lrot_N2 = (4.65e-39)*ne[i][j]*(0.78*nd2[i][j])*(Te[i][j]-Te0[i][j])/sqrt(Te[i][j]);
	  Lrot_O2 = (1.11e-38)*ne[i][j]*(0.21*nd2[i][j])*(Te[i][j]-Te0[i][j])/sqrt(Te[i][j]);
	  Lvib_N2 = (4.79e-37)*ne[i][j]*(0.78*nd2[i][j])*exp(f_N2*(Te[i][j]-2000.0)/(2000.0*Te[i][j])) * (1-exp(-gg*(Te[i][j]-Te0[i][j])/(Te[i][j]*Te0[i][j])));
	  Lvib_O2 = (8.32e-38)*ne[i][j]*(0.21*nd2[i][j])*exp(f_O2*(Te[i][j]-700.0)/(700.0*Te[i][j])) * (1-exp(-2700.0*(Te[i][j]-Te0[i][j])/(Te[i][j]*Te0[i][j])));
	  Le = Lelast_N2 + Lrot_N2 + Lvib_N2+ Lelast_O2 + Lrot_O2 + Lvib_O2;

	  Te[i][j] = Te[i][j] + dt * (2.0/(3.0*ne[i][j]*KB)) * (JdotE[i][j] - Le);

	  // update collision frequency
	  nue[i][j] = 1.6 * ( 2.33e-17*(0.78*nd2[i][j])*(1.0 - 1.25e-4*Te[i][j])*Te[i][j] + 1.82e-16*(0.21*nd2[i][j])*(1.0 + 3.6e-2*sqrt(Te[i][j]))*sqrt(Te[i][j]) );

	}
      }
    } // if do transmitter



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

    // compute DFTs along the ground at specified frequencies
    if (doDFT) {

      for (int j = 0; j < hh; j++) {
	  for (int m = 0; m < numDFTfreqs; m++) {

	    DFTfieldsEr[j][2*m] += Er[nground+1][j] * sin(2*PI*DFTfreqs[m]*(t+1/2)*dt) * dt;
            DFTfieldsEr[j][2*m+1] += Er[nground+1][j] * cos(2*PI*DFTfreqs[m]*(t+1/2)*dt) * dt;
	    DFTfieldsEt[j][2*m] += Et[nground+1][j] * sin(2*PI*DFTfreqs[m]*(t+1/2)*dt) * dt;
            DFTfieldsEt[j][2*m+1] += Et[nground+1][j] * cos(2*PI*DFTfreqs[m]*(t+1/2)*dt) * dt;
	    DFTfieldsEp[j][2*m] += Ep[nground+1][j] * sin(2*PI*DFTfreqs[m]*(t+1/2)*dt) * dt;
            DFTfieldsEp[j][2*m+1] += Ep[nground+1][j] * cos(2*PI*DFTfreqs[m]*(t+1/2)*dt) * dt;
	    
	    DFTfieldsHr[j][2*m] += Hr[nground+1][j] * sin(2*PI*DFTfreqs[m]*t*dt) * dt;
	    DFTfieldsHr[j][2*m+1] += Hr[nground+1][j] * cos(2*PI*DFTfreqs[m]*t*dt) * dt;
	    DFTfieldsHt[j][2*m] += Ht[nground+1][j] * sin(2*PI*DFTfreqs[m]*t*dt) * dt;
            DFTfieldsHt[j][2*m+1] += Ht[nground+1][j] * cos(2*PI*DFTfreqs[m]*t*dt) * dt;
	    DFTfieldsHp[j][2*m] += Hp[nground+1][j] * sin(2*PI*DFTfreqs[m]*t*dt) * dt;
            DFTfieldsHp[j][2*m+1] += Hp[nground+1][j] * cos(2*PI*DFTfreqs[m]*t*dt) * dt;
	}
      }

    }
    

    ////////////////
    
    if (t % (int)(tsteps/numfiles) == 0) {
      partialtime = 100 * t/tsteps;   
      logfile.open("log.txt",std::fstream::app); 
      logfile << "t = " << t << ": " << partialtime << " % Done...  ";

      DATATYPE Emax = 0;
      for (int i = 0; i < rr; i++) {
	for (int j = 0; j < hh; j++) {
	  if (fabs(Er[i][j]) > Emax) {
	    Emax = fabs(Er[i][j]);
	  }
	}
      }
      logfile << "Maximum abs(Er) is " << Emax << "\n";
      logfile.close();

      // decimate all of these fields!

      for (int i = 0; i < drr; i++) {
	for (int j = 0; j < dhh; j++) {
	  Erdec[i][j] = Er[i*decfactor][j*decfactor];
	  Etdec[i][j] = Et[i*decfactor][j*decfactor];
	  Epdec[i][j] = Ep[i*decfactor][j*decfactor];
	  Hrdec[i][j] = Hr[i*decfactor][j*decfactor];
	  Htdec[i][j] = Ht[i*decfactor][j*decfactor];
	  Hpdec[i][j] = Hp[i*decfactor][j*decfactor];
	  Jerdec[i][j] = Jer[i*decfactor][j*decfactor];
	  Jetdec[i][j] = Jet[i*decfactor][j*decfactor];
	  Jepdec[i][j] = Jep[i*decfactor][j*decfactor];
	  Jirdec[i][j] = Jir[i*decfactor][j*decfactor];
	  Jitdec[i][j] = Jit[i*decfactor][j*decfactor];
	  Jipdec[i][j] = Jip[i*decfactor][j*decfactor];

	  Eeffdec[i][j] = Eeff[i*decfactor][j*decfactor];
	  Ekdec[i][j] = Ek[i*decfactor][j*decfactor];
	  heatdec[i][j] = heat[i*decfactor][j*decfactor];
	  Sdec[i][j] = S[i*decfactor][j*decfactor];
	  Tedec[i][j] = Te[i*decfactor][j*decfactor];
	  nedec[i][j] = ne[i*decfactor][j*decfactor];
	  nOmdec[i][j] = nOm[i*decfactor][j*decfactor];
	  nuedec[i][j] = nue[i*decfactor][j*decfactor];
	  nN21Pdec[i][j] = nN21P[i*decfactor][j*decfactor];
	  nN22Pdec[i][j] = nN22P[i*decfactor][j*decfactor];
	  nN2P1Ndec[i][j] = nN2P1N[i*decfactor][j*decfactor];
	  nN2PMdec[i][j] = nN2PM[i*decfactor][j*decfactor];
	  nO2P1Ndec[i][j] = nO2P1N[i*decfactor][j*decfactor];
	  nOreddec[i][j] = nOred[i*decfactor][j*decfactor];
	  nOgrndec[i][j] = nOgrn[i*decfactor][j*decfactor];
	}
      }

      // write to files

      FILE * filePtr;
 
      if (savefields[0]) {
	filePtr = fopen("output_E.dat","ab");
	fwrite(Erdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Etdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Epdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fclose(filePtr);
      }

      if (savefields[1]) {
	filePtr = fopen("output_J.dat","ab");
	fwrite(Jerdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Jetdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Jepdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Jirdec,sizeof(DATATYPE),drr*dhh,filePtr);      
	fwrite(Jitdec,sizeof(DATATYPE),drr*dhh,filePtr);      
	fwrite(Jipdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fclose(filePtr);
      }

      if (savefields[2]) {
	filePtr = fopen("output_H.dat","ab");
	fwrite(Hrdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Htdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Hpdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fclose(filePtr);
      }

      if (savefields[3]) {
	filePtr = fopen("output_K.dat","ab");
	fwrite(Eeffdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Ekdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(heatdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Sdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(Tedec,sizeof(DATATYPE),drr*dhh,filePtr);
	fclose(filePtr);
      }
      
      if (savefields[4]) {
	filePtr = fopen("output_D.dat","ab");
	fwrite(nedec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(nOmdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(nuedec,sizeof(DATATYPE),drr*dhh,filePtr);
	fclose(filePtr);
      }

      if (savefields[5]) {
	filePtr = fopen("output_O.dat","ab");
	fwrite(nN21Pdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(nN22Pdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(nN2P1Ndec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(nN2PMdec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(nO2P1Ndec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(nOreddec,sizeof(DATATYPE),drr*dhh,filePtr);
	fwrite(nOgrndec,sizeof(DATATYPE),drr*dhh,filePtr);
	fclose(filePtr);
      }

    } // if tstep to write out
    
    // end giant time loop
  }

 endofbigloop:

  /////////////////////////////////////////////////////////////////////////////

  // write the elve

  if (doDFT) {
    FILE * dftFile;
    dftFile = fopen("dft.dat","wb");
    fwrite(&numDFTfreqs,sizeof(int),1,dftFile);
    fwrite(&DFTfreqs,sizeof(double),numDFTfreqs,dftFile);
    fwrite(&DFTfieldsEr,sizeof(DATATYPE),hh*numDFTfreqs*2,dftFile);
    fwrite(&DFTfieldsEt,sizeof(DATATYPE),hh*numDFTfreqs*2,dftFile);
    fwrite(&DFTfieldsEp,sizeof(DATATYPE),hh*numDFTfreqs*2,dftFile);
    fwrite(&DFTfieldsHr,sizeof(DATATYPE),hh*numDFTfreqs*2,dftFile);
    fwrite(&DFTfieldsHt,sizeof(DATATYPE),hh*numDFTfreqs*2,dftFile);
    fwrite(&DFTfieldsHp,sizeof(DATATYPE),hh*numDFTfreqs*2,dftFile);
    fclose(dftFile);
  }

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
    fwrite(&elveOred,sizeof(DATATYPE),numpixels*elvesteps,elveFile);
    fwrite(&elveOgrn,sizeof(DATATYPE),numpixels*elvesteps,elveFile);
    fclose(elveFile);

  }

  FILE * ProbeFile;
  ProbeFile = fopen("Probe.dat","wb");
  fwrite(&nprobes,sizeof(int),1,ProbeFile);
  fwrite(&prober,sizeof(int),nprobes,ProbeFile);
  fwrite(&probet,sizeof(int),nprobes,ProbeFile);
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

  logfile << "Total run time = " << totalruntime << " minutes.\n";
  logfile.close();

  // end program - do not type below this!

}

//////////////////////////////////////////////////////////

