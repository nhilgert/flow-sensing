#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mersenne_twister.hpp"

using namespace std;

#define vnaught 10.0 // equivalent to peclet
#define beta 10.0
// #define T 100.0
#define R 10.0
#define alpha 0.75
#define s 0.01
#define Nr 100
#define Nc 100

int main(int argc, char *argv[]) {
  printf("test1\n");
  int seed = atoi(argv[1]);
  double T = atof(argv[2]);
  MTRand rg(seed);

  vector<double> posx, posy, posz;

  // slice up R^3 space into wedges by dr, dcostheta
  double dr = (R - 1.0)/Nr;
  double dc = 2.0 / Nc;
  int record = (int) T; // last seconds recorded
  // initialize array for binning
  int count [record][(int) Nc];
  for (int k = 0; k < record; k++) {
      for (int j = 0; j< Nc; j++) {
          count[k][j] = 0;
      }
  }

  int tsteps = (int) T;
  double deltat = 0.1 * fmin(dr*dr,fmin(1.0/vnaught, fmin(s/alpha, 1/(4*M_PI*beta))));
  // make sure deltat isn't zero
  if (abs(deltat - 0.0) < 0.001) {
      deltat = 0.001;
  }
  deltat = 0.001;
  // probability of death event
  double pa = alpha / s * deltat;
  double conc [record][(int) Nc];
  int captured [(int) Nc];
  for (int j = 0; j < (int) Nc; j++) {
      captured[j] = 0;
    for (int i = 0; i < record; i++) {
      conc[i][j] = 0;
    }
  }
  
  std::string cons1 = "constants";
  std::string dist1 = "distribution";
  std::string capt1 = "captured";
  std::string idx = argv[1];
  std::string time = argv[2];
  std::string txt = ".txt";
  std::string cons = cons1 + idx + "t" + time + txt;
  std::string dist = dist1 + idx + "t" + time + txt;
  std::string capt = capt1 + idx + "t" + time + txt;
  char consstring[cons.size() + 1];
  char diststring[dist.size() + 1];
  char captstring[capt.size() + 1];
 
  strcpy(consstring, cons.c_str()); 
  strcpy(diststring, dist.c_str());
  strcpy(captstring, capt.c_str());

  FILE *fpc;
  fpc = fopen(captstring,"w");

  for (int i = 0; i < (int) (2.0 * M_PI * beta/(1.0 + alpha) * (R*R - 1.0)); i++) {
    double ctheta = 2.0 * rg.rand() - 1;
    double phi = 2.0 * M_PI * rg.randExc();
    double r = sqrt(1.0 + rg.rand() * (R*R - 1.0));
    posx.push_back(r * sqrt(1.0 - ctheta * ctheta) * cos(phi));
    posy.push_back(r * sqrt(1.0 - ctheta * ctheta) * sin(phi));
    posz.push_back(r * ctheta);
  }
  // main loop
  int tidx = 0;
  for (double t = 0.0; t < (double) tsteps; t += deltat) {
    if (rg.rand() <= 4*M_PI*beta*deltat) { // birth event
       double ctheta = 2.0 * rg.rand() - 1;
       double phi = 2.0 * M_PI * rg.randExc();
        posx.push_back(sqrt(1.0 - ctheta * ctheta) * cos(phi));
        posy.push_back(sqrt(1.0 - ctheta * ctheta) * sin(phi));
        posz.push_back(ctheta);
      }
    for (int i = 0; i < posx.size(); i++) {
      // diffusion
      double dx = rg.randNorm(0.0,sqrt(2.0 * deltat));
      double dy = rg.randNorm(0.0,sqrt(2.0 * deltat));
      double dz = rg.randNorm(0.0,sqrt(2.0 * deltat));

      // drift; coordinates
      double r = sqrt(posx.at(i)*posx.at(i) + posy.at(i)*posy.at(i) + posz.at(i)*posz.at(i));
      double ctheta = posz.at(i)/r;
      double stheta = sqrt(1.0 - ctheta * ctheta);
      double cphi = 0.0;
      double sphi = 0.0;
      // watch out for sign of sine
      if (stheta > 0.0) {
          cphi = posx.at(i)/(r*stheta);
          sphi = posy.at(i)/(r*stheta);
        }
      double zeta = 5 * pow(10,-4);
      double bigz = 1 + 3*zeta + 3*zeta*zeta;
      double ur = 1 - bigz/(r*r*r) + (3*zeta)/(r*r)*(1 + (zeta/r))*exp((1 - r)/zeta);
      double utheta = 1 + bigz/(2*r*r*r) - 3/(2 * r)*(1 + zeta/r + (zeta*zeta)/(r*r))*exp((1 - r)/zeta);
      // velocity components
      double vx = vnaught * ctheta * stheta * cphi * (ur - utheta);
      double vy = vnaught * ctheta * stheta * sphi * (ur - utheta);
      double vz = vnaught * (ctheta*ctheta*ur + stheta*stheta*utheta);
      posx.at(i) += dx + vx*deltat;
      posy.at(i) += dy + vy*deltat;
      posz.at(i) += dz + vz*deltat;
    }

      int i = 0;

      while (i < posx.size()) {
        double r = sqrt(posx.at(i)*posx.at(i) + posy.at(i)*posy.at(i) + posz.at(i)*posz.at(i));

	// deletion
	if (r >= R) {
          posx.erase(posx.begin() + i);
          posy.erase(posy.begin() + i);
          posz.erase(posz.begin() + i);
	}

        // absorption
        
        else if ((r < 1.0 + s && rg.rand() <= pa)) {
            double ctheta = posz.at(i)/r;
            posx.erase(posx.begin() + i);
            posy.erase(posy.begin() + i);
            posz.erase(posz.begin() + i);
            if (ctheta == 1.0) {
                ctheta -= dc/2.0;
            }
	    captured[(int) ((ctheta + 1.0)/dc)]++;
          } else if (r < 1.0) { // reflect
            posx.at(i) /= r*r;
            posy.at(i) /= r*r;
            posz.at(i) /= r*r;
	    i++;
	  } else { i++; }

        } } 
  for (int c = 0; c < (int) Nc; c++) {
    fprintf(fpc,"%d ", captured[c]);
  }
