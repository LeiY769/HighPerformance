#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "utils.h"

#if defined(_OPENMP)
#include <omp.h>
#define GET_TIME() (omp_get_wtime()) // wall time
#else
#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time
#endif

#define GET(data, i, j) ((data)->values[(data)->nx * (j) + (i)])
#define SET(data, i, j, val) ((data)->values[(data)->nx * (j) + (i)] = (val))

double interpolate_data(const data *data, double x, double y,double invdx,double invdy)
{
  int i = (int)(x * invdx);
  int j = (int)(y * invdy);

  double tx = (x - i / invdx) * invdx;
  double ty = (y - j / invdy) * invdy;
  if(i < 0) i = 0;
  else if(i > data->nx - 2) i = data->nx - 2; // Be sure to not go out of bounds
  if(j < 0) j = 0;
  else if(j > data->ny - 2) j = data->ny - 2; // Be sure to not go out of bounds

  // Get the corners
  double coord_ij = GET(data, i, j);
  double coord_i1j = GET(data, i + 1, j);
  double coord_ij1 = GET(data, i, j + 1);
  double coord_i1j1 = GET(data, i + 1, j + 1);

  //Interpolation
  double val = (1 - tx) * (1 - ty) * coord_ij + tx * (1 - ty) * coord_i1j + (1 - tx) * ty * coord_ij1 + tx * ty * coord_i1j1;
  return val;
}

int main(int argc, char **argv)
{
  if(argc != 2) {
    printf("Usage: %s parameter_file\n", argv[0]);
    return 1;
  }

  parameters param;
  if(read_parameters(&param, argv[1])) return 1;
  print_parameters(&param);

  data h;
  if(read_data(&h, param.input_h_filename)) return 1;
  // Variable declaration
  double A = 5;
  double f = 1. / 20.;
  double c1 = param.dt * param.g;
  double c2 = param.dt * param.gamma;
  double alpha = 2* M_PI *f;
  double invdx = 1./param.dx;
  double invdy = 1./param.dy;
  double c1_dx = c1 / param.dx;
  double c1_dy = c1 / param.dy;

  // infer size of domain from input elevation data
  double hx = h.nx * h.dx;
  double hy = h.ny * h.dy;
  int nx = floor(hx * invdx);
  int ny = floor(hy * invdy);
  if(nx <= 0) nx = 1;
  if(ny <= 0) ny = 1;
  int nt = floor(param.max_t / param.dt);

  printf(" - grid size: %g m x %g m (%d x %d = %d grid points)\n",
         hx, hy, nx, ny, nx * ny);
  printf(" - number of time steps: %d\n", nt);

  data eta, u, v;
  init_data(&eta, nx, ny, param.dx, param.dy, 0.);
  init_data(&u, nx + 1, ny, param.dx, param.dy, 0.);
  init_data(&v, nx, ny + 1, param.dx, param.dy, 0.);

  // interpolate bathymetry
  data h_interp;
  init_data(&h_interp, nx, ny, param.dx, param.dy, 0.);
  for(int j = 0; j < ny; j++) {
    for(int i = 0; i < nx; i++) {
      double x = i * param.dx;
      double y = j * param.dy;
      double val = interpolate_data(&h, x, y,invdx,invdy);
      SET(&h_interp, i, j, val);
    }
  }

  double start = GET_TIME();


  // time loop
  for(int n = 0; n < nt; n++) {

    if(n && (n % (nt / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta = (nt - n) * time_sofar / n;
      printf("Computing step %d/%d (ETA: %g seconds)     \r", n, nt, eta);
      fflush(stdout);
    }

    // output solution
    if(param.sampling_rate && !(n % param.sampling_rate)) {
      write_data_vtk(&eta, "water elevation", param.output_eta_filename, n);
      //write_data_vtk(&u, "x velocity", param.output_u_filename, n);
      //write_data_vtk(&v, "y velocity", param.output_v_filename, n);
    }

    // impose boundary conditions
    double t = n * param.dt;
    double sin_t = sin(alpha * t);
    if(param.source_type == 1) {
      // sinusoidal velocity on top boundary
      for(int j = 0; j < ny; j++) {
        for(int i = 0; i < nx; i++) {
          SET(&u, 0, j, 0.);
          SET(&u, nx, j, 0.);
          SET(&v, i, 0, 0.);
          SET(&v, i, ny, A * sin_t);
        }
      }
    }
    else if(param.source_type == 2) {
      // sinusoidal elevation in the middle of the domain
      SET(&eta, nx / 2, ny / 2, A * sin_t);
    }
    else {
      // TODO: add other sources
      printf("Error: Unknown source type %d\n", param.source_type);
      exit(0);
    }

    // update eta
    for(int j = 0; j < ny ; j++) {
      for(int i = 0; i < nx; i++) {
        // TODO: this does not evaluate h at the correct locations
        double h_ij = GET(&h_interp, i, j);
        double c1_eta = param.dt * h_ij;
        double eta_ij = GET(&eta, i, j)
          - c1_eta * invdx * (GET(&u, i + 1, j) - GET(&u, i, j))
          - c1_eta * invdy * (GET(&v, i, j + 1) - GET(&v, i, j));
        SET(&eta, i, j, eta_ij);
      }
    }

    // update u and v
    for(int j = 0; j < ny; j++) {
      for(int i = 0; i < nx; i++) {
        double eta_ij = GET(&eta, i, j);
        double eta_imj = GET(&eta, (i == 0) ? 0 : i - 1, j);
        double eta_ijm = GET(&eta, i, (j == 0) ? 0 : j - 1);
        double u_ij = (1. - c2) * GET(&u, i, j)
          - (c1_dx) * (eta_ij - eta_imj);
        double v_ij = (1. - c2) * GET(&v, i, j)
          - (c1_dy) * (eta_ij - eta_ijm);
        SET(&u, i, j, u_ij);
        SET(&v, i, j, v_ij);
      }
    }
  }

  write_manifest_vtk(param.output_eta_filename,
                     param.dt, nt, param.sampling_rate);
  //write_manifest_vtk("x velocity", param.output_u_filename,
  //                   param.dt, nt, param.sampling_rate);
  //write_manifest_vtk("y velocity", param.output_v_filename,
  //                   param.dt, nt, param.sampling_rate);

  double time = GET_TIME() - start;
  printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
         1e-6 * (double)eta.nx * (double)eta.ny * (double)nt / time);

  free_data(&h_interp);
  free_data(&eta);
  free_data(&u);
  free_data(&v);

  return 0;
}