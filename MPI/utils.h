#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>

/**
 * Structure to hold parameters.
 * @param dx,dy: grid spacing in x and y directions
 * @param dt: time step
 * @param max_t: maximum time
 * @param g: gravitational acceleration
 * @param gamma: dissipation coefficient
 * @param source_type: type of source
 * @param sampling_rate: rate at which to sample the solution
 * @param input_h_filename: name of the input bathymetry file
 * @param output_eta_filename: name of the output elevation file
 * @param output_u_filename: name of the output x velocity file
 * @param output_v_filename: name of the output y velocity file
 */
typedef struct parameters_t {
  double dx, dy, dt, max_t;
  double g, gamma;
  int source_type;
  int sampling_rate;
  char input_h_filename[256];
  char output_eta_filename[256];
  char output_u_filename[256];
  char output_v_filename[256];
}parameters;
/**
 * Structure to hold data.
 * @param nx,ny: number of grid points in x and y directions
 * @param dx,dy: grid spacing in x and y directions
 * @param values: array of values
 */
typedef struct data_t {
  int nx, ny;
  double dx, dy;
  double *values;
}data;

/**
 * Read parameters from a file.
 * @param param pointer to a parameters structure
 * @param filename name of the file to read
 * @return 0 if successful, 1 otherwise
 */
int read_parameters(parameters *param, const char *filename);
/**
 * Print parameters to the console.
 * @param param pointer to a parameters structure
 */
void print_parameters(const parameters *param);
/**
 * Write the data into .dat file format.
 * @param data pointer to a data structure
 * @param filename name of the file to write
 * @param step step number
 * @return 0 if successful, 1 otherwise
 */
int write_data(const data *data, const char *filename, int step);
/**
 * Read data from a file.
 * @param data pointer to a data structure
 * @param filename name of the file to read
 * @return 0 if successful, 1 otherwise
 */
int read_data(data *data, const char *filename);
/**
 * Write the data into .vtk file the content of at the time step
 * @param data pointer to a data structure
 * @param name name of the file to hold for summary the different files
 * @param filename name of the file to write
 * @param step The current step
 * @return 0 if successful, 1 otherwise
 */
int write_data_vtk(const data *data, const char *name,const char *filename, int step);
/**
 * Write the global manifest file for the vtk files .pwd
 * @param filename name of the file to write
 * @param dt time step
 * @param nt number of time steps
 * @param sampling_rate rate at which to sample the solution
 * @return 0 if successful, 1 otherwise
 */
int write_manifest_vtk(const char *filename,double dt, int nt, int sampling_rate);

/**
 * Initialize the data structure.
 * @param data pointer to a data structure
 * @param nx number of grid points in x direction
 * @param ny number of grid points in y direction
 * @param dx grid spacing in x direction
 * @param dy grid spacing in y direction
 * @param val initial value
 */
int init_data(data *data, int nx, int ny, double dx, double dy,double val);

/**
 * Free the data structure.
 * @param data pointer to a data structure
 */
void free_data(data *data);



#endif