#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>

#if defined(_OPENMP)
#include <omp.h>
#define GET_TIME() (omp_get_wtime()) // wall time
#else
#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time
#endif

struct parameters {
  double dx, dy, dt, max_t;
  double g, gamma;
  int source_type;
  int sampling_rate;
  char input_h_filename[256];
  char output_eta_filename[256];
  char output_u_filename[256];
  char output_v_filename[256];
};


struct data {
  int nx, ny;
  double dx, dy;
  double *values;
};

static inline double get_value(const struct data *data, int i ,int j){
    return data->values[(data->nx * j) + i];}

static inline void set_value(struct data *data, int i, int j , double val){
    data->values[(data->nx * j) + i] = val;}

int read_parameters(struct parameters *param, const char *filename)
{
  FILE *fp = fopen(filename, "r");
  if(!fp) {
    printf("Error: Could not open parameter file '%s'\n", filename);
    return 1;
  }
  int ok = 1;
  if(ok) ok = (fscanf(fp, "%lf", &param->dx) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->dy) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->dt) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->max_t) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->g) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->gamma) == 1);
  if(ok) ok = (fscanf(fp, "%d", &param->source_type) == 1);
  if(ok) ok = (fscanf(fp, "%d", &param->sampling_rate) == 1);
  if(ok) ok = (fscanf(fp, "%256s", param->input_h_filename) == 1);
  if(ok) ok = (fscanf(fp, "%256s", param->output_eta_filename) == 1);
  if(ok) ok = (fscanf(fp, "%256s", param->output_u_filename) == 1);
  if(ok) ok = (fscanf(fp, "%256s", param->output_v_filename) == 1);
  fclose(fp);
  if(!ok) {
    printf("Error: Could not read one or more parameters in '%s'\n", filename);
    return 1;
  }
  return 0;
}

void print_parameters(const struct parameters *param)
{
  printf("Parameters:\n");
  printf(" - grid spacing (dx, dy): %g m, %g m\n", param->dx, param->dy);
  printf(" - time step (dt): %g s\n", param->dt);
  printf(" - maximum time (max_t): %g s\n", param->max_t);
  printf(" - gravitational acceleration (g): %g m/s^2\n", param->g);
  printf(" - dissipation coefficient (gamma): %g 1/s\n", param->gamma);
  printf(" - source type: %d\n", param->source_type);
  printf(" - sampling rate: %d\n", param->sampling_rate);
  printf(" - input bathymetry (h) file: '%s'\n", param->input_h_filename);
  printf(" - output elevation (eta) file: '%s'\n", param->output_eta_filename);
  printf(" - output velocity (u, v) files: '%s', '%s'\n",
         param->output_u_filename, param->output_v_filename);
}

int read_data(struct data *data, const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  if(!fp) {
    printf("Error: Could not open input data file '%s'\n", filename);
    return 1;
  }
  int ok = 1;
  if(ok) ok = (fread(&data->nx, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fread(&data->ny, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fread(&data->dx, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fread(&data->dy, sizeof(double), 1, fp) == 1);
  if(ok) {
    int N = data->nx * data->ny;
    if(N <= 0) {
      printf("Error: Invalid number of data points %d\n", N);
      ok = 0;
    }
    else {
      data->values = (double*)malloc(N * sizeof(double));
      if(!data->values) {
        printf("Error: Could not allocate data (%d doubles)\n", N);
        ok = 0;
      }
      else {
        ok = (fread(data->values, sizeof(double), N, fp) == N);
      }
    }
  }
  fclose(fp);
  if(!ok) {
    printf("Error reading input data file '%s'\n", filename);
    return 1;
  }
  return 0;
}

int write_data(const struct data *data, const char *filename, int step)
{
  char out[512];
  if(step < 0)
    sprintf(out, "%s.dat", filename);
  else
    sprintf(out, "%s_%d.dat", filename, step);
  FILE *fp = fopen(out, "wb");
  if(!fp) {
    printf("Error: Could not open output data file '%s'\n", out);
    return 1;
  }
  int ok = 1;
  if(ok) ok = (fwrite(&data->nx, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fwrite(&data->ny, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fwrite(&data->dx, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fwrite(&data->dy, sizeof(double), 1, fp) == 1);
  int N = data->nx * data->ny;
  if(ok) ok = (fwrite(data->values, sizeof(double), N, fp) == N);
  fclose(fp);
  if(!ok) {
    printf("Error writing data file '%s'\n", out);
    return 1;
  }
  return 0;
}

int write_data_vtk(const struct data *data, const char *name,
                   const char *filename, int step)
{
  char out[512];
  if(step < 0)
    sprintf(out, "%s.vti", filename);
  else
    sprintf(out, "%s_%d.vti", filename, step);

  FILE *fp = fopen(out, "wb");
  if(!fp) {
    printf("Error: Could not open output VTK file '%s'\n", out);
    return 1;
  }

  uint64_t num_points = data->nx * data->ny;
  uint64_t num_bytes = num_points * sizeof(double);

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" "
          "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" "
          "Spacing=\"%lf %lf 0.0\">\n",
          data->nx - 1, data->ny - 1, data->dx, data->dy);
  fprintf(fp, "    <Piece Extent=\"0 %d 0 %d 0 0\">\n",
          data->nx - 1, data->ny - 1);

  fprintf(fp, "      <PointData Scalars=\"scalar_data\">\n");
  fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
          "format=\"appended\" offset=\"0\">\n", name);
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </PointData>\n");

  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </ImageData>\n");

  fprintf(fp, "  <AppendedData encoding=\"raw\">\n_");

  fwrite(&num_bytes, sizeof(uint64_t), 1, fp);
  fwrite(data->values, sizeof(double), num_points, fp);

  fprintf(fp, "  </AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");

  fclose(fp);
  return 0;
}

int write_manifest_vtk(const char *filename, double dt, int nt,
                       int sampling_rate)
{
  char out[512];
  sprintf(out, "%s.pvd", filename);

  FILE *fp = fopen(out, "wb");
  if(!fp) {
    printf("Error: Could not open output VTK manifest file '%s'\n", out);
    return 1;
  }

  fprintf(fp, "<VTKFile type=\"Collection\" version=\"0.1\" "
          "byte_order=\"LittleEndian\">\n");
  fprintf(fp, "  <Collection>\n");
  for(int n = 0; n < nt; n++) {
    if(sampling_rate && !(n % sampling_rate)) {
      double t = n * dt;
      fprintf(fp, "    <DataSet timestep=\"%g\" file='%s_%d.vti'/>\n", t,
              filename, n);
    }
  }
  fprintf(fp, "  </Collection>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
  return 0;
}

int init_data(struct data *data, int nx, int ny, double dx, double dy,double val)
{
  data->nx = nx;
  data->ny = ny;
  data->dx = dx;
  data->dy = dy;
  data->values = (double*)malloc(nx * ny * sizeof(double));
  if(!data->values){
    printf("Error: Could not allocate data\n");
    return 1;
  }
  for(int i = 0; i < nx * ny; i++) data->values[i] = val;
  return 0;
}

int check_input(const struct data *data, const struct parameters *param){
  if(data->nx <= 0 || data->ny <= 0) 
  {
      printf("Error: Invalid grid size (%d, %d)\n", data->nx, data->ny);
      return 1;
  }
  if(data->dx <= 0 || data->dy <= 0) 
  {
      printf("Error: Invalid grid spacing (%f, %f)\n", data->dx, data->dy);
      return 1;
  }
  if(param->dt <= 0 || param->max_t <= 0) 
  {
      printf("Error: Invalid time parameters (dt: %f, max_t: %f)\n", param->dt, param->max_t);
      return 1;
  }
  if(param->g <= 0 || param->gamma < 0) 
  {
      printf("Error: Invalid physical parameters (g: %3.f, gamma: %3.f)\n", param->g, param->gamma);
      return 1;
  }
  if(param->sampling_rate < 0) 
  {
      printf("Error: Invalid sampling rate %d\n", param->sampling_rate);
      return 1;
  }
  return 0;
}

void free_data(struct data *data)
{
    if(data->values) {
        free(data->values);
        data->values = NULL; 
    }
    return;
}

void free_all(struct data *data,struct data *data2, struct data *data3, struct data *data4,struct data *data5,struct data *data6)
{
    free_data(data);
    free_data(data2);
    free_data(data3);
    free_data(data4);
    free_data(data5);
    free_data(data6);
    return;
}


double interpolate_data(const struct data *data, double x, double y,double invdx,double invdy, double dx, double dy)
{
  int i = (int)(floor(x * invdx));
  int j = (int)(floor(y * invdy));

  double tx = (x - i * dx) * invdx;
  double ty = (y - j * dy) * invdy;
  if(i < 0) i = 0;
  else if(i > data->nx - 2) i = data->nx - 2; // Be sure to not go out of bounds
  if(j < 0) j = 0;
  else if(j > data->ny - 2) j = data->ny - 2; // Be sure to not go out of bounds

  // Get the corners
  double coord_ij = get_value(data, i, j);
  double coord_i1j = get_value(data, i + 1, j);
  double coord_ij1 = get_value(data, i, j + 1);
  double coord_i1j1 = get_value(data, i + 1, j + 1);
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

  struct parameters param;
  if(read_parameters(&param, argv[1])) return 1;
  print_parameters(&param);

  struct data h;
  if(read_data(&h, param.input_h_filename)) return 1;

  if(check_input(&h, &param))
  {
    free_data(&h);
    printf("Input data check failed.\n");
    return 1;
  }
  // Variable declaration
  double A = 5;
  double f = 1. / 20.;
  double c1 = param.dt * param.g;
  double c2 = param.dt * param.gamma;
  double alpha = 2* M_PI *f;

  // Constant factor to get the multiplication instead of division
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

  struct data eta, u, v;
    if (init_data(&eta, nx, ny, param.dx, param.dy, 0.))
    {
        free_data(&h);
        return 1;
    }
    if (init_data(&u, nx + 1, ny, param.dx, param.dy, 0.)) 
    {
        free_all(&eta, &h, NULL, NULL, NULL,NULL);
        return 1;
    }
    if (init_data(&v, nx, ny + 1, param.dx, param.dy, 0.)) 
    {
        free_all(&eta, &u, NULL, NULL, NULL,&h);
        return 1;
    }
    // Interpolation initialization
    struct data h_u;
    if (init_data(&h_u, nx + 1, ny, param.dx, param.dy, 0.)) 
    {
        free_all(&eta, &u, &v, NULL, NULL,&h);
        return 1;
    }
  struct data h_v;
    if (init_data(&h_v, nx, ny + 1, param.dx, param.dy, 0.)) 
    {
        free_all(&eta, &u, &v, &h_u, NULL,&h);
        return 1;
    }
  for(int j = 0; j < ny; j++) 
  {
    for(int i = 0; i < nx+1; i++) 
    {
      double h_xu = i * param.dx;
      double h_yu = (j+0.5) * param.dy;
      double val = interpolate_data(&h, h_xu, h_yu,invdx,invdy,param.dx,param.dy);
      set_value(&h_u, i, j, val);
    }
  }
  for(int j = 0; j < ny+1; j++) 
  {
    for(int i = 0; i < nx; i++) 
    {
      double h_xv = (i+0.5) * param.dx;
      double h_yv = j * param.dy;
      double val = interpolate_data(&h, h_xv, h_yv,invdx,invdy,param.dx,param.dy);
      set_value(&h_v, i, j, val);
    }
  }


  double start = GET_TIME();


  // time loop
  for(int n = 0; n < nt; n++) 
  {

    if(n && (n % (nt / 10)) == 0) 
    {
      double time_sofar = GET_TIME() - start;
      double eta = (nt - n) * time_sofar / n;
      printf("Computing step %d/%d (ETA: %g seconds)     \r", n, nt, eta);
      fflush(stdout);
    }

    // output solution
    if(param.sampling_rate && !(n % param.sampling_rate)) 
    {
        write_data_vtk(&eta, "water elevation", param.output_eta_filename, n);
        //write_data_vtk(&u, "x velocity", param.output_u_filename, n);
        //write_data_vtk(&v, "y velocity", param.output_v_filename, n);
  
    }

    // impose boundary conditions
    double t = n * param.dt;
    double sin_t = sin(alpha * t);
    if(param.source_type == 1) 
    {
      // sinusoidal velocity on top boundary
      for(int j = 0; j < nx; j++) 
      {
        set_value(&u, 0, j, 0.);
        set_value(&u, nx, j, 0.);

      }
      for(int i = 0; i < ny; i++) 
      {
        set_value(&v, i, 0, 0.);
        set_value(&v, i, ny, A * sin_t);
      }
    }
    else if(param.source_type == 2) 
    {
      // sinusoidal elevation in the middle of the domain
      set_value(&eta, nx / 2, ny / 2, A * sin_t);
    }
    else if(param.source_type == 3)
    {
      for(int j = 0; j < nx; j++) 
      {
        set_value(&u, 0, j, 0.);
        set_value(&u, nx, j, 0.);

      }
      for(int i = 0; i < ny; i++) 
      {
        set_value(&v, i, 0, -A * sin_t);
        set_value(&v, i, ny, A * sin_t);
      }

      set_value(&eta, 3* nx / 4, ny / 2, A * sin_t);
      set_value(&eta, nx / 4, ny / 2, A * sin_t);

    }
    else 
    {
      printf("Error: Unknown source type %d\n", param.source_type);
      free_all(&eta, &u, &v, &h_u, &h_v, &h);
      exit(0);
    }

    // update eta
    for(int j = 0; j < ny ; j++) 
    {
      for(int i = 0; i < nx; i++) 
      {
        double current_eta = get_value(&eta, i, j);
        double eta_ij = current_eta
          - param.dt * invdx * (get_value(&h_u, i + 1, j) * get_value(&u, i+1,j) - get_value(&h_u, i, j) * get_value(&u, i, j))
          - param.dt * invdy * (get_value(&h_v, i, j + 1) * get_value(&v,i,j+1) - get_value(&h_v, i, j) * get_value(&v, i, j));
        set_value(&eta, i, j, eta_ij);
      }
    }

    // update u and v
    for(int j = 0; j < ny; j++) 
    {
      for(int i = 0; i < nx; i++) 
      {
        double eta_ij = get_value(&eta, i, j);
        double eta_imj = get_value(&eta, (i == 0) ? 0 : i - 1, j);
        double eta_ijm = get_value(&eta, i, (j == 0) ? 0 : j - 1);
        double u_ij = (1. - c2) * get_value(&u, i, j)
          - (c1_dx) * (eta_ij - eta_imj);
        double v_ij = (1. - c2) * get_value(&v, i, j)
          - (c1_dy) * (eta_ij - eta_ijm);
        set_value(&u, i, j, u_ij);
        set_value(&v, i, j, v_ij);
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
  printf("\nDone: %g seconds (%g MUpdates/s)\n", time,1e-6 * (double)eta.nx * (double)eta.ny * (double)nt / time);
    free_all(&eta, &u, &v, &h_u, &h_v,&h);
         


  return 0;
}