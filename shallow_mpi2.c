#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <mpi.h>

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
  int offsetx, offsety; // offsets for subdomain in the global grid
  double dx, dy;
  double *values;
};

#define GET(data, i, j) ((data)->values[(data)->nx * (j) + (i)])
#define SET(data, i, j, val) ((data)->values[(data)->nx * (j) + (i)] = (val))

// Rows and columns in the subdomain
#define TOP_ROW(data) (&(data)->values[0])
#define INTERNAL_FIRST_ROW(data) (&(data)->values[(data)->nx])
#define INTERNAL_LAST_ROW(data) (&(data)->values[((data)->ny - 2) * (data)->nx])
#define BOTTOM_ROW(data) (&(data)->values[((data)->ny - 1) * (data)->nx])

#define LEFT_COL(data, j) ((data)->values[(j) * (data)->nx])
#define INTERNAL_FIRST_COL(data, j) ((data)->values[(j) * (data)->nx + 1])
#define INTERNAL_LAST_COL(data, j) ((data)->values[(j) * (data)->nx + (data)->nx - 2])
#define RIGHT_COL(data, j) ((data)->values[(j) * (data)->nx + (data)->nx - 1])


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
    const char *filename, int step, int rank)
{
char out[512];
if(step < 0)
sprintf(out, "%s_rank%d.vti", filename, rank);
else
sprintf(out, "%s_rank%d_%d.vti", filename, rank, step);

FILE *fp = fopen(out, "wb");
if(!fp) {
printf("Error: Could not open output VTK file '%s'\n", out);
return 1;
}

uint64_t num_points = data->nx * data->ny;
uint64_t num_bytes = num_points * sizeof(double);

fprintf(fp, "<?xml version=\"1.0\"?>\n"
"<VTKFile"
" type=\"ImageData\""
" version=\"1.0\""
" byte_order=\"LittleEndian\""
" header_type=\"UInt64\""
">\n"
"  <ImageData"
" WholeExtent=\"0 %d 0 %d 0 0\""
" Spacing=\"%lf %lf 0.0\""
" Origin=\"%lf %lf 0\""
">\n"
"    <Piece Extent=\"0 %d 0 %d 0 0\">\n"
"      <PointData Scalars=\"scalar_data\">\n"
"        <DataArray"
" type=\"Float64\""
" Name=\"%s\""
" format=\"appended\""
" offset=\"0\""
">\n"
"        </DataArray>\n"
"      </PointData>\n"
"    </Piece>\n"
"  </ImageData>\n"
"  <AppendedData encoding=\"raw\">\n_",
data->nx - 1, data->ny - 1,
data->dx, data->dy,
data->offsetx * data->dx, data->offsety * data->dy,
data->nx - 1, data->ny - 1,
name);

fwrite(&num_bytes, sizeof(uint64_t), 1, fp);
fwrite(data->values, sizeof(double), num_points, fp);

fprintf(fp, "  </AppendedData>\n"
"</VTKFile>\n");

fclose(fp);

return 0;
}

int write_manifest_vtk(const char *filename, double dt, int nt,
        int sampling_rate, int numranks)
{
char out[512];
sprintf(out, "%s.pvd", filename);

FILE *fp = fopen(out, "wb");

if(!fp) {
printf("Error: Could not open output VTK manifest file '%s'\n", out);
return 1;
}

fprintf(fp, "<VTKFile"
" type=\"Collection\""
" version=\"0.1\""
" byte_order=\"LittleEndian\">\n"
"  <Collection>\n");

for(int n = 0; n < nt; n++) {
if(sampling_rate && !(n % sampling_rate)) {
double t = n * dt;
for (int rank = 0; rank < numranks; rank++) {
fprintf(fp, "    <DataSet timestep=\"%g\" part=\"%d\" file='%s_rank%d_%d.vti'/>\n",
 t, rank, filename, rank, n);
}
}
}

fprintf(fp, "  </Collection>\n"
"</VTKFile>\n");

fclose(fp);

return 0;
}

int init_data(struct data *data, int nx, int ny, int offsetx, int offsety, double dx, double dy,
              double val)
{
  data->nx = nx;
  data->ny = ny;
  data->offsetx = offsetx;
  data->offsety = offsety;
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

void free_data(struct data *data)
{
    if(data->values) {
        free(data->values);
        data->values = NULL; 
    }
    return;
}

int check_input(const struct data *data, const struct parameters *param){
    if(data->nx <= 0 || data->ny <= 0) 
    {
        printf("Error: Invalid grid size (%d, %d)\n", data->nx, data->ny);
        return 1;
    }
    if(data->dx <= 0 || data->dy <= 0) 
    {
        printf("Error: Invalid grid spacing (%g, %g)\n", data->dx, data->dy);
        return 1;
    }
    if(!data->values) 
    {
        printf("Error: Data values not initialized\n");
        return 1;
    }
    if(param->dt <= 0 || param->max_t <= 0) 
    {
        printf("Error: Invalid time parameters (dt: %g, max_t: %g)\n", param->dt, param->max_t);
        return 1;
    }
    if(param->g <= 0 || param->gamma < 0) 
    {
        printf("Error: Invalid physical parameters (g: %g, gamma: %g)\n", param->g, param->gamma);
        return 1;
    }
    if(param->sampling_rate <= 0) 
    {
        printf("Error: Invalid sampling rate %d\n", param->sampling_rate);
        return 1;
    }
    return 0;
}
void free_all(struct data *data,struct data *data2, struct data *data3, struct data *data4,struct data *data5)
{
    free_data(data);
    free_data(data2);
    free_data(data3);
    free_data(data4);
    free_data(data5);
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
  double coord_ij = GET(data, i, j);
  double coord_i1j = GET(data, i + 1, j);
  double coord_ij1 = GET(data, i, j + 1);
  double coord_i1j1 = GET(data, i + 1, j + 1);

  //Interpolation
  double val = (1 - tx) * (1 - ty) * coord_ij + tx * (1 - ty) * coord_i1j + (1 - tx) * ty * coord_ij1 + tx * ty * coord_i1j1;
  return val;
}
/* MPI Implementation */

void neighbor_ranks(int rank,int px, int py,int *up_rank, int *down_rank,int *left_rank, int *right_rank)
{
    int row = rank / px;
    int col = rank % px;
    if (row > 0) 
    {
        *up_rank = rank - px; // rank above
    } 
    else 
    {
        *up_rank = MPI_PROC_NULL; // no rank above
    }
    if(row < py - 1) 
    {
        *down_rank = rank + px; // rank below
    } 
    else 
    {
        *down_rank = MPI_PROC_NULL; // no rank below
    }
    if(col > 0) 
    {
        *left_rank = rank - 1; // rank to the left
    } 
    else 
    {
        *left_rank = MPI_PROC_NULL; // no rank to the left
    }
    if(col < px - 1) 
    {
        *right_rank = rank + 1; // rank to the right
    } 
    else 
    {
        *right_rank = MPI_PROC_NULL; // no rank to the right
    }
}


void update(struct data *data, int left_rank, int right_rank,int top_rank, int down_rank)
{
    int ny = data->ny;

    double *send_buffer = (double *)calloc(ny, sizeof(double));
    double *recv_buffer = (double *)calloc(ny, sizeof(double));

    if (left_rank != MPI_PROC_NULL) 
    {
        for (int i = 0; i < ny; i++) 
        {
            send_buffer[i] = INTERNAL_FIRST_COL(data, i);
        }
        MPI_Sendrecv(send_buffer, ny, MPI_DOUBLE, left_rank, 0, recv_buffer, ny,MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < ny; i++) 
        {
            LEFT_COL(data, i) = recv_buffer[i];
        }
    }
    if (right_rank != MPI_PROC_NULL) 
    {
        for (int i = 0; i < ny; i++) 
        {
            send_buffer[i] = INTERNAL_LAST_COL(data, i);
        }
        MPI_Sendrecv(send_buffer, ny, MPI_DOUBLE, right_rank, 0, recv_buffer, ny,MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < ny; i++) 
        {
            RIGHT_COL(data, i) = recv_buffer[i];
        }
    }

    int nx = data->nx;
    if(top_rank != MPI_PROC_NULL) 
    {
        MPI_Sendrecv(INTERNAL_FIRST_ROW(data), nx, MPI_DOUBLE, top_rank, 0,TOP_ROW(data), nx, MPI_DOUBLE, top_rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    if(down_rank != MPI_PROC_NULL) 
    {
        MPI_Sendrecv(INTERNAL_LAST_ROW(data), nx, MPI_DOUBLE, down_rank, 0,BOTTOM_ROW(data), nx, MPI_DOUBLE, down_rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    free(send_buffer);
    free(recv_buffer);
}


void number_partitions(int world_size, int nx, int ny, int *Px, int *Py) 
{
    double best_ratio = DBL_MAX;
    *Px = *Py = 0;

    for (int px = 1; px <= world_size; px++) {
        if (world_size % px != 0) continue;
        int py = world_size / px;
        
        // Calculate aspect ratio of subdomains
        double rx = (double)nx / px;
        double ry = (double)ny / py;
        double ratio = fmax(rx, ry) / fmin(rx, ry);

        if (ratio < best_ratio) {
            best_ratio = ratio;
            *Px = px;
            *Py = py;
        }
    }
}


int main(int argc, char **argv)
{
    int world_size;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(argc != 2) 
    {
        if(rank == 0) printf("Usage: %s parameter_file\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    struct parameters param;
    if(read_parameters(&param, argv[1]))
    {
        MPI_Finalize();
        return 1;
    }
    if (rank == 0) 
    {
        print_parameters(&param);
        printf("World size: %d\n", world_size);
    }

    struct data h;
    if(read_data(&h, param.input_h_filename))
    { 
        MPI_Finalize();
        return 1;
    }

    if(check_input(&h, &param))
    {
        free_data(&h);
        printf("Input data check failed.\n");
        MPI_Finalize();
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

    if (rank == 0) 
    {
        printf(" - grid size: %g m x %g m (%d x %d = %d grid points)\n",hx, hy, nx, ny, nx * ny);
        printf(" - number of time steps: %d\n", nt);
    }
    
    // Calculate optimal domain decomposition
    int Px = 0;
    int Py = 0;
    number_partitions(world_size, nx, ny, &Px, &Py);
    
    if (Py == 0 || Px == 0)
    {
        if (rank == 0) 
        {
            printf("Error: Cannot determine partitions for %d ranks\n", world_size);
        }
        MPI_Finalize();
        return 1;
    }
    else
    {
        if (rank == 0) 
        {
            printf(" - number of ranks: %d (%d x %d)\n", world_size, Px, Py);
        }
    }
    
    // Calculate local domain sizes with remainder distribution
    int base_nx = nx / Px;
    int rem_nx = nx % Px;
    int base_ny = ny / Py;
    int rem_ny = ny % Py;

    int col = rank % Px;
    int row = rank / Px;

    int local_nx = base_nx + (col < rem_nx ? 1 : 0);
    int local_ny = base_ny + (row < rem_ny ? 1 : 0);

    int offset_x = col < rem_nx ? 
                  col * (base_nx + 1) : 
                  rem_nx * (base_nx + 1) + (col - rem_nx) * base_nx;
                  
    int offset_y = row < rem_ny ? 
                  row * (base_ny + 1) : 
                  rem_ny * (base_ny + 1) + (row - rem_ny) * base_ny;

    int top_rank = MPI_PROC_NULL;
    int bottom_rank = MPI_PROC_NULL;
    int left_rank = MPI_PROC_NULL;
    int right_rank = MPI_PROC_NULL;
    neighbor_ranks(rank, Px, Py, &top_rank, &bottom_rank, &left_rank, &right_rank);

    if (rank == 0) {
        printf("Domain decomposition:\n");
        printf(" - Global size: %d x %d\n", nx, ny);
        printf(" - Base cell size: %d x %d\n", base_nx, base_ny);
        printf(" - Remainder: %d in x, %d in y\n", rem_nx, rem_ny);
    }
    printf("Rank %d: local size = %d x %d, offset = (%d, %d)\n", 
           rank, local_nx, local_ny, offset_x, offset_y);

    struct data eta, u, v;
    init_data(&eta, local_nx, local_ny, offset_x, offset_y, param.dx, param.dy, 0.);
    init_data(&u, local_nx + 1, local_ny, offset_x, offset_y, param.dx, param.dy, 0.);
    init_data(&v, local_nx, local_ny + 1, offset_x, offset_y, param.dx, param.dy, 0.);
    
    // Interpolation initialization
    struct data h_u;
    struct data h_v;
    init_data(&h_u, local_nx + 1, local_ny, offset_x, offset_y, param.dx, param.dy, 0.);
    init_data(&h_v, local_nx, local_ny + 1, offset_x, offset_y, param.dx, param.dy, 0.);

    // Interpolate bathymetry to u and v grids
    for(int j = 0; j < local_ny; j++) 
    {
        for(int i = 0; i < local_nx + 1; i++) 
        {
            double h_xu = (offset_x + i) * param.dx;
            double h_yu = (offset_y + j + 0.5) * param.dy;
            double val = interpolate_data(&h, h_xu, h_yu, invdx, invdy, param.dx, param.dy);
            SET(&h_u, i, j, val);
        }
    }
    for(int j = 0; j < local_ny + 1; j++) 
    {
        for(int i = 0; i < local_nx; i++) 
        {
            double h_xv = (offset_x + i + 0.5) * param.dx;
            double h_yv = (offset_y + j) * param.dy;
            double val = interpolate_data(&h, h_xv, h_yv, invdx, invdy, param.dx, param.dy);
            SET(&h_v, i, j, val);
        }
    }
    update(&h_u, left_rank, right_rank, top_rank, bottom_rank);
    update(&h_v, left_rank, right_rank, top_rank, bottom_rank);

    double start = GET_TIME();

    // time loop
    for(int n = 0; n < nt; n++) 
    {
        if(n && (n % (nt / 10)) == 0) 
        {
        double time_sofar = GET_TIME() - start;
        double eta = (nt - n) * time_sofar / n;
        if(rank == 0)
        {
            printf("Computing step %d/%d (ETA: %g seconds)     \r", n, nt, eta);
            fflush(stdout);
        }
        }

        // output solution
        if(param.sampling_rate && !(n % param.sampling_rate)) 
        {
            write_data_vtk(&eta, "water elevation", param.output_eta_filename, n, rank);
        }

        // impose boundary conditions
        double t = n * param.dt;
        double sin_t = sin(alpha * t);
        if(param.source_type == 1) 
        {
            // sinusoidal velocity on top boundary
            if (offset_y + local_ny == ny) 
            { // top global boundary
                for(int i = 0; i < local_nx; i++) 
                {
                    SET(&v, i, local_ny, A * sin_t);
                }
            }
        }
        else if(param.source_type == 2) 
        {
            // sinusoidal elevation in the middle of the domain
            int mid_x = nx / 2;
            int mid_y = ny / 2;
            if (offset_x <= mid_x && mid_x < offset_x + local_nx &&
                offset_y <= mid_y && mid_y < offset_y + local_ny) 
            {
                int i = mid_x - offset_x;
                int j = mid_y - offset_y;
                SET(&eta, i, j, A * sin_t);
            }
        }
        else 
        {
            if(rank == 0) printf("Error: Unknown source type %d\n", param.source_type);
            MPI_Finalize();
            exit(0);
        }

        // update eta
        for(int j = 0; j < local_ny ; j++) 
        {
            for(int i = 0; i < local_nx; i++) 
            {
                double current_eta = GET(&eta, i, j);
                double eta_ij = current_eta
                - param.dt * invdx * (GET(&h_u, i + 1, j) * GET(&u, i+1,j) - GET(&h_u, i, j) * GET(&u, i, j))
                - param.dt * invdy * (GET(&h_v, i, j + 1) * GET(&v,i,j+1) - GET(&h_v, i, j) * GET(&v, i, j));
                SET(&eta, i, j, eta_ij);
            }
        }
        update(&eta, left_rank, right_rank, top_rank, bottom_rank);

        // update u and v
        for(int j = 0; j < local_ny; j++) 
        {
            for(int i = 0; i < local_nx; i++) 
            {
                double eta_ij = GET(&eta, i, j);
                double eta_imj = (i == 0) ? GET(&eta, 0, j) : GET(&eta, i-1, j);
                double eta_ijm = (j == 0) ? GET(&eta, i, 0) : GET(&eta, i, j-1);
                
                double u_ij = (1. - c2) * GET(&u, i, j)
                - (c1_dx) * (eta_ij - eta_imj);
                double v_ij = (1. - c2) * GET(&v, i, j)
                - (c1_dy) * (eta_ij - eta_ijm);
                
                SET(&u, i, j, u_ij);
                SET(&v, i, j, v_ij);
            }
        }
        update(&u, left_rank, right_rank, top_rank, bottom_rank);
        update(&v, left_rank, right_rank, top_rank, bottom_rank);
    }

    if (rank == 0) {
        write_manifest_vtk(param.output_eta_filename, param.dt, nt, param.sampling_rate, world_size);
    }
    
    if (rank == 0)
    {
        double time = GET_TIME() - start;
        printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
               1e-6 * (double)nx * (double)ny * (double)nt / time);
    }
    
    free_all(&eta, &u, &v, &h_u, &h_v);
    free_data(&h);

    MPI_Finalize();
    return 0;
}