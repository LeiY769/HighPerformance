#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/**
 * Linearly interpolates between two values.
 * @param start = start value
 * @param end = end value
 * @param interpolation_factor = interpolation factor (0 <= t <= 1)
 */

double lerp(double start, double end, double interpolation_factor) 
{
    return start + interpolation_factor * (end - start);
}

/**
 * Fade function as defined by Ken Perlin. This eases coordinate values
 */
double fade(double t)
{
    return t * t * t * (t * (t * 6 - 15) + 10);
}

/**
 * Generates gradient vectors for each pixel in the heightmap.
 * @param width = width of the heightmap
 * @param height = height of the heightmap
 * @param grad_x = output array for x gradients
 * @param grad_y = output array for y gradients
 */
void generate_gradients(int width, int height, double* grad_x, double* grad_y) 
{
    for (int i = 0; i < height; i++) 
    {
        for (int j = 0; j < width; j++) 
        {
            double angle = (double)rand() / RAND_MAX * 2 * M_PI;
            grad_x[i * width + j] = cos(angle);
            grad_y[i * width + j] = sin(angle);
        }
    }
}

/**
 * Computes Perlin noise at a given 2D coordinate.
 * @param x = x coordinate
 * @param y = y coordinate
 * @param grid_width = width of the noise grid
 * @param grid_height = height of the noise grid
 * @param gradient_x = x gradients
 * @param gradient_y = y gradients
 */

double perlin_noise(double x, double y, int grid_width, int grid_height, double* gradient_x, double* gradient_y) 
{
    int cell_x = (int)floor(x) & (grid_width - 1);
    int cell_y = (int)floor(y) & (grid_height - 1);
    
    double frac_x = x - floor(x);
    double frac_y = y - floor(y);
    

    double fade_x = fade(frac_x);
    double fade_y = fade(frac_y);
    
    int bottom_left = (cell_y * grid_width + cell_x);
    int bottom_right = (cell_y * grid_width + ((cell_x + 1) & (grid_width - 1)));
    int top_left = (((cell_y + 1) & (grid_height - 1)) * grid_width + cell_x);
    int top_right = (((cell_y + 1) & (grid_height - 1)) * grid_width + ((cell_x + 1) & (grid_width - 1)));
    
    double dot_bottom_left = gradient_x[bottom_left] * frac_x + gradient_y[bottom_left] * frac_y;
    double dot_bottom_right = gradient_x[bottom_right] * (frac_x - 1) + gradient_y[bottom_right] * frac_y;
    double dot_top_left = gradient_x[top_left] * frac_x + gradient_y[top_left] * (frac_y - 1);
    double dot_top_right = gradient_x[top_right] * (frac_x - 1) + gradient_y[top_right] * (frac_y - 1);
    
    double interpolated_bottom = lerp(dot_bottom_left, dot_bottom_right, fade_x);
    double interpolated_top = lerp(dot_top_left, dot_top_right, fade_x);
    
    return lerp(interpolated_bottom, interpolated_top, fade_y);
}
/**
 * Generates a heightmap using Perlin noise.
 * @param map = output heightmap
 * @param width = width of the heightmap
 * @param height = height of the heightmap
 * @param octaves = number of octaves for noise generation
 * @param persistence = persistence value for noise generation
 * @param lacunarity = lacunarity value for noise generation
 * @param normalize = normalization factor for the heightmap
 */

void generate_map(double* map, int width, int height, double octaves, double persistence, double lacunarity,int normalize) 
{
    double* grad_x = malloc(width * height * sizeof(double));
    double* grad_y = malloc(width * height * sizeof(double));
    generate_gradients(width, height, grad_x, grad_y);

    double max_val = -1.0;
    double min_val = 1.0;

    for (int i = 0; i < height; i++) 
    {
        for (int j = 0; j < width; j++) 
        {
            float amplitude = 1.0;
            float frequency = 1.0;
            float value = 0.0;

            for (int k = 0; k < octaves; k++) 
            {
                value += perlin_noise(j * frequency / 64.0, i * frequency / 64.0, width, height,grad_x, grad_y) * amplitude;
                amplitude *= persistence;
                frequency *= lacunarity;
            }

            map[i * width + j] = value;
            if (value > max_val) max_val = value;
            if (value < min_val) min_val = value;
        }
    }

    // Normalization between 0 and normalize
    float range = max_val - min_val;
    for (int i = 0; i < width * height; i++) {
        map[i] = ((map[i] - min_val) / range) * normalize;
    }

    free(grad_x);
    free(grad_y);
}


int main(int argc, char *argv[])
{

    char *file_name = "map.dat";
    char *input = argv[1];
    int height;
    int width;
    int seed;
    double dx;
    double dy;
    double octaves;
    double persistence;
    double lacunarity;
    int normalize;

    FILE *fp = fopen(input, "r");
    if (!fp) 
    {
        fprintf(stderr, "Error opening file %s\n", input);
        return 1;
    }
    fscanf(fp, "%d", &height);
    fscanf(fp, "%d", &width);
    fscanf(fp, "%d", &seed);
    fscanf(fp, "%lf", &dx);
    fscanf(fp, "%lf", &dy);
    fscanf(fp, "%lf", &octaves);
    fscanf(fp, "%lf", &persistence);
    fscanf(fp, "%lf", &lacunarity);
    fscanf(fp, "%d", &normalize);
    fclose(fp);

    printf("Height: %d\n", height);
    printf("Width: %d\n", width);
    printf("Seed: %d\n", seed);
    printf("dx: %lf\n", dx);
    printf("dy: %lf\n", dy);
    printf("Random Octaves: %lf\n", octaves);
    printf("Random Persistence: %lf\n", persistence);
    printf("Random Lacunarity: %lf\n", lacunarity);
    printf("Normalize: %d\n", normalize);

    //Init the random seeder generator
    srand(seed);

    double *map = malloc(width * height * sizeof(double));
    if (!map) 
    {
        printf("Error to allocate memory\n");
        return 1;
    }

    generate_map(map, width, height, octaves, persistence, lacunarity, normalize);

    FILE* file = fopen(file_name, "wb");
    if (!file) 
    {
        printf("Error: Unable to create file %s\n", file_name);
        free(map);
        return 1;
    }
    fwrite(&width, sizeof(int), 1, file);
    fwrite(&height, sizeof(int), 1, file);
    fwrite(&dx, sizeof(double), 1, file);
    fwrite(&dy, sizeof(double), 1, file);
    fwrite(map, sizeof(double), width * height, file);

    free(map);
    fclose(file);

    return 0;
}



