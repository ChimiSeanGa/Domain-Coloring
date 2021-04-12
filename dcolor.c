#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#define MIN_X -5
#define MAX_X 5
#define SAMPLE 1500
#define SEQ_LEN 50
#define CONTOUR_BASE 1
#define MAN_BOUND 2
#define MAN_ITER 500
#define ZETA_ITER 100

// Mandelbrot constants
double man_magnify = 1;
complex man_center = -0.7076715936 + 0.3528062781 * I;

// Weierstrass constant
complex om1 = I + 1;
complex om2 = I;

// Struct for storing RGB values.
typedef struct {
   int r;
   int g;
   int b;
} tRGB;

// Struct for storing HLS values.
typedef struct {
   double h;
   double l;
   double s;
} tHLS;

complex *seq;

// Prints a complex number.
void print_complex(complex z) {
   printf("%.2f + %.2fi\n", creal(z), cimag(z));
}

// Simple factorial function.
int factorial(int n) {
   if (n == 0) {
      return 1;
   }
   return n * factorial(n - 1);
}

// Factor for the Blaschke product.
complex b_factor(complex x, complex a) {
   return (cabs(a) / a) * (a - x) / (1 - conj(a) * x);
}

// Blaschke product function.
complex b(complex x) {
   complex res = 1;
   for (int i = 0; i < SEQ_LEN; i++) {
      res *= b_factor(x, seq[i]);
   }
   return res;
}

// Iterate the Mandelbrot function.
complex m_iter(complex z, complex c) {
   return cpow(z, 2) + c;
}

// Colors the complex plane via the Mandelbrot set.
complex m(complex c) {
   complex c_mod = (c / man_magnify) + man_center;
   complex cur_iter = m_iter(0, c_mod);
   int iter_num = 0;
   for (int i = 0; i < MAN_ITER; i++) {
      if (cabs(cur_iter) > MAN_BOUND) {
         return cpow(M_E, ((complex) iter_num) / MAN_ITER * 4 * M_PI * I);
      }
      cur_iter = m_iter(cur_iter, c_mod);
      iter_num++;
   }
   return 0;
}

// Weierstrass p function.
complex wp(complex z) {
   int up_bound = 50;
   complex a = cpow((M_PI / (2 * om1)), 2);
   complex sum1 = 0;
   complex sum2 = 0;
   complex s1, s2, s3, s4;

   for (int i = 0; i < up_bound; i++) {
      s1 = 1 / cpow(csin((z - 2 * i * om2) / (2 * om1) * M_PI), 2);
      s2 = 1 / cpow(csin((z - 2 * -i * om2) / (2 * om1) * M_PI), 2);

      if (!isnan(creal(s1)) && !isnan(cimag(s1))) {
         sum1 += s1;
      }
      if (!isnan(creal(s2)) && !isnan(cimag(s2))) {
         sum1 += s2;
      }

      s3 = 1 / cpow(csin((i * om2) / om1 * M_PI), 2);
      s4 = 1 / cpow(csin((-i * om2) / om1 * M_PI), 2);

      if (!isnan(creal(s3)) && !isnan(cimag(s3))) {
         sum2 += s3;
      }
      if (!isnan(creal(s4)) && !isnan(cimag(s4))) {
         sum2 += s4;
      }
   }
   return a * (-1 / 3 + sum1 - sum2);
}

// Derivative of the Weierstrass p function.
complex wp_prime(complex z) {
   return 4 * (wp(z) - wp(om1 / 2)) * (wp(z) - wp(om2 / 2))
      * (wp(z) - wp((om1 + om2) / 2));
}

// Adds noise to the spiral function.
complex spiral_noiser(double r) {
   int N = 10;
   double sum = 0;
   for (int i = 1; i <= N; i++) {
      sum += sin(r) * sin(i * r) + cos(r) * cos(i * r);
   }
   return sum * 1 / sqrt(N);
}

// Composes the spiral noiser with the spiral function.
complex noise_composer(complex z) {
   return cabs(z) * cpow(M_E, I * (spiral_noiser(cabs(z)) + carg(z)));
   // return creal(z) + I * (cimag(z) + sin(creal(z)));
}

// Tie dye spiral function.
complex tie_dye_spiral(complex z) {
   double c = 3;
   // return cpow(z, 1 + c * I);
   double r = cabs(z);
   double t = carg(z);
   return r * cpow(M_E, -c * t) * cpow(M_E, I * (c * log(r) + t));
}

// Simple tetration function.
complex tetrate(complex z, int n) {
   complex res = z;
   for (int i = 0; i < n; i++) {
      res = cpow(z, res);
   }
   return res;
}

// Riemann zeta function.
complex zeta(complex z) {
   complex res = 0;
   for (int i = 1; i <= ZETA_ITER; i++) {
      res += 1 / cpow(i, z);
   }
   return res;
}

// Coefficients for the Eisenstein series.
complex eisen_coeff(int n, int k) {
   int res = 0;
   for (int m = 1; m <= n; m++) {
      if (n % m == 0) {
         res += pow(m, k - 1);
      }
   }
   return res;
}

// Returns the weight k Eisenstein series.
complex eisenstein(complex z, int k) {
   complex res = 0;
   complex q = cpow(M_E, 2 * M_PI * I * z);
   for (int n = 1; n < ZETA_ITER; n++) {
      res += eisen_coeff(n, k) * cpow(q, n);
   }

   return 2 * zeta(k) + 2 * cpow(2 * M_PI * I, k) / factorial(k - 1) * res;
}

// Binet function.
complex binet(complex z) {
   complex phi = (1 + sqrt(5)) / 2;
   return (cpow(phi, z) - ccos(z * M_PI) * cpow(phi, -z)) / sqrt(5);
}

// Weierstrass sigma function.
complex weierstrass_sigma(complex z) {
   int up_bound = 50;
   complex tau = I + 1;
   complex q = cpow(M_E, 2 * M_PI * I * tau);
   complex prod = 1;

   for (int n = 1; n < up_bound; n++) {
      prod *= (1 - cpow(M_E, 2 * M_PI * I * z) * cpow(q, n));
      prod *= (1 - cpow(M_E, 2 * M_PI * I * -z) * cpow(q, n));
      prod /= cpow(1 - cpow(q, n), 2);
   }

   prod /= (2 * M_PI * I);
   prod *= cpow(M_E, 0.5 * eisenstein(tau, 2) * cpow(z, 2));
   prod *= (cpow(M_E, M_PI * I * z) - cpow(M_E, -M_PI * I * z));

   return prod;
}

// Complex valued function.
complex f(complex x) {
   // return cpow(x, 3) - 1;
   // return 3 * cpow(x, 3) - 2;
   // return (cpow(x, 2) - 1) * cpow((x - 2 - I), 2) / (cpow(x, 2) + 2 + 2 * I);
   // return creal(x) + cimag(x) + creal(x) * cimag(x);
   // return x;

   // complex w = -sqrt(2) + sqrt(2) * I;
   // return (w - x) / (1 - conj(w) * x);

   // return b(x);

   // return cpow(x, I);

   // return cpow(x, 7) + cexp(x);

   // return m(x);

   // return wp_prime(x);

   // return tie_dye_spiral(x);

   // return tie_dye_spiral(noise_composer(x));

   // return noise_composer(tie_dye_spiral(x));

   // return m(b(x));

   // return tetrate(x, 40);

   // return zeta(x);

   // return eisenstein(x, 4);

   // return b(zeta(x));

   // return -1 / x;

   // return binet(x);

   return weierstrass_sigma(x);
}

void get_hls(complex z, tHLS *hls) {
   double t, hue, m, power, r0, r1, r, lightness, saturation;

   if (z == 0) {
      hls->h = 0;
      hls->l = 0;
      hls->s = 0;
      return;
   }

   // Calculate hue based on argument (angle) of input.
   t = carg(z);
   if (t < 0) {
      t += 2.0 * M_PI;
   }
   hls->h = t;

   // Calculate lightness based on absolute value of input.
   m = cabs(z);
   power = log(m) / log(CONTOUR_BASE);
   r0 = pow(CONTOUR_BASE, floor(power));
   r1 = pow(CONTOUR_BASE, ceil(power));
   if (r1 == r0) {
      hls->l = 0.6;
   }
   else {
      r = (m - r0) / (r1 - r0);
      hls->l = 0.3 * r + 0.3;
   }

   // Set saturation to 1.
   hls->s = 1.0;
}

// Converts HLS color coordinates into RGB.
void hls_to_rgb(tHLS *hls, tRGB *rgb) {
   double h, l, s, hp, c, x, m, r0, g0, b0;
   h = hls->h;
   l = hls->l;
   s = hls->s;

   // Custom scale for the h coordinate of the color map.
   hp = h * 3 / M_PI;

   // Color offset constants.
   c = (1 - fabs(2 * l - 1)) * s;
   x = c * (1 - fabs(fmod(hp, 2) - 1));
   m = l - c / 2;

   // Cases for conversion.
   if (hp >= 0 && hp <= 1) {
      r0 = c;
      g0 = x;
      b0 = 0;
   }
   else if (hp >= 1 && hp <= 2) {
      r0 = x;
      g0 = c;
      b0 = 0;
   }
   else if (hp >= 2 && hp <= 3) {
      r0 = 0;
      g0 = c;
      b0 = x;
   }
   else if (hp >= 3 && hp <= 4) {
      r0 = 0;
      g0 = x;
      b0 = c;
   }
   else if (hp >= 4 && hp <= 5) {
      r0 = x;
      g0 = 0;
      b0 = c;
   }
   else if (hp >= 5 && hp <= 6) {
      r0 = c;
      g0 = 0;
      b0 = x;
   }
   else {
      r0 = 0;
      g0 = 0;
      b0 = 0;
   }

   rgb->r = (int) ((r0 + m) * 255);
   rgb->g = (int) ((g0 + m) * 255);
   rgb->b = (int) ((b0 + m) * 255);
}

// Create domain of numbers.
void create_domain(double minx, double maxx, int sample, double* domain) {
   double step  = (maxx - minx) / sample;
   for (int i = 0; i < sample; i++) {
      domain[i] = minx + i * step;
   }
}

// Create domain of complex numbers.
void create_complex_domain(int sample, double *domain, complex *cdomain) {
   for (int i = 0; i < sample; i++) {
      for (int j = 0; j < sample; j++) {
         cdomain[i + sample * j] = domain[i] + I * domain[j];
      }
   }
}

// Return a random sequence of numbers for a random Blaschke product.
void random_blaschke() {
   time_t t;
   double a, b;

   seq = (complex *) malloc(sizeof(complex) * SEQ_LEN);
   srand((unsigned) time(&t));

   for (int i = 0; i < SEQ_LEN; i++) {
      a = 2 * ((double) rand() - (double) RAND_MAX / 2) / (double) RAND_MAX;
      b = 2 * ((double) rand() - (double) RAND_MAX / 2) / (double) RAND_MAX;
      seq[i] = a + b * I;
   }
}

// Draw image to a .ppm file.
void draw_to_file(char *path_name) {
   FILE *fp;
   double *domain = (double *) malloc(sizeof(double) * SAMPLE);
   complex *cdomain = (complex *) malloc(sizeof(complex) * SAMPLE * SAMPLE);
   tRGB *colors = (tRGB *) malloc(sizeof(tRGB) * SAMPLE * SAMPLE);

   // Call Blaschke sequence generator (if using Blaschke function)
   random_blaschke();

   create_domain(MIN_X, MAX_X, SAMPLE, domain);
   create_complex_domain(SAMPLE, domain, cdomain);

   tHLS hls;
   tRGB rgb;

   for (int i = 0; i < SAMPLE * SAMPLE; i++) {
      int is_grid_point = 0;

      // Draw grid (optional)
      // double eps = 0.00001;
      // for (int j = 0; j < MAX_X; j++) {
      //    if (fabs(creal(cdomain[i]) - j) < eps
      //       || fabs(cimag(cdomain[i]) - j) < eps
      //       || fabs(creal(cdomain[i]) + j) < eps
      //       || fabs(cimag(cdomain[i]) + j) < eps) {
      //       get_hls(0, &hls);
      //       is_grid_point = 1;
      //    }
      // }

      if (!is_grid_point) {
         get_hls(f(cdomain[i]), &hls);
      }
      hls_to_rgb(&hls, &rgb);
      colors[i] = rgb;
   }

   fp = fopen(path_name, "w");
   fprintf(fp, "P3\n");
   fprintf(fp, "%d %d\n", SAMPLE, SAMPLE);
   fprintf(fp, "255\n");
   for (int i = SAMPLE - 1; i >= 0; i--) {
      for (int j = 0; j < SAMPLE; j++) {
         rgb = colors[SAMPLE * i + j];
         fprintf(fp, "%d %d %d ", rgb.r, rgb.g, rgb.b);
      }
      fprintf(fp, "\n");
   }

   fclose(fp);

   free(domain);
   free(cdomain);
   free(colors);
   free(seq);
}

// Generates a series of images to be converted to a GIF for a Mandelbrot zoom.
void generate_man_zoom() {
   double zoom_factor = 2;
   for (int i = 501; i <= 600; i++) {
      man_magnify = zoom_factor * pow(1.05, i);
      char path[25];
      snprintf(path, 25, "./man_zoom/pic%d.ppm", i);
      draw_to_file(path);
   }
}

// Main function.
int main() {
   // Initialize random number generator
   // time_t t;
   // srand((unsigned) time(&t));

   // man_magnify = 2 * pow(1.2, 38);
   draw_to_file("dc.ppm");
   // generate_man_zoom();

   return 0;
}
