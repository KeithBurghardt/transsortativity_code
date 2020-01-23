/* Calculate the probability of majority in a zero-one network */
/* P(x=1) is specified and assigned randomly initially         */
/* then swaped to higher xk correlation                        */
/* Add step width, space optimization */
/* Compile by g++ -O3 zeroone_opt.c -o zeroone */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include "factorial.h"

int sum_array(int a[], int n);
double stdev_array(int a[], int n);
int compare (const void *a, const void *b);

int main(int argc, char *argv[]) {
  srand(time(NULL));
  int n = 10000; /* Number of nodes in the network */
  double x_prob = 0.01; /* P(x=1) */
  int step = 5;

  bool **A = new bool*[n];
  for (int i=0; i<n; i++) A[i] = new bool[n];
  for (int i=0; i<n; i++) memset(A[i], 0, sizeof(bool)*n);
  int i, j;
  // infile
  string infile = "../Preprocess/Sets24corr/24mn07_c30";
  // outfile
  string of = "zeroone.out";
  FILE *edgelist;
  edgelist = fopen(infile, "r"); /* Input file */
  char line[16];
  while (fgets(line, 16, edgelist))
    if (line[0] != '#' && sscanf(line, "%d %d", &i, &j) == 2) {
      A[i][j]++;
      A[j][i]++;
    }
  fclose(edgelist);
  FILE *outfile;
  outfile = fopen(of, "w"); /* Output file */

  /* assortativity */
  double Etotal = 0.0, c = 0.0, s = 0.0;
  int *k = new int[n];
  memset(k, 0.0, sizeof(int)*n);
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (A[i][j]) k[i]++;
  double k_avg = (double) sum_array(k, n)/n;
  double k_std = stdev_array(k, n);
  int *k_sorted = new int[n];
  memcpy(k_sorted, k, sizeof(int)*n);
  qsort(k_sorted, n, sizeof(int), compare);
  int k_max = k_sorted[0];

  int *k_hist = new int[k_max+1];
  memset(k_hist, 0, sizeof(int)*(k_max+1));
  for (int i=0; i<n; i++) k_hist[k[i]]++;
  double **E = new double*[k_max+1];
  for (int i=0; i<k_max+1; i++) E[i] = new double[k_max+1];
  for (int i=0; i<k_max+1; i++) memset(E[i], 0, sizeof(double)*(k_max+1));
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (A[i][j]) E[k[i]][k[j]]++;
  double *Q = new double[k_max+1];
  memset(Q, 0.0, sizeof(double)*(k_max+1));
  for (int i=0; i<k_max+1; i++) for (int j=0; j<k_max+1; j++) Q[i] += E[i][j];
  for (int i=0; i<k_max+1; i++) Etotal += Q[i];
  for (int i=0; i<k_max+1; i++) for (int j=0; j<k_max+1; j++) E[i][j] /= Etotal;
  for (int i=0; i<k_max+1; i++) Q[i] /= Etotal;
  for (int i=0; i<k_max+1; i++) for (int j=0; j<k_max+1; j++) {
    c += i*j*(E[i][j] - Q[i]*Q[j]);
    if (i == j) s += i*j*(Q[i] - Q[i]*Q[j]);
    else s -= i*j*Q[i]*Q[j];
  }
  double r = c/s;
  printf("Network assortativity: %f\n", r);
  fprintf(outfile, "Network assortativity: %f\n", r);
  fprintf(outfile, "Max Degree: %d\n", k_max);

  /* k-x correlation */
  double x_avg = x_prob;
  double x_std = sqrt(x_prob*(1-x_prob));
  int *x = new int[n];
  memset(x, 0, sizeof(int)*n);
  while (sum_array(x, n) < x_prob*n) x[rand()%n] = 1;
  double **p_xk = new double*[2];
  for (int i=0; i<2; i++) p_xk[i] = new double[n];
  for (int i=0; i<n; i++) p_xk[x[i]][k[i]] += 1.0/n;
  int xk_sum = 0;
  for (int i=0; i<n; i++) if (x[i]) xk_sum += k[i];
  double xk_avg = (double) xk_sum/n;
  double corr = (xk_avg - x_avg * k_avg)/x_std/k_std;

  int target = sum_array(k_sorted, (int)(x_avg * n));
  double kx_corr = 0.0;

  /* swap x=1 to higher degrees to achieve higher rho_kx */
  while (xk_sum < target) {
    int u, v;
    for (int i=0; i<step; i++) {
      while (true) {
        u = rand()%n;
        v = rand()%n;
        if (x[u] && !x[v] && k[v] > k[u]) break;
      }
      x[u]--;
      x[v]++;
      p_xk[1][k[u]] -= 1.0/n;
      p_xk[0][k[u]] += 1.0/n;
      p_xk[0][k[v]] -= 1.0/n;
      p_xk[1][k[v]] += 1.0/n;
      xk_sum += k[v] - k[u];
      if (xk_sum >= target) break;
    }
    xk_avg = (double) xk_sum/n;
    corr = (xk_avg - x_avg * k_avg)/x_std/k_std;
    printf("One edges: %d/%d, Correlation: %f\n", xk_sum, target, corr);

    int *nn_active = new int[n];
    memset(nn_active, 0, sizeof(int)*n);
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (A[i][j] && x[j]) nn_active[i]++;
    int paradox = 0;
    for (int i=0; i<n; i++) if (nn_active[i] > (double) k[i]/2.0) paradox++;

    double *P_f = new double[k_max+1];
    memset(P_f, 0.0, sizeof(double)*(k_max+1));
    for (int i=0; i<k_max+1; i++) for (int j=0; j<k_max+1; j++)
      if (Q[i] != 0 && p_xk[0][j] + p_xk[1][j] != 0) 
        P_f[i] += p_xk[1][j]/(p_xk[0][j] + p_xk[1][j])*E[i][j]/Q[i];
    //if (corr > kx_corr) {
    //  for (int i=0; i<k_max+1; i++) fprintf(outfile, "%f, ", P_f[i]);
    //  fprintf(outfile, "\n");
    //  kx_corr += 0.2;
    //}
    
#if 1
    double **P_ff = new double*[k_max+1];
    for (int i=0; i<k_max+1; i++) P_ff[i] = new double[3];
    for (int i=0; i<k_max+1; i++) memset(P_ff[i], 0, sizeof(double)*3);
    for (int i=0; i<n; i++) {
      double ph = (double) nn_active[i]/ k[i];
      P_ff[k[i]][0] += ph / (double) k_hist[k[i]];
      P_ff[k[i]][1] += ph * ph / (double) k_hist[k[i]];
    }
    for (int i=0; i<k_max+1; i++) P_ff[i][2] = sqrt(P_ff[i][1] - P_ff[i][0] * P_ff[i][0]);
#endif

#if 1
    double P_maj = 0.0;
    P_maj += (double) k_hist[1]/n * P_ff[1][0];
    P_maj += (double) k_hist[2]/n * (2*P_ff[2][0]*P_ff[2][0] + 2*P_ff[2][2]*P_ff[2][2] - P_ff[2][0]);
    for (int i=3; i<k_max+1; i++) if (P_ff[i][0] > 0.0 && P_ff[i][0] < 1.0) {
		double z;
		if (i%2) // to avoid half weight of middle bin problem
		  z = (0.5 - P_ff[i][0]) / P_ff[i][2];
		else
		  z = (0.5 + 0.5/i - P_ff[i][0]) / P_ff[i][2];
        //double z = (0.5-P_f[i]) / sqrt(P_f[i]*(1-P_f[i])) * sqrt(i);
        //double z = (0.5+0.5/(i+1)-P_ff[i][0]) / sqrt(P_ff[i][0]*(1-P_ff[i][0])) * sqrt(i);
        P_maj += (double) k_hist[i]/n * (0.5 - 0.5 * erf(z/sqrt(2))); /* Normal CDF tail integral */
    }
#else
    double P_maj = 0.0, log_prob;
    for (int i=1; i<k_max+1; i++) for (int j = floor(i/2)+1; j<i+1; j++)
      if (P_f[i] > 0.0 && P_f[i] < 1.0) {
        log_prob = log_fac[i] - log_fac[j] - log_fac[i-j] + j * log10(P_f[i]) + (i - j) * log10(1.0 - P_f[i]);
        P_maj += (double) k_hist[i]/n * pow(10,log_prob);
    }
#endif
    printf("Paradox: %d/%d, Prob of Maj: %f\n", paradox, n, P_maj);
    fprintf(outfile, "[%f, %d, %f], ", corr, paradox, P_maj);
  }
  fprintf(outfile, "\n");
  fclose(outfile);
  return 0;
}

int sum_array(int a[], int n) {
  int sum = 0;
  for (int i=0; i<n; i++) sum += a[i];
  return sum;
}

double stdev_array(int a[], int n) {
  double sqsum = 0.0, sum = 0.0;
  for (int i=0; i<n; i++) sum += a[i];
  double mu = sum/n;
  for (int i=0; i<n; i++) sqsum += a[i]*a[i];
  return sqrt(sqsum/n - mu * mu);
}

int compare (const void *a, const void *b) {
  return *(int*)b - *(int*)a;
}

double assortativity(int n, int **A) {
  double Etotal = 0.0, c = 0.0, s = 0.0;
  int *k = new int[n];
  memset(k, 0, sizeof(int)*n);
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) k[i] += A[i][j];
  double **E = new double*[n];
  for (int i=0; i<n; i++) E[i] = new double[n];
  for (int i=0; i<n; i++) memset(E[i], 0, sizeof(double)*n);
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (A[i][j]) E[k[i]][k[j]]++;
  double *Q = new double[n];
  memset(Q, 0, sizeof(double)*n);
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) Q[i] += E[i][j];
  for (int i=0; i<n; i++) Etotal += Q[i];
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) E[i][j] /= Etotal;
  for (int i=0; i<n; i++) Q[i] /= Etotal;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    c += i*j*(E[i][j] - Q[i]*Q[j]);
    if (i == j) s += i*j*(Q[i] - Q[i]*Q[j]);
    else s -= i*j*Q[i]*Q[j];
  }
  return c/s;
}
