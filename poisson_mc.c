#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

#define N 30
#define STEPS 100000000
#define EPS 0.0000001  // Stddev for Gaussian update
#define BETA 1000
#define MAX_SAMPLES (STEPS/1000000)
// It is interesting to see how EPS must scale with teh larttice size

double energies[MAX_SAMPLES]; 
double f[N][N], u[N][N], ux[N][N], uy[N][N];
double phi_x[N][N], phi_y[N][N];
double best_phi_x[N][N], best_phi_y[N][N];
double best_energy = 1e20;
int best_iter = -1;

gsl_rng *rng;

double compute_max_norm() {
    double maxv = 0.0;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        double vx = ux[i][j] + phi_x[i][j];
        double vy = uy[i][j] + phi_y[i][j];
        double mag = hypot(vx, vy);
        if (mag > maxv) maxv = mag;
    }
    return maxv;
}

void solve_poisson() {
    fftw_complex *fhat = fftw_malloc(sizeof(fftw_complex)*N*(N/2+1));
    fftw_complex *uhat = fftw_malloc(sizeof(fftw_complex)*N*(N/2+1));
    double *in = fftw_malloc(sizeof(double)*N*N);
    double *out = fftw_malloc(sizeof(double)*N*N);
    fftw_plan pf = fftw_plan_dft_r2c_2d(N,N,in,fhat,FFTW_ESTIMATE);
    fftw_plan pb = fftw_plan_dft_c2r_2d(N,N,uhat,out,FFTW_ESTIMATE);

    for (int i=0;i<N;i++)for(int j=0;j<N;j++) in[i*N+j]=f[i][j];
    fftw_execute(pf);

    for(int i=0;i<N;i++) {
        int kx = (i<=N/2)?i:i-N;
        for(int j=0;j<=N/2;j++){
            int ky=j;
            int idx = i*(N/2+1)+j;
            if(kx==0&&ky==0) uhat[idx][0]=uhat[idx][1]=0.0;
            else {
                double den = -4*M_PI*M_PI*(kx*kx+ky*ky);
                uhat[idx][0] = fhat[idx][0]/den;
                uhat[idx][1] = fhat[idx][1]/den;
            }
        }
    }

    fftw_execute(pb);
    for (int i=0;i<N;i++) for(int j=0;j<N;j++)
         u[i][j] = out[i*N+j]/(N*N);

    fftw_destroy_plan(pf);
    fftw_destroy_plan(pb);
    fftw_free(fhat); fftw_free(uhat);
    fftw_free(in); fftw_free(out);
}

void compute_grad() {
    for (int i=0;i<N;i++)for(int j=0;j<N;j++){
        int ip=(i+1)%N, jp=(j+1)%N;
        ux[i][j] = u[ip][j] - u[i][j];
        uy[i][j] = u[i][jp] - u[i][j];
    }
}

void step_mh(int iter) {
    int i = gsl_rng_uniform_int(rng, N);
    int j = gsl_rng_uniform_int(rng, N);
    double E_old = compute_max_norm();
    double delta = gsl_ran_gaussian(rng, EPS);

    phi_x[i][j] += delta;
    phi_x[i][(j+1)%N] -= delta;
    phi_y[i][j] -= delta;
    phi_y[(i+1)%N][j] += delta;

    double E = compute_max_norm();

    if (E < best_energy) {
        best_energy = E;
        best_iter = iter;
        // copy best
        for(int a=0;a<N;a++)for(int b=0;b<N;b++){
            best_phi_x[a][b]=phi_x[a][b];
            best_phi_y[a][b]=phi_y[a][b];
        }
    } else {
        double dE = E - E_old;
        if (gsl_rng_uniform(rng) > exp(-BETA*dE/2)) {
            // revert
            phi_x[i][j] -= delta;
            phi_x[i][(j+1)%N] += delta;
            phi_y[i][j] += delta;
            phi_y[(i+1)%N][j] -= delta;
        }
    }
}

void save_csv(const char* file, double mat[N][N]) {
    FILE *f=fopen(file,"w");
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            fprintf(f,"%.8f", mat[i][j]);
            if (j+1<N) fputc(',',f);
        }
        fputc('\n',f);
    }
    fclose(f);
}

int main(int argc, char** argv){
int sample_freq = 1000000;  // sample every 100 steps
int sample_idx = 0;

    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));

    // Example: delta mass in center
    for(int i=0;i<N;i++)for(int j=0;j<N;j++)
        f[i][j] = (i==N/2 && j==N/2) ? 1.0 : 0.0;

    solve_poisson();
    compute_grad();

    // Initialize phi to zero
    for(int i=0;i<N;i++)for(int j=0;j<N;j++)
        phi_x[i][j]=phi_y[i][j]=0.0;

    for(int it=0;it<STEPS;it++){
        step_mh(it);
        
            if (it % sample_freq == 0 && sample_idx < MAX_SAMPLES) {
    	    energies[sample_idx++] = compute_max_norm();
    }   
            if (it % 10000 == 0)
            printf("Iter %d, best ||∇u+Φ||∞ = %.6f at %d\n",
                   it, best_energy, best_iter);
    }
    printf("Done. Best energy=%.6f at iter %d\n", best_energy, best_iter);
    
    FILE *fe = fopen("energies.csv", "w");
for(int i=0; i<sample_idx; i++) {
    fprintf(fe, "%.8f\n", energies[i]);
}
fclose(fe);

    // Save best results
    save_csv("u.csv", u);
    save_csv("ux.csv", ux);
    save_csv("uy.csv", uy);
    save_csv("phi_x.csv", best_phi_x);
    save_csv("phi_y.csv", best_phi_y);

    gsl_rng_free(rng);
    return 0;
}

