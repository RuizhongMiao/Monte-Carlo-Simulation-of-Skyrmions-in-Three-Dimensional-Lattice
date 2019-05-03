#ifndef SkyrmionM_lattice_h
#define SkyrmionM_lattice_h

#include <iostream>
#include <fstream>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#define BOLTZMANN_K 1.0

//\vec{a} \cdot \vec{b}
double innerProduct(double a[3], double b[3]){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

//\vec{a} \times \vec{b} \cdot \hat{n}
double mixtureProduct(double a[3], double b[3], char direction){
    double result = 0;
    switch (direction) {
        case 'x':
            result = a[1] * b[2] - a[2] * b[1];
            break;
        case 'y':
            result = a[2] * b[0] - a[0] * b[2];
            break;
        case 'z':
            result = a[0] * b[1] - a[1] * b[0];
            break;
        default:
            break;
    }
    return result;
}

class lattice{
private:
    const int l_x ;                   //length of the lattice
    const int l_y ;                   //width of the lattice
    const int l_z ;                   //height of the lattice
    double ****spin;
    double **B;
    double T = 0.0;
    double J = 1.0;                 //interaction coefficient
    double K = tan(M_PI / 5.0) * J; //DM interaction
    double Jprime = J / 16.0;
    double Kprime = K / 16.0;
    
    double H(int i, int j, int k, double direction[3]){
        double result = 0;
        
        result -= J * (innerProduct(direction, spin[(i + 1) % l_x][j][k]) +
                       innerProduct(direction, spin[(i + l_x - 1) % l_x][j][k]) +
                       innerProduct(direction, spin[i][(j + 1) % l_y][k]) +
                       innerProduct(direction, spin[i][(j + l_y - 1) % l_y][k]) +
                       innerProduct(direction, spin[i][j][(k + 1) % l_z]) +
                       innerProduct(direction, spin[i][j][(k + l_z - 1) % l_z]));
        //result -= B * innerProduct(direction, Bhat);
        result -= B[i][j] * direction[2];
        result -= K * (mixtureProduct(direction, spin[(i + 1) % l_x][j][k], 'x') -
                       mixtureProduct(direction, spin[(i + l_x - 1) % l_x][j][k], 'x') +
                       mixtureProduct(direction, spin[i][(j + 1) % l_y][k], 'y') -
                       mixtureProduct(direction, spin[i][(j + l_y - 1) % l_y][k], 'y') +
                       mixtureProduct(direction, spin[i][j][(k + 1) % l_z], 'z') -
                       mixtureProduct(direction, spin[i][j][(k + l_z - 1) % l_z], 'z'));
        
        return result;
    }
    
    double H(int i, int j, int k){
        return H(i, j, k, spin[i][j][k]);
    }
    
    double Hprime(int i, int j, int k, double direction[3]){
        
        //alternative:
        //omit second nearest-neighbor interaction
        //return 0;
        
        double result = 0;
        
        result += Jprime * (innerProduct(direction, spin[(i + 2) % l_x][j][k]) +
                            innerProduct(direction, spin[(i + l_x - 2) % l_x][j][k]) +
                            innerProduct(direction, spin[i][(j + 2) % l_y][k]) +
                            innerProduct(direction, spin[i][(j + l_y - 2) % l_y][k]) +
                            innerProduct(direction, spin[i][j][(k + 2) % l_z]) +
                            innerProduct(direction, spin[i][j][(k + l_z - 2) % l_z]));
        result += Kprime * (mixtureProduct(direction, spin[(i + 2) % l_x][j][k], 'x') -
                            mixtureProduct(direction, spin[(i + l_x - 2) % l_x][j][k], 'x') +
                            mixtureProduct(direction, spin[i][(j + 2) % l_y][k], 'y') -
                            mixtureProduct(direction, spin[i][(j + l_y - 2) % l_y][k], 'y') +
                            mixtureProduct(direction, spin[i][j][(k + 2) % l_z], 'z') -
                            mixtureProduct(direction, spin[i][j][(k + l_z - 2) % l_z], 'z'));
        
        return result;
    }
    
    double Hprime(int i, int j, int k){
        return Hprime(i, j, k, spin[i][j][k]);
    }
    
public:
    double MCS(){
        return l_x * l_y * l_z;
    }
    
    lattice(int l1, int l2, int l3):l_x(l1), l_y(l2), l_z(l3){
        
        //allocate the memory for spin
        double *rawMem = new double [l_x * l_y * l_z * 3];
        
        double **z = new double *[l_x * l_y * l_z];
        for (int i = 0; i < l_x * l_y * l_z; i++) {
            z[i] = &(rawMem[3*i]);
        }
        
        double ***y = new double **[l_x * l_y];
        for (int i = 0; i < l_x * l_y; i++) {
            y[i] = &(z[l_z * i]);
        }
        
        spin = new double ***[l_x];
        for (int i = 0; i < l_x; i++) {
            spin[i] = &(y[l_y * i]);
        }
        /*for (int i = 0; i < l_x; i++) {
            spin[i] = new double **[l_y];
            for (int j = 0; j < l_y; j++) {
                spin[i][j] = new double *[l_z];
                for (int k = 0; k < l_z; k++) {
                    spin[i][i][k] = &(rawMem[3 * (k + l_z * j + l_y * l_z * i)]);
                }
            }
        }*/
        
        //allocate the memory for B
        B = new double*[l_x];
        for (int i = 0; i < l_x; i++) {
            B[i] = new double[l_y];
        }
    }
    
    lattice &operator=(const lattice &rhs){
        if (this == &rhs) {
            return *this;
        }
        
        this->T = rhs.T;
        this->J = rhs.J;
        this->K = rhs.K;
        
        //allocate the memory for spin
        spin = new double***[l_x];
        for (int i = 0; i < l_x; i++) {
            spin[i] = new double **[l_y];
            for (int j = 0; j < l_y; j++) {
                spin[i][j] = new double*[l_z];
                for (int k = 0; k < l_z; k++) {
                    spin[i][i][k] = new double[3];
                }
            }
        }
        
        //allocate the memory for B
        B = new double*[l_x];
        for (int i = 0; i < l_x; i++) {
            B[i] = new double[l_y];
        }
        
        for (int i = 0; i < l_x; i++) {
            for (int j = 0; j < l_y; j++) {
                for (int k = 0; k < l_z; k++) {
                    this->spin[i][j][k][0] = rhs.spin[i][j][k][0];
                    this->spin[i][j][k][1] = rhs.spin[i][j][k][1];
                    this->spin[i][j][k][2] = rhs.spin[i][j][k][2];
                }
            }
        }
        
        return *this;
    }
    
    void setSpin(char mode){
        double theta, phi;
        switch (mode) {
            case 'r':
                for (int i = 0; i < l_x; i++) {
                    for (int j =0 ; j < l_y; j++) {
                        for (int k = 0; k < l_z; k++) {
                            theta = 2 * M_PI * (double)rand() / (double)RAND_MAX;
                            phi = M_PI * (double)rand() / (double)RAND_MAX;
                            spin[i][j][k][0] = cos(theta) * sin(phi);
                            spin[i][j][k][1] = sin(theta) * sin(phi);
                            spin[i][j][k][2] = cos(phi);
                        }
                    }
                }
                break;
                
            case 'u':
                for (int i = 0; i < l_x; i++) {
                    for (int j =0 ; j < l_y; j++) {
                        for (int k = 0; k < l_z; k++) {
                            spin[i][j][k][0]=0;
                            spin[i][j][k][1]=0;
                            spin[i][j][k][2]=1;
                        }
                    }
                }
                break;
                
            case 'd':
                for (int i = 0; i < l_x; i++) {
                    for (int j =0 ; j < l_y; j++) {
                        for (int k = 0; k < l_z; k++) {
                            spin[i][j][k][0]=0;
                            spin[i][j][k][1]=0;
                            spin[i][j][k][2]=-1;
                        }
                    }
                }
                
                /*for (int i = 10; i < 20; i++) {
                    for (int j =10 ; j < 20; j++) {
                        for (int k = 10; k < 20; k++) {
                            spin[i][j][k][0]=0;
                            spin[i][j][k][1]=0;
                            spin[i][j][k][2]=-1;
                        }
                    }
                }*/
                break;
                
            default:
                break;
        }
    }
    
    void setTemperature(double t){
        T = t;
    }
    
    void setMagneticFieldUniform(double b){
        for (int i = 0 ; i < l_x; i++) {
            for (int j = 0; j < l_y; j++) {
                B[i][j] = b;
            }
        }
    }
    
    //change the magnetic field within a rectangular area
    void setMagneticFieldWindow(int x, int y, int l, int w, double b){
        for (int i = x ; i < x+l; i++) {
            for (int j = y; j < y+w; j++) {
                B[i][j] = b;
            }
        }
    }
    
    void oneFlip(void){
        int i, j, k;
        double theta, phi, p, deltaE;
        double flippedSpin[3];
        
        //pick up an random site
        i = rand() % l_x;
        j = rand() % l_y;
        k = rand() % l_z;
        
        //choose a random direction
        theta = 2 * M_PI * (double)rand() / (double)RAND_MAX;
        phi = M_PI * (double)rand() / (double)RAND_MAX;
        flippedSpin[0] = cos(theta) * sin(phi);
        flippedSpin[1] = sin(theta) * sin(phi);
        flippedSpin[2] = cos(phi);
        
        //energy change
        deltaE = H(i, j, k, flippedSpin) - H(i, j, k) +
                 Hprime(i, j, k, flippedSpin) - Hprime(i, j, k);
        
        if (deltaE < 0) {
            spin[i][j][k][0] = flippedSpin[0];
            spin[i][j][k][1] = flippedSpin[1];
            spin[i][j][k][2] = flippedSpin[2];
        }
        else{
            //flipping probability regarding the energy change
            p = exp(- 1.0 * deltaE / (BOLTZMANN_K * T));
            
            //flip the spin
            if (((double)rand() / (double)RAND_MAX) < p) {
                spin[i][j][k][0] = flippedSpin[0];
                spin[i][j][k][1] = flippedSpin[1];
                spin[i][j][k][2] = flippedSpin[2];
            }
        }
    }
    
    void writeSpin(char *fn){
        std::ofstream fp(fn);
        
        for (int i = 0; i < l_x; i++) {
            for (int j = 0; j < l_y; j++) {
                for (int k = 0; k < l_z; k++) {
                    fp << i << "\t" << j << "\t" << k << "\t" << spin[i][j][k][0] << "\t" << spin[i][j][k][1] << "\t" << spin[i][j][k][2] << "\n";
                }
            }
        }
        
        fp.close();
    }
};

#endif
