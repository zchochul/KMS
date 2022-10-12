#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;
//input n, a, f, R, L
int main() {
    //std::cout << "Hello World!\n";
    //opening files
    ofstream outputfile;
    ifstream inputfile;
    ofstream outputfile2; //saving momentum
    
    outputfile.open("Lab1dane.txt");
    inputfile.open("Lab1parametry.txt");
    outputfile2.open("Lab1Momentum.txt");

    if(!inputfile) {
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    //reading parameters
    double param[5]; //variable for parameters
    for (int i=0; i<=5; i++){
        inputfile >> param[i];
    }
    
    //asigning values to parameters
    double n = param[0]; 
    double a = param[1];
    double f = param[2];
    double R = param[3];
    double L = param[4];
    double m = param[5];
    double T_0 = 100;
    
    //x,y,z & p
    double b0[3] = {a, 0, 0};
    double b1[3] = {a/2, a*std::pow(3, 1/2.)/2, 0};
    double b2[3] = {a/2, a*std::pow(3, 1/2.)/6, a*std::pow(0.6666666666, 1/2.)};
    
    //cout << "b0 = " << b0[0] << "\t" << b0[1] << "\t" << b0[2] << endl;
    //cout << "b1 = " << b1[0] << "\t" << b1[1] << "\t" << b1[2] << endl;
    //cout << "b2 = " << b2[0] << "\t" << b2[1] << "\t" << b2[2] << endl;
    
    //double Ekin_t0[number_of_particles][3];
    /* initialize random seed: */
    srand (time(NULL));
    int number_of_particles = pow(n, 3); //Number of particles == N
    double r_t0[number_of_particles][3];
    double p_t0[number_of_particles][3];
    double Ekinx_t0, Ekiny_t0, Ekinz_t0;
    double lambdax, lambday, lambdaz;
    double signx, signy, signz;
    double k = 8.31e-3 ;//* pow(10, -23) ; //1.38 * 10^-23
    double Ekin_constant = - 0.5 * k * T_0;

    outputfile << number_of_particles << endl << endl;
    for(int i0 = 0; i0<= n-1; i0++){
        for(int i1 = 0; i1<= n-1; i1++){
            for(int i2 = 0; i2<= n-1; i2++){
                int particle_number = i0 + i1*n + i2*pow(n,2);
                double s1 = i0 - (n-1)/2;
                double s2 = i1 - (n-1)/2;
                double s3 = i2 - (n-1)/2;
                // x, y, z
                r_t0[particle_number][0] = (s1)*b0[0] + (s2)*b1[0] + (s3)*b2[0];
                r_t0[particle_number][1] = (s1)*b0[1] + (s2)*b1[1] + (s3)*b2[1];
                r_t0[particle_number][2] = (s1)*b0[2] + (s2)*b1[2] + (s3)*b2[2];
                outputfile << "Ar" << "\t" << r_t0[particle_number][0] << "\t" << r_t0[particle_number][1] << "\t" << r_t0[particle_number][2] << endl;

                //momentum
                //different lambda for different coordinate
                lambdax = (float) rand()/RAND_MAX;
                lambday = (float) rand()/RAND_MAX;
                lambdaz = (float) rand()/RAND_MAX;

                Ekinx_t0 = Ekin_constant * log(lambdax);
                Ekiny_t0 = Ekin_constant * log(lambday);
                Ekinz_t0 = Ekin_constant * log(lambdaz);

                double p_const = pow(2*m, 0.5);
                //cout << Ekin_t0 [particle_number] << endl;
                signx = (float) rand()/RAND_MAX;
                signy = (float) rand()/RAND_MAX;
                signz = (float) rand()/RAND_MAX;

                if(signx < 0.5) signx = -1;
                else signx = 1;
                if(signy < 0.5) signy = -1;
                else signy = 1;
                if(signz < 0.5) signz = -1;
                else signz = 1;

                //cout << "E_kin = " << pow(Ekinx_t0, 0.5) << endl;
                //cout << "E_kin = " << Ekinx_t0 << endl;
                p_t0[particle_number][0] = signx * p_const* pow(abs((double)Ekinx_t0), 0.5);
                p_t0[particle_number][1] = signy * p_const* pow(abs((double)Ekiny_t0), 0.5);
                p_t0[particle_number][2] = signz * p_const* pow(abs((double)Ekinz_t0), 0.5);
                
                double Px = 0;
                double Py = 0;
                double Pz = 0;

                for(int i=0; i<number_of_particles; i++){
                    Px +=p_t0[i][0];
                    Py +=p_t0[i][1];
                    Pz +=p_t0[i][2];
                }
                
                //cout << "P_x = " << Px << endl;

                //cout << "p = " << p_t0[particle_number][0]  -1/number_of_particles*Px <<  "\t" << p_t0[particle_number][1] -1/number_of_particles*Py<< "\t" << p_t0[particle_number][0] -1/number_of_particles*Pz<< endl;
                p_t0[particle_number][0] = p_t0[particle_number][0] -1/number_of_particles*Px;
                outputfile2 << p_t0[particle_number][0] << endl;
                outputfile2 << p_t0[particle_number][1] << endl;
                outputfile2 << p_t0[particle_number][2] << endl;
            }
        }
    }
    
    for (int i = 0; i< number_of_particles ; i++){

    }
    
    outputfile2.close();
    //histogram

    //closing files
    outputfile.close();
    inputfile.close();
    std::cout << "Finished!\n";
    return 0;
}