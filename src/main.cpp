#include <cstdlib>
#include <iostream>
#include <armadillo>
#include <math.h>
#include <complex>

arma::cx_vec inverse_fourier_1d(arma::vec q, std::string myfunction, double l);
void print_fr_1d(arma::cx_vec fr, std::string myfunction, double l);

int main(int argc, char *argv[]){
    if(argc!=4){
        std::cout << "Error: require system size, length scale (lambda), and dimension (1 or 3)." << std::endl;
        exit(0);
    }
    int nx = std::stoi(argv[1]);
    double l = std::stod(argv[2]);
    int dim = std::stoi(argv[3]);
    std::string myfunction = "inverse-exp";

    if(dim==1){
        arma::vec q;
        q.zeros(nx);
        for(int i=0; i<nx; i++) q(i) = 2*M_PI*i/nx;
        arma::cx_vec fr = inverse_fourier_1d(q, myfunction, l);
        print_fr_1d(fr, myfunction, l);
    }

    return 0;
}

arma::cx_vec inverse_fourier_1d(arma::vec q, std::string myfunction, double l){

    using namespace std::complex_literals; //to get 1i

    int nx = q.size();//int(arma::size(q));
    std::cout << "nx: " << nx << std::endl;
    arma::vec fq;
    arma::cx_vec fr;
    fq.zeros(nx);
    fr.zeros(nx);
    if(myfunction=="inverse-exp"){
        for(int i=0; i<nx/2+1; i++){
            double q_sq = q(i)*q(i);
            fq(i) = nx/(1+l*l*q_sq);
            if(i>0) fq(nx-i) = fq(i);
            std:: cout << fq(i) << std::endl;
        }
    }
    for(int i=0; i<nx; i++){
        for(int k1=0; k1<nx; k1++){
            std::complex<double> update = fq(k1)*std::exp(2*M_PI*1i*(1.0*i*k1/nx));
            fr(i) += update;
        }
    }
    std::cout << fr << std::endl;
    return fr;
}

void print_fr_1d(arma::cx_vec fr, std::string myfunction, double l){
    int nx = fr.size();
    std::string output = myfunction + "_nx=" + std::to_string(nx) + "_l=" + std::to_string(l) + ".txt";
    std::ofstream ofile;
    ofile.open(output);
    for(int i=0; i<nx; i++){
        ofile << 1.0*i << " " << std::real(fr(i)) << std::endl;
    }
}