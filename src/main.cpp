#include <cstdlib>
#include <iostream>
#include <armadillo>
#include <math.h>
#include <complex>

arma::cx_vec inverse_fourier_1d(arma::vec q, std::string myfunction, double l);
void print_fr_1d(arma::cx_vec fr, std::string myfunction, double l);
arma::cx_cube inverse_fourier_3d(int nx, int ny, int nz, std::string myfunction, double l);
void print_fr_3d(int nx, arma::cx_cube fr, std::string myfunction, double l);

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
    else if(dim==3){
        arma::cx_cube fr = inverse_fourier_3d(nx, nx, nx, myfunction, l);
        print_fr_3d(nx, fr, myfunction, l);
    }
    else{
        std::cout << "Dimension not supported." << std::endl;
        exit(0);
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
        //Compute inverse FT of exponential
        for(int i=0; i<nx/2+1; i++){
            double q_sq = q(i)*q(i);
            fq(i) = nx/(1+l*l*q_sq);
            if(i>0) fq(nx-i) = fq(i);
            std:: cout << fq(i) << std::endl;
        }
    }
    //Get FT
    for(int i=0; i<nx; i++){
        for(int k1=0; k1<nx; k1++){
            std::complex<double> update = fq(k1)*std::exp(2*M_PI*1i*(1.0*i*k1/nx));
            fr(i) += update;
        }
    }
    std::cout << fr << std::endl;
    return fr;
}

arma::cx_cube inverse_fourier_3d(int nx, int ny, int nz, std::string myfunction, double l){

    using namespace std::complex_literals; //to get 1i

    arma::cube fq;
    arma::cx_cube fr;
    fq.zeros(nx, ny, nz);
    fr.zeros(nx, ny, nz);
    arma::cube is_filled(nx, ny, nz, arma::fill::zeros);
    if(myfunction=="inverse-exp"){
        //Compute inverse FT of exponential
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                for(int k=0; k<nz; k++){
                    if(is_filled((nx-i)%nx,(ny-j)%ny,(nz-k)%nz)==1){                    
                        fq(i,j,k) = fq((nx-i)%nx,(ny-j)%ny,(nz-k)%nz);
                        is_filled(i,j,k)=1.0;
                    }
                    else{
                        double q_sq = 4*M_PI*M_PI*(i*i/(1.0*nx*nx) + j*j/(1.0*ny*ny) + k*k/(1.0*nz*nz));
                        fq(i,j,k) = 1.0*nx*ny*nz/pow(1+l*l*q_sq, 2.0);
                        //std::cout << q_sq << std::endl;
                        //std::cout << (1+l*l*q_sq)*(1+l*l*q_sq) << std::endl;
                        is_filled(i,j,k)=1.0;
                    }
                }
            }
        }
    }
    //Get FT
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                for(int q1=0; q1<nx; q1++){
                    for(int q2=0; q2<ny; q2++){
                        for(int q3=0; q3<nz; q3++){
                            std::complex<double> update = fq(q1,q2,q3)*std::exp(2*M_PI*1i*(
                                (1.0*i*q1)/(1.0*nx) +
                                (1.0*j*q2)/(1.0*ny) +
                                (1.0*k*q3)/(1.0*nz)));
                            //Check that imaginary part is zero
                            //if(update.imag()>1e-3) std::cout << "imaginary part:" << update.imag() << std::endl;
                            fr(i,j,k) += update;
                        }
                    }
                }
                fr(i,j,k) *= 1.0/(nx*ny*nz);
            }
        }
    }

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

void print_fr_3d(int nx, arma::cx_cube fr, std::string myfunction, double l){
    std::string output = myfunction + "_3d_nx=" + std::to_string(nx) + "_l=" + std::to_string(l) + ".txt";
    std::ofstream ofile;
    ofile.open(output);
    for(int i=0; i<nx; i++){
        for(int j=0; j<nx; j++){
            for(int k=0; k<nx; k++){
                ofile << 1.0*i << " " << 1.0*j << " " << 1.0*k << " " << std::real(fr(i,j,k)) << std::endl;
            }
        }
    }
}