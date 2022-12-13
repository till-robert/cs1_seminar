#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <array>
#include <random>
#include <string>


//define constants
// const int dE_edge_to_antiparallel = 2;                      //000 -> 100 or 111 -> 011
// const int dE_edge_to_parallel = -dE_edge_to_antiparallel;   //100 -> 000 or 011 -> 111

// const int dE_to_antiparallel = 4                            //000 -> 010 or 111 -> 101
// const int dE_to_parallel = -dE_to_antiparallel;             //010 -> 000 or 101 -> 111

// const int dE_balanced = 0                                   //100 <-> 110 or 001 <-> 011


//std::random_device rd;  //Will be used to obtain a seed for the random number engine
const int seed = 42;
std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> real_distrib(0.0,1.0);
std::uniform_int_distribution<> bool_distrib(0,1);


struct Spin{ //underlying structure is a boolean, but casting/setting with -1,1,
    bool spin_bool;
    
    inline operator signed char const (){
        return (signed char) spin_bool * 2 - 1;
    }
    // signed char operator=(signed char s){
    //     if(!(s == -1 || s == 1)) throw std::invalid_argument( "received invalid spin value" );
    //     spin_bool = s+1;
    //     return s;
    // }
    inline bool flip(){
        return spin_bool = !spin_bool;
    }
};
class SpinSystem{
public:
    int N_spins;
    bool periodic = false;

    std::vector<Spin> spins;
    std::array<double,5> flip_probabilities;
    std::uniform_int_distribution<>  random_spin_choice_distrib;

    int E_J = 0; //energy without magnetic field contribution: E_tot = E_J - h*M
    int M = 0;
    
    inline static double glauber_prb(double beta,signed char E){
        return exp(-beta*E)/(1+exp(-beta*E));
    }

    
    operator std::string (){
        std::string s = "";
        for(int i = 0; i < N_spins; i++){
            //s.append(((spins[i].spin_bool) ? std::string("ʌ") : std::string("v")) + " "); 
            s.append(((spins[i].spin_bool) ? std::string(" ") : std::string("")) + std::to_string(spins[i]) + " "); 
        }
        return s;
    }
    inline signed char spinAt(int i){
        if(periodic){
            if(i == N_spins) return spins[0];
            else if(i == -1) return spins[N_spins-1];
            else return spins[i];
        }
        //free boundary conditions: return 0 for spins out of bounds
        return (i < N_spins && i >= 0) ? spins[i] : 0;
    }
    // void setBeta(double beta){
    //     flip_probabilities = {exp(-beta*2),exp(-beta*4)};
    // }
    // signed char flip_spin(int i){ //flip a spin, adjust M,E

    //     signed char dE = 2 * (spinAt(i-1)*spinAt(i) + spinAt(i)*spinAt(i+1)); //if i is out of bounds, [] operator returns 0
    //     signed char dM = - 2 * spinAt(i);

    //     spins[i].flip();
    //     E_J += dE;
    //     M += dM;


    //     return spinAt(i);
    // }
    signed char flip_if_criterion(int i){ //flip a spin according to criterion

        signed char dE =   2 * spinAt(i) * (spinAt(i-1) + spinAt(i+1)); //if i is out of bounds, spinAt(i) returns 0
        signed char dM = - 2 * spinAt(i);

        unsigned char dE_index = dE/2 + 2;
        double rand = real_distrib(gen);
        if(rand < flip_probabilities[dE_index]){ //only flip with probability exp(-beta*dE)
            spins[i].flip();
            E_J += dE;
            M += dM;
        }
    

        return spinAt(i);
    }
    void sweep(){
        for(int i = 0; i < N_spins; i++){
            int randomSpin = random_spin_choice_distrib(gen);
            flip_if_criterion(randomSpin);
        }
    }

    SpinSystem(int N_spins, double beta, std::string init_as = "random", std::string update_rule = "metropolis", std::string boundary = "free")
        :   N_spins(N_spins),
            spins(std::vector<Spin>(N_spins)),
            random_spin_choice_distrib(std::uniform_int_distribution<>(0,N_spins-1))
    {
        if(update_rule == "glauber") for(int i = 0; i < 5; i++) flip_probabilities[i] = glauber_prb(beta,(i-2)*2);
        else flip_probabilities = {1,1,1,exp(-beta*2),exp(-beta*4)};

        if(boundary == "periodic") periodic = true;
        else periodic = false;

        if(init_as == "ordered"){
            for(int i = 0; i < N_spins; i++){//initialize ordered spins
                spins[i].spin_bool = true;
            }
            E_J = - N_spins + 1;
            M = N_spins;
        }
        else{
            for(int i = 0; i < N_spins; i++){//initialize random spins, calculate initial M and E
                spins[i].spin_bool = bool_distrib(gen);
                M += spinAt(i);
                if(i > 0) E_J += -(spinAt(i-1) * spinAt(i));
            }
        }   
    }

};

int main(int argc, char* argv[]){
    //SpinSystem system(5,1);

    // std::cout << (std::string)system << std::endl;
    // std::cout << "E = " << system.E_J << ", M = " << system.M << "\n" << std::endl;

    //system.sweep();
    
    std::array<int,4> N_spins_arr = {20,100,1000,10000};
    //std::array<int,1> N_spins_arr = {20};

    std::array<double,21> beta_range;
    for(int i = 0; i < 21; i++) beta_range[i] = (double) i / 10; //define range for beta \in (0, 0.1, ..., 2)


    for(int N_spins : N_spins_arr){
        std::string filename("out/");
        filename.append(std::to_string(N_spins)).append(".txt");
        std::ofstream out(filename.c_str());

        out << "# #beta\t<E_J>\t<M>\tC\tchi\t\tN_spins = " << N_spins << "\t\t#Measurement_sweeps = 10000" << std::endl;
        std::cout << "N_spins = " << N_spins << std::endl;

        for(double beta : beta_range){


            SpinSystem system_random(N_spins,beta,"random","glauber", "periodic");
            SpinSystem system_ordered(N_spins,beta,"ordered","glauber", "periodic");

            //thermalize spin system until random and ordered system cross

            int N_sweeps_until_thermalized = 0;
            std::cout << "\nbeta = " << beta << "\nthermalizing ..." << std::flush;

            //no need for thermalization if beta = 0
            if(beta != 0) for(int i = 0; i < 1000; i++){ //assume thermalization when |E_rand - E_ord| < 1
                system_random.sweep();
                system_ordered.sweep();
                N_sweeps_until_thermalized++;
            }
            std::cout << " done! N_sweeps_until_thermalized = "<< N_sweeps_until_thermalized << std::endl;

            //begin measurement of average
            

            // not used (//rule of thumb: 10% of sweeps used for thermalization, 90% for measurement)
            // int N_sweeps_measurement = 9 * N_sweeps_until_thermalized;

            int N_sweeps_measurement = 10000;

            SpinSystem& system = system_random; //use only one system for measurement

            long sum_E_J = 0;
            long sum_M = 0;

            long sum_E_J_square = 0;
            long sum_M_square = 0;

            std::cout << "measuring: 0 %" << std::flush;
            for(int i = 0; i < N_sweeps_measurement; i++){
                sum_E_J += system.E_J;
                sum_M += abs(system.M);
                sum_E_J_square += system.E_J*system.E_J;
                sum_M_square += system.M*system.M;
                system.sweep();
                std::cout << "\rmeasuring... " << i * 100 / N_sweeps_measurement << " %" << std::flush;
            }
            std::cout << "\rmeasuring... done!" << std::endl;

            double avg_E_J = (double) sum_E_J / N_sweeps_measurement;
            double avg_M = (double) sum_M / N_sweeps_measurement;

            double avg_E_J_squared = (double) sum_E_J_square / N_sweeps_measurement;
            double avg_M_squared = (double) sum_M_square / N_sweeps_measurement;

            double C = beta*beta * (avg_E_J_squared - (avg_E_J*avg_E_J));
            double chi = beta * avg_M_squared;

            out << beta << " " <<  avg_E_J << " " << avg_M << " " << C << " " << chi << "\n"; //print average to file

        }
        out.close();
        std::cout << std::endl;
    }
    
}