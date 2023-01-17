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
const int seed = 41;
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
            E_J = - N_spins + !periodic; //if free b. c., there is 1 less interaction
            M = N_spins;
        }
        else{
            //add 1 if periodic to account for additional interaction
            for(int i = 0; i < N_spins + periodic; i++){//initialize random spins, calculate initial M and E
                spins[i].spin_bool = bool_distrib(gen);
                M += spinAt(i);
                if(i > 0) E_J += -(spinAt(i-1) * spinAt(i));
            }
        }   
    }

};

int main(int argc, char* argv[]){

    
    std::ofstream out("out/timeseries.txt");
    int L = 100;
    double beta = 1.0;

    SpinSystem system(L,beta,"random","metropolis","free");
        for(int sweep = 0; sweep < 10000; sweep++){
            system.sweep();
            out << system.E_J << " " << system.M << "\n";

        }
    out.close();
}