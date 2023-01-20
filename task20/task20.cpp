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
std::uniform_int_distribution<> neighbor_distrib(0,3);
struct RandomNeighbor //used for 2D positions
{
    int x, y;

    RandomNeighbor() = default;
    RandomNeighbor(int p_x, int p_y) { //direction can be 0: LEFT, 1: UP, 2: RIGHT, 3: DOWN
            int dir = neighbor_distrib(gen);
            switch (dir)
            {
            case 0:
                x = p_x-1;
                y = p_y;
                break;
            case 1:
                x = p_x;
                y = p_y+1;
                break;
            case 2:
                x = p_x+1;
                y = p_y;
                break;
            case 3:
                x = p_x;
                y = p_y-1;
                break;
            }
    }
};

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
    inline Spin flip(){
        spin_bool = !spin_bool;
        return *this;
    }
};
class SpinSystem{
public:
    int L;
    bool periodic = false;
    bool kawasaki = false;

    std::vector<std::vector<Spin>> spins;

    std::vector<double> flip_probabilities;
    std::uniform_int_distribution<>  random_spin_choice_distrib;

    int E_J = 0; //energy without magnetic field contribution: E_tot = E_J - h*M
    int M = 0;
    int acceptance_per_sweep = 0;
    
    inline static double glauber_prb(double beta,signed char E){
        return 1/(1+exp(beta*E));
    }

    
    operator std::string (){
        std::string s = "";
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){

            //s.append(((spins[i].spin_bool) ? std::string("ʌ") : std::string("v")) + " "); 
            s.append(((spins[i][j].spin_bool) ? std::string(" ") : std::string("")) + std::to_string(spins[i][j]) + " "); 
            }
            s.append("\n");
        }
        return s;
    }
    inline signed char spinAt(int i, int j){
        if(periodic){
            if(i <= -1) i += L;
            else if(i >= L) i -= L;

            if(j <= -1) j += L;
            else if(j >= L) j -= L;

            return spins[i][j];
        }
        //free boundary conditions: return 0 for spins out of bounds
        return (i < L && i >= 0 && j < L && j >= 0) ? spins[i][j] : 0;
    }
    inline signed char safe_flip(int i, int j){
        if(periodic){
            if(i <= -1) i += L;
            else if(i >= L) i -= L;

            if(j <= -1) j += L;
            else if(j >= L) j -= L;

            return spins[i][j].flip();
        }
        //free boundary conditions: return 0 for spins out of bounds
        return (i < L && i >= 0 && j < L && j >= 0) ? spins[i][j].flip() : 0;
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
    signed char flip_if_criterion(int i, int j){ //flip a spin according to criterion
        if(!kawasaki){
            signed char dE =   2 * spinAt(i,j) * (spinAt(i-1,j) + spinAt(i+1,j)+ spinAt(i,j-1) + spinAt(i,j+1)); //if i is out of bounds, spinAt(i) returns 0
            signed char dM = - 2 * spinAt(i,j);

            unsigned char dE_index = dE/2 + 4;
            double rand = real_distrib(gen);
            if(rand < flip_probabilities[dE_index]){ //only flip with probability exp(-beta*dE)
                spins[i][j].flip();
                acceptance_per_sweep ++;
                E_J += dE;
                M += dM;
            }
            return spinAt(i,j);
        }
        else{//kawasaki => M=const
            RandomNeighbor rn(i,j);
            signed char dE =   2 * spinAt(i,j) * (spinAt(i-1,j) + spinAt(i+1,j)+ spinAt(i,j-1) + spinAt(i,j+1)); //if i is out of bounds, spinAt(i) returns 0
                             + 2 * spinAt(rn.x,rn.y) * (spinAt(rn.x-1,rn.y) + spinAt(rn.x+1,rn.y)+ spinAt(rn.x,rn.y-1) + spinAt(rn.x,rn.y+1)); //if i is out of bounds, spinAt(i) returns 0
            signed char dM = - 2 * (spinAt(i,j) + spinAt(rn.x,rn.y));

            unsigned char dE_index = dE/2 + 8;
            double rand = real_distrib(gen);
            if(rand < flip_probabilities[dE_index] && dM == 0){ //only flip with probability exp(-beta*dE) and the spins are opposite
                spins[i][j].flip();
                safe_flip(rn.x,rn.y); //this method checks for out of bounds first
                acceptance_per_sweep ++;
                E_J += dE;
            }
            return spinAt(i,j);

        }
    }
    void sweep(){
        acceptance_per_sweep = 0;
        for(int i = 0; i < L*L; i++){
            int random_i = random_spin_choice_distrib(gen);
            int random_j = random_spin_choice_distrib(gen);
            flip_if_criterion(random_i,random_j);
        }
    }

    SpinSystem(int L, double beta, std::string init_as = "random", std::string update_rule = "metropolis", std::string boundary = "free")
        :   L(L),
            random_spin_choice_distrib(std::uniform_int_distribution<>(0,L-1))
    {
        std::vector<Spin> row(L);

        if(update_rule == "glauber") for(int i = 0; i < 9; i++) flip_probabilities[i] = glauber_prb(beta,(i-4)*2);
        else if(update_rule == "kawasaki"){
            flip_probabilities = {1,1,1,1,1,1,1,1,1,exp(-beta*2),exp(-beta*4),exp(-beta*6),exp(-beta*8),exp(-beta*10),exp(-beta*12),exp(-beta*14),exp(-beta*16)};
            kawasaki = true;
        }
        else flip_probabilities = {1,1,1,1,1,exp(-beta*2),exp(-beta*4),exp(-beta*6),exp(-beta*8)}; //"normal" metropolis

        if(boundary == "periodic") periodic = true;
        else periodic = false;

        if(init_as == "ordered"){

            for(int i = 0; i < L; i++){//initialize ordered spins
                spins.push_back(row);
                spins[i].reserve(L);
                for(int j = 0; j < L; j++){
                    spins[i][j].spin_bool = true;
                }
            }
            E_J = - L*L + 2*L*!periodic; // if free b. c. there are less interactions
            
            M = L*L;
        }
        else{
            //add 1 if periodic to account for additional interaction
            for(int i = 0; i < L+periodic; i++){//initialize random spins, calculate initial M and E
                spins.push_back(row);
                spins[i].reserve(L);
                for(int j = 0; j < L+periodic; j++){
                    spins[i][j].spin_bool = bool_distrib(gen);
                    M += spinAt(i,j);
                    if(i > 0 && j > 0) E_J += -spinAt(i,j)*(spinAt(i-1,j) + spinAt(i,j-1));
                }
            }
        }   
    }

};

int main(int argc, char* argv[]){
    int N_sweeps_measurement = 10000;
    double beta = 0.881;
    int L = 32;


    std::ofstream out("out/avg_timeseries.txt");
    std::ofstream snapshots("out/snapshots.txt");
    
    int N_systems = 5;
    std::vector<SpinSystem> systems;
    systems.reserve(N_systems);
    for(int j = 0; j < N_systems; j++){
        systems.emplace_back(L,beta,"random","kawasaki", "periodic");
    }

    //begin measurement



    std::cout << "measuring: 0 %" << std::flush;
    for(int i = 0; i < N_sweeps_measurement; i++){
        int sum_E = 0;
        int sum_acceptance = 0;
        if(i == 9 || i == 99 || i == 999 || i == 9999){
            snapshots << (std::string) systems[0] << "\n";
        }
        for(int j = 0; j < N_systems; j++){
            sum_E += systems[j].E_J;
            sum_acceptance += systems[j].acceptance_per_sweep;
            systems[j].sweep();
        }
        out << (double) sum_E / N_systems << " " <<  (double) sum_acceptance / N_systems << "\n"; //print average to file

        std::cout << "\rmeasuring... " << i * 100 / N_sweeps_measurement << " %" << std::flush;
    }
    std::cout << "\rmeasuring... done!" << std::endl;
    
    out.close();
    snapshots.close();
}