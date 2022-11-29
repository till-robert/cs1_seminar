#include <iostream>
#include <random>

std::mt19937_64 rng(54321);
std::uniform_real_distribution<double>rnd(-1, 1);
int N = 1e4;

int main(int argc, char* arcv[]){
    int hits = 0;
    for(int i = 0; i < N; ++i){
        double x = rnd(rng);
        double y = rnd(rng);
        if(x*x + y*y < 1){
            hits++;
        }
    }

    std::cout << "hits / trials = " << (double)hits / N << std::endl;
}