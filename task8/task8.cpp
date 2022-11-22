#include <fstream>

int prev_number = 1;
double LCG(int a = 5, int c = 0, int m = 2048){
    prev_number = ((a*prev_number+c) % m);
    return (double) prev_number / m;
}
int main(int argc, char* argv[]){
    std::ofstream out("out/lcg.txt");
    for(int i = 0; i < 1e5; ++i){
        out << LCG() << " ";
    }
    out << std::flush;
}