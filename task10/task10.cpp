#include <fstream>
#include <math.h>

long long prev_number = 1LL;
long a = 25214903917;
int c = 11;
long long m = 1LL << 48;

long double LCG(){
    prev_number = fabs((a*prev_number+c) % m);
    return (double) (prev_number) / m;
}
int main(int argc, char* argv[]){
    std::ofstream out("out/lcg.txt");
    for(int i = 0; i < 1e5; ++i){
        out << LCG() << " ";
    }
    out << std::flush;
}