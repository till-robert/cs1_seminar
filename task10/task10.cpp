#include <fstream>
#include <string>
#include <math.h>
long prev_number = 1;
long a = 25214903917;
int c = 11;
unsigned long long  m = pow(2.0,48.0);

long double LCG(){
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