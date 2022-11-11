#include <fstream>
#include <string>
#include <math.h>

std::ofstream cpp("build/task5_loops.cpp");
std::ofstream header("build/task5_loops.h");


void genLoops(int i, int N){
    std::string indent = "";
    for(int j = 0; j < i+1 ; ++j){
        indent += "\t";
    }

    if(i < N){
        cpp << indent << "for(int s" << i << " = 0; s" << i << " < 2; ++s" << i << "){" << std::endl;
        genLoops(i+1,N);
        cpp << indent << "}" << std::endl;
    }
    if(i == N){
        cpp << indent << "Z += exp( s0*s1 ";
        for(int j = 1; j < N-1; ++j){
            cpp << "+ s" << j << "*s" << j+1 << " ";
        }
        cpp << ");" << std::endl;
    }
    return;
}

void genZ(int N){
    cpp << "\ndouble Z" << N << "(){\n\tdouble Z = 0;\n" << std::endl;
    header << "\ndouble Z" << N << "();" << std::flush;
    genLoops(0, N);
    cpp << "\n\treturn Z;\n}" << std::flush;
}

int main(int argc, char* argv[]){
    cpp << "#include <math.h>" << std::endl;
    header << "#pragma once" << std::endl;

    genZ(5);
    genZ(10);
    genZ(20);
}