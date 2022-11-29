#include <math.h>
#include <vector>
#include <fstream>


int i2E(int i, int N){
    return (i+0.5)*2 - N;
} // 0 => -4, 1 => -2, 2 => 0, 3 => 2, 4 => 4
int E2i(int E, int N){
    return (E+N)/2 + 0.5;
}
int i2M(int i, int N){
    return i*2 - N;
} // 0 => -5, 1 => -3, 2 => -1, 3 => 1, 4 => 3, 5 => 5
int M2i(int M, int N){
    return (M+N)/2;
}
std::vector<std::vector<int>> omega5(){
    int E = 0;
	int M = 0;

	//Number of columns
	int num_M = 6;

	// Number of rows
	int num_E = 5;

	// Initializing a single row
	std::vector<int> row(num_E, 0);

    std::vector<std::vector<int>> omega (num_M,row);

	for(int s0 = -1; s0 < 2; s0 += 2){
		for(int s1 = -1; s1 < 2; s1 += 2){
			for(int s2 = -1; s2 < 2; s2 += 2){
				for(int s3 = -1; s3 < 2; s3 += 2){
					for(int s4 = -1; s4 < 2; s4 += 2){
                        E = s0*s1 + s1*s2 + s2*s3 + s3*s4;
						M = s0+s1+s2+s3+s4;
                        omega[M2i(M,5)][E2i(E,5)]++;
					}
				}
			}
		}
	}
	return omega;
}
std::vector<std::vector<int>> omega20(){
    int E = 0;
	int M = 0;

	//Number of columns
	int num_M = 21;

	// Number of rows
	int num_E = 20;

	// Initializing a single row
	std::vector<int> row(num_E, 0);

    std::vector<std::vector<int>> omega (num_M,row);

	for(int s0 = -1; s0 < 2; s0 += 2){
		for(int s1 = -1; s1 < 2; s1 += 2){
			for(int s2 = -1; s2 < 2; s2 += 2){
				for(int s3 = -1; s3 < 2; s3 += 2){
					for(int s4 = -1; s4 < 2; s4 += 2){
						for(int s5 = -1; s5 < 2; s5 += 2){
							for(int s6 = -1; s6 < 2; s6 += 2){
								for(int s7 = -1; s7 < 2; s7 += 2){
									for(int s8 = -1; s8 < 2; s8 += 2){
										for(int s9 = -1; s9 < 2; s9 += 2){
											for(int s10 = -1; s10 < 2; s10 += 2){
												for(int s11 = -1; s11 < 2; s11 += 2){
													for(int s12 = -1; s12 < 2; s12 += 2){
														for(int s13 = -1; s13 < 2; s13 += 2){
															for(int s14 = -1; s14 < 2; s14 += 2){
																for(int s15 = -1; s15 < 2; s15 += 2){
																	for(int s16 = -1; s16 < 2; s16 += 2){
																		for(int s17 = -1; s17 < 2; s17 += 2){
																			for(int s18 = -1; s18 < 2; s18 += 2){
																				for(int s19 = -1; s19 < 2; s19 += 2){
																					E = s0*s1 + s1*s2 + s2*s3 + s3*s4 + s4*s5 + s5*s6 + s6*s7 + s7*s8 + s8*s9 + s9*s10 + s10*s11 + s11*s12 + s12*s13 + s13*s14 + s14*s15 + s15*s16 + s16*s17 + s17*s18 + s18*s19;
																					M = s0+s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19;
																					omega[M2i(M,20)][E2i(E,20)]++;
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return omega;
}

int main(int argc, char* argv[]){

	std::ofstream out("out/omega.txt");
	auto omega = omega20();

	const int height = omega.size();
	const int width = omega[0].size();
	//std::cout << width << " " << height << std::endl;
	for(int col = 0; col < width; ++col){
		for(int row = 0; row < height; ++row){
			out << omega[row][col] << " ";
		}
		out << "\n";
	}
	out << std::flush;
}