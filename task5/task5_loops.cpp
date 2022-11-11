#include <math.h>

double Z5(){
	double Z = 0;

	for(int s0 = 0; s0 < 2; ++s0){
		for(int s1 = 0; s1 < 2; ++s1){
			for(int s2 = 0; s2 < 2; ++s2){
				for(int s3 = 0; s3 < 2; ++s3){
					for(int s4 = 0; s4 < 2; ++s4){
						Z += exp( s0*s1 + s1*s2 + s2*s3 + s3*s4 );
					}
				}
			}
		}
	}

	return Z;
}
double Z10(){
	double Z = 0;

	for(int s0 = 0; s0 < 2; ++s0){
		for(int s1 = 0; s1 < 2; ++s1){
			for(int s2 = 0; s2 < 2; ++s2){
				for(int s3 = 0; s3 < 2; ++s3){
					for(int s4 = 0; s4 < 2; ++s4){
						for(int s5 = 0; s5 < 2; ++s5){
							for(int s6 = 0; s6 < 2; ++s6){
								for(int s7 = 0; s7 < 2; ++s7){
									for(int s8 = 0; s8 < 2; ++s8){
										for(int s9 = 0; s9 < 2; ++s9){
											Z += exp( s0*s1 + s1*s2 + s2*s3 + s3*s4 + s4*s5 + s5*s6 + s6*s7 + s7*s8 + s8*s9 );
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

	return Z;
}
double Z20(){
	double Z = 0;

	for(int s0 = 0; s0 < 2; ++s0){
		for(int s1 = 0; s1 < 2; ++s1){
			for(int s2 = 0; s2 < 2; ++s2){
				for(int s3 = 0; s3 < 2; ++s3){
					for(int s4 = 0; s4 < 2; ++s4){
						for(int s5 = 0; s5 < 2; ++s5){
							for(int s6 = 0; s6 < 2; ++s6){
								for(int s7 = 0; s7 < 2; ++s7){
									for(int s8 = 0; s8 < 2; ++s8){
										for(int s9 = 0; s9 < 2; ++s9){
											for(int s10 = 0; s10 < 2; ++s10){
												for(int s11 = 0; s11 < 2; ++s11){
													for(int s12 = 0; s12 < 2; ++s12){
														for(int s13 = 0; s13 < 2; ++s13){
															for(int s14 = 0; s14 < 2; ++s14){
																for(int s15 = 0; s15 < 2; ++s15){
																	for(int s16 = 0; s16 < 2; ++s16){
																		for(int s17 = 0; s17 < 2; ++s17){
																			for(int s18 = 0; s18 < 2; ++s18){
																				for(int s19 = 0; s19 < 2; ++s19){
																					Z += exp( s0*s1 + s1*s2 + s2*s3 + s3*s4 + s4*s5 + s5*s6 + s6*s7 + s7*s8 + s8*s9 + s9*s10 + s10*s11 + s11*s12 + s12*s13 + s13*s14 + s14*s15 + s15*s16 + s16*s17 + s17*s18 + s18*s19 );
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

	return Z;
}