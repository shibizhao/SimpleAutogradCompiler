//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case7(float (&dB)[16][32], float (&dA)[32][16]) {
  for(int j = 0; j < 32; ++j){
    for(int i = 0; i < 16; ++i){
      float tmp_1 = 0;
      tmp_1 = (tmp_1 + ((float) 0  + dB[i][j]));
      dA[j][i] = tmp_1;
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
