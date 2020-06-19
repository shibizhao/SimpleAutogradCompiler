//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case4(float (&B)[16][32], float (&C)[32][32], float (&dA)[16][32], float (&dB)[16][32], float (&dC)[32][32]) {
  for(int i = 0; i < 16; ++i){
    for(int k = 0; k < 32; ++k){
      float tmp_1 = 0;
      for(int j = 0; j < 32; ++j){
        tmp_1 = (tmp_1 + ((float) 0  + (dA[i][j] * C[k][j])));
      }
      dB[i][k] = tmp_1;
    }
  }
  for(int k = 0; k < 32; ++k){
    for(int j = 0; j < 32; ++j){
      float tmp_2 = 0;
      tmp_2 = (tmp_2 + ((float) 0  + dC[k][j]));
      float tmp_3 = 0;
      for(int i = 0; i < 16; ++i){
        tmp_3 = (tmp_3 + ((float) 0  + (B[i][k] * dA[i][j])));
      }
      dC[k][j] = (tmp_2 + tmp_3);
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
