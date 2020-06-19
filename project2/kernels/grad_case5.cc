//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case5(float (&C)[32][32], float (&D)[4][32], float (&dA)[16][32], float (&dB)[16][32][4]) {
  for(int i = 0; i < 16; ++i){
    for(int k = 0; k < 32; ++k){
      for(int l = 0; l < 4; ++l){
        float tmp_1 = 0;
        for(int j = 0; j < 32; ++j){
          tmp_1 = (tmp_1 + ((float) 0  + ((dA[i][j] * C[k][j]) * D[l][j])));
        }
        dB[i][k][l] = tmp_1;
      }
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
