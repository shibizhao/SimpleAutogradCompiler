//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case3(float (&B)[16][16], float (&dC)[4][16], float (&dA)[4][16]) {
  for(int i = 0; i < 4; ++i){
    for(int k = 0; k < 16; ++k){
      float tmp_1 = 0;
      for(int j = 0; j < 16; ++j){
        tmp_1 = (tmp_1 + ((float) 0  + (dC[i][j] * B[k][j])));
      }
      dA[i][k] = tmp_1;
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
