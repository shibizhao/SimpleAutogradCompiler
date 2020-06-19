//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case2(float (&A)[4][16], float (&dB)[4][16], float (&dA)[4][16]) {
  for(int i = 0; i < 4; ++i){
    for(int j = 0; j < 16; ++j){
      float tmp_1 = 0;
      tmp_1 = (tmp_1 + ((float) 0  + (dB[i][j] * A[i][j])));
      float tmp_2 = 0;
      tmp_2 = (tmp_2 + ((float) 0  + (A[i][j] * dB[i][j])));
      float tmp_3 = 0;
      tmp_3 = (tmp_3 + (float) 0 );
      dA[i][j] = ((tmp_1 + tmp_2) + tmp_3);
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
