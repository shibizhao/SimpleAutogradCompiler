//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case1(float (&B)[4][16], float (&dC)[4][16], float (&dA)[4][16]) {
  for(int i = 0; i < 4; ++i){
    for(int j = 0; j < 16; ++j){
      float tmp_1 = 0;
      tmp_1 = (tmp_1 + ((float) 0  + (dC[i][j] * B[i][j])));
      float tmp_2 = 0;
      tmp_2 = (tmp_2 + (float) 0 );
      dA[i][j] = (tmp_1 + tmp_2);
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
