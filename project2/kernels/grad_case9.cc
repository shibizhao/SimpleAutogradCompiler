//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case9(float (&dB)[4][6], float (&dA)[4]) {
  for(int i = 0; i < 4; ++i){
    float tmp_1 = 0;
    for(int j = 0; j < 6; ++j){
      tmp_1 = (tmp_1 + ((float) 0  + dB[i][j]));
    }
    dA[i] = tmp_1;
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
