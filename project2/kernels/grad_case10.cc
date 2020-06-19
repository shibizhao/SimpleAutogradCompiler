//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case10(float (&dA)[8][8], float (&dB)[10][8]) {
  for(int t_2 = 0; t_2 < 10; ++t_2){
    for(int j = 0; j < 8; ++j){
      if ((j >= 0 && j < 8)) {
        float tmp_1 = 0;
        if (((j >= 0 && j < 8) && ((t_2-2) >= 0 && (t_2-2) < 8))) {
          tmp_1 = (tmp_1 + ((float) 0  + dA[(t_2-2)][j]));
        } else {
        }
        float tmp_2 = 0;
        if (((j >= 0 && j < 8) && ((t_2-2)+1 >= 0 && (t_2-2)+1 < 8))) {
          tmp_2 = (tmp_2 + ((float) 0  + dA[(t_2-2)+1][j]));
        } else {
        }
        float tmp_3 = 0;
        if (((j >= 0 && j < 8) && (t_2 >= 0 && t_2 < 8))) {
          tmp_3 = (tmp_3 + ((float) 0  + dA[t_2][j]));
        } else {
        }
        float tmp_4 = 0;
        tmp_4 = (tmp_4 + ((float) 0  + (((tmp_1 + tmp_2) + tmp_3) / (float) 3 )));
        dB[t_2][j] = tmp_4;
      } else {
      }
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
