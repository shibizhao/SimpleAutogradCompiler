//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case6(float (&C)[8][16][3][3], float (&dA)[2][8][5][5], float (&dB)[2][16][7][7]) {
  for(int n = 0; n < 2; ++n){
    for(int c = 0; c < 16; ++c){
      for(int z_1 = 0; z_1 < 7; ++z_1){
        for(int z_2 = 0; z_2 < 7; ++z_2){
          if (((n >= 0 && n < 2) && (c >= 0 && c < 16))) {
            float tmp_1 = 0;
            for(int k = 0; k < 8; ++k){
              for(int r = 0; r < 3; ++r){
                for(int s = 0; s < 3; ++s){
                  if (((s >= 0 && s < 3) && ((r >= 0 && r < 3) && ((c >= 0 && c < 16) && (((z_2-s) >= 0 && (z_2-s) < 5) && (((z_1-r) >= 0 && (z_1-r) < 5) && ((k >= 0 && k < 8) && (n >= 0 && n < 2)))))))) {
                    tmp_1 = (tmp_1 + ((float) 0  + (dA[n][k][(z_1-r)][(z_2-s)] * C[k][c][r][s])));
                  } else {
                  }
                }
              }
            }
            dB[n][c][z_1][z_2] = tmp_1;
          } else {
          }
        }
      }
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
