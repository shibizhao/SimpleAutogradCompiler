//Produced at : Fri Jun 19 01:00:39 2020
#include "../run2.h" 
void grad_case8(float (&dB)[32], float (&dA)[2][16]) {
  for(int z_1 = 0; z_1 < 2; ++z_1){
    for(int z_2 = 0; z_2 < 16; ++z_2){
      float tmp_1 = 0;
      if (((z_1*16+z_2) >= 0 && (z_1*16+z_2) < 32)) {
        tmp_1 = (tmp_1 + ((float) 0  + dB[(z_1*16+z_2)]));
      } else {
      }
      dA[z_1][z_2] = tmp_1;
    }
  }
}
//Finished at : Fri Jun 19 01:00:39 2020
