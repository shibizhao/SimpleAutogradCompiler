//Produced at : Fri Jun 19 01:00:37 2020

 #include "../run.h" 
void kernel_case9(float (&B)[16][32][8], float (&C)[32][32], float (&D)[8][32], float (&A)[16][32]) {
  for(int i = 0; i < 16; ++i){
    for(int j = 0; j < 32; ++j){
      float tmp_1 = 0;
      for(int k = 0; k < 32; ++k){
        for(int l = 0; l < 8; ++l){
          tmp_1 = (tmp_1 + ((B[i][k][l] * C[k][j]) * D[l][j]));
        }
      }
      A[i][j] = (A[i][j] + tmp_1);
    }
  }
}
//Finished at : Fri Jun 19 01:00:37 2020
