//Produced at : Fri Jun 19 01:00:37 2020

 #include "../run.h" 
void kernel_case5(float (&B)[16][32], float (&C)[32][32], float (&D)[16][32], float & alpha, float & beta, float (&A)[16][32]) {
  for(int i = 0; i < 16; ++i){
    for(int j = 0; j < 32; ++j){
      float tmp_1 = 0;
      for(int k = 0; k < 32; ++k){
        tmp_1 = (tmp_1 + (alpha * (B[i][k] * C[k][j])));
      }
      A[i][j] = (A[i][j] + tmp_1);
    }
  }
  for(int i = 0; i < 16; ++i){
    for(int j = 0; j < 32; ++j){
      A[i][j] = (A[i][j] + (beta * D[i][j]));
    }
  }
}
//Finished at : Fri Jun 19 01:00:37 2020
