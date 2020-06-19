//Produced at : Fri Jun 19 01:00:37 2020

 #include "../run.h" 
void kernel_case3(int (&B)[16][32], int (&C)[16][32], int (&A)[16][32]) {
  for(int i = 0; i < 16; ++i){
    for(int j = 0; j < 32; ++j){
      A[i][j] = (B[i][j] + C[i][j]);
    }
  }
}
//Finished at : Fri Jun 19 01:00:37 2020
