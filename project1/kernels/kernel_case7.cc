//Produced at : Fri Jun 19 01:00:37 2020

 #include "../run.h" 
void kernel_case7(float (&A)[32][16], float (&B)[16][32]) {
  for(int i = 0; i < 16; ++i){
    for(int j = 0; j < 32; ++j){
      B[i][j] = A[j][i];
    }
  }
}
//Finished at : Fri Jun 19 01:00:37 2020
