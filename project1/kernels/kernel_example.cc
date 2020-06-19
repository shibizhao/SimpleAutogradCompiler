//Produced at : Fri Jun 19 01:00:37 2020

 #include "../run.h" 
void kernel_example(float (&B)[32][16], float (&C)[32][16], float (&A)[32][16]) {
  for(int i = 0; i < 32; ++i){
    for(int j = 0; j < 16; ++j){
      A[i][j] = (C[i][j] * B[i][j]);
    }
  }
}
//Finished at : Fri Jun 19 01:00:37 2020
