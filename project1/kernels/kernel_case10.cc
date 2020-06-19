//Produced at : Fri Jun 19 01:00:37 2020

 #include "../run.h" 
void kernel_case10(float (&B)[10][10], float (&A)[8][8]) {
  for(int i = 0; i < 8; ++i){
    for(int j = 0; j < 8; ++j){
      if (((i+1 >= 0 && i+1 < 10) && (i+2 >= 0 && i+2 < 10))) {
        A[i][j] = (((B[i][j] + B[i+1][j]) + B[i+2][j]) / (float) 3 );
      } else {
      }
    }
  }
}
//Finished at : Fri Jun 19 01:00:37 2020
