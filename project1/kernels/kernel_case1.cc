//Produced at : Fri Jun 19 01:00:37 2020

 #include "../run.h" 
void kernel_case1(float (&A)[32][16]) {
  for(int i = 0; i < 32; ++i){
    for(int j = 0; j < 16; ++j){
      A[i][j] = (float) 2 ;
    }
  }
}
//Finished at : Fri Jun 19 01:00:37 2020
