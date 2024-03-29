#include <stdlib.h>
#include <stdio.h>
#include <FLAME.h>

#define PRINT_DATA

// Declaration of local prototypes.
int ZeroMatrix_D_unb(FLA_Obj A);
void SetZero_D(double* x,int m, int n);

int main( int argc, char *argv[] ) {
  int m_A, n_A;
  FLA_Obj  A;

  printf("Please enter the size of the matrix m_A and n_A: \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A); //FLA_INT
  ZeroMatrix_D_unb(A);
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
#endif
  printf("End of Program\n" );
  FLA_Finalize();
  return 0;
}

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int ZeroMatrix_D_unb(FLA_Obj A)
{
  FLA_Obj AL,AR,A0,a1,A2;
  double  * buff_a1;
  int m_A, n_A, ldim_A;
  int j,count;

  count = 0;  
  FLA_Part_1x2( A,&AL,&AR,0,FLA_LEFT);

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    buff_a1 = ( double * ) FLA_Obj_buffer_at_view( a1 );
    m_A    = FLA_Obj_length(a1);
    n_A    = FLA_Obj_width (a1);
    SetZero_D(buff_a1,m_A,n_A);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}

void SetZero_D(double * x,int m, int n){
  int i,j;
  int count;
  count=0;
  if(n==1){
    for(j=0;j<m;++j){
      x[j]=(double)count;
      ++count;
    }
  }else{
    for(i=0;i<n;++i){
      x[i]=(double)count;
      ++count;
    }
  }
}
