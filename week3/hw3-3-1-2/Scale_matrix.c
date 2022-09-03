#include <stdio.h>
#include <FLAME.h>
#define PRINT_DATA

// Declaration of local prototypes.
void matrix_generate( FLA_Obj A );
void ScaleColumnMatrix(double* x,int n,double alpha);
int Scale_matrix_unb( FLA_Obj alpha, FLA_Obj A );

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A,alpha;
  double* buff_alpha;
  printf("Please enter the size of the matrix m_A and n_A: \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&alpha);
  buff_alpha = (double *) FLA_Obj_buffer_at_view(alpha);
  printf("Please enter the scalar : \n");
  scanf("%lf",buff_alpha);
  matrix_generate(A);
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", alpha, "%le", " ];" );
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
#endif
  Scale_matrix_unb(alpha,A);
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
#endif
  printf("End of Program\n");
  FLA_Obj_free(&A);
  FLA_Obj_free(&alpha);
  FLA_Finalize();
  return 0;
}
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int Scale_matrix_unb( FLA_Obj alpha, FLA_Obj A )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;
  double value_alpha;
  double* buff_a1;
  int m_a1;
  
  value_alpha = *((double *) FLA_Obj_buffer_at_view(alpha));
  m_a1       =            FLA_Obj_length(A);
  
  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    buff_a1  =  (double *) FLA_Obj_buffer_at_view(a1);
    ScaleColumnMatrix(buff_a1,m_a1,value_alpha);
    /*------------------------------------------------------------*/
    
    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );
    
  }
  
  return FLA_SUCCESS;
}
//
void matrix_generate( FLA_Obj A ) {
  double  * buff_A;
  int     m_A, n_A, ldim_A;
  int     i, j, num;

  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );
  num = 1;
  for ( j = 0; j < n_A; j++ ) {
    for ( i = 0; i < m_A; i++ ) {
      buff_A[ i + j * ldim_A ] = ( double ) num;
      num++;
    }
  }
}
void ScaleColumnMatrix(double* x,int m,double alpha){
  int i;
  for(i=0;i<m;++i){
    x[i]=alpha*x[i];
  }
}
