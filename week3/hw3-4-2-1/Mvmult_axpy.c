#include <stdio.h>
#include <FLAME.h>
#include "axpy.h"

#define PRINT_DATA

// Declaration of local prototypes.
void matrix_generate( FLA_Obj A );
int Mvmult_axpy_unb(FLA_Obj A,FLA_Obj x,FLA_Obj y);

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A,x,y;
  printf("Please enter the size of the matrix m_A and n_A: \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  FLA_Obj_create(FLA_DOUBLE,n_A,1,0,0,&x);
  FLA_Obj_create(FLA_DOUBLE,m_A,1,0,0,&y);
  matrix_generate(A);
  matrix_generate(x);
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " X = [ ", x, "%le", " ];" );
  FLA_Obj_show( " Y = [ ", y, "%le", " ];" );
#endif
  Mvmult_axpy_unb(A,x,y);
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " X = [ ", x, "%le", " ];" );
  FLA_Obj_show( " Y = [ ", y, "%le", " ];" );
#endif
  printf("End of Program\n");
  FLA_Obj_free(&A);
  FLA_Obj_free(&x);
  FLA_Obj_free(&y);
  FLA_Finalize();
  return 0;
}
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int Mvmult_axpy_unb( FLA_Obj A, FLA_Obj x, FLA_Obj y )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj xT,              x0,
          xB,              chi1,
                           x2;
  
  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_TOP );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                        /* ** */            /* **** */
                                              &chi1, 
                           xB,                &x2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    // y = y + chi1*a1
    AXPY_unb(chi1,a1,y);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &xT,                x0, 
                                                  chi1, 
                            /* ** */           /* **** */
                              &xB,                x2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}
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
