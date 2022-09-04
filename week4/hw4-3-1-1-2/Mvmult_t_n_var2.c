#include <stdio.h>
#include <FLAME.h>
#include "AXPY_Row_to_Column.h"

#define PRINT_DATA

// Declaration of local prototypes.
void matrix_generate( FLA_Obj A );
double spendtimeinsecond();
int Mvmult_t_unb_var2( FLA_Obj A, FLA_Obj x, FLA_Obj y );

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A,x,y;
  double t1,t2;
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
  t1=spendtimeinsecond();
  Mvmult_t_unb_var2(A,x,y);
  t2=spendtimeinsecond();
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " X = [ ", x, "%le", " ];" );
  FLA_Obj_show( " Y = [ ", y, "%le", " ];" );
#endif
  printf("The spent time for computing  A*x in second : %lf\n",t2-t1);
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
int Mvmult_t_unb_var2( FLA_Obj A, FLA_Obj x, FLA_Obj y )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj xT,              x0,
          xB,              chi1,
                           x2;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* *** */
                                              &a1t, 
                           AB,                &A2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                        /* ** */            /* **** */
                                              &chi1, 
                           xB,                &x2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    // y = alpha*x1t + y
    AXPY_Row_to_Column_unb(chi1,a1t,y);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  a1t, 
                            /* ** */           /* *** */
                              &AB,                A2,     FLA_TOP );

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
//
double spendtimeinsecond()
{
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
