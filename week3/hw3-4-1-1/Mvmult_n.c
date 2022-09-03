#include <stdio.h>
#include <FLAME.h>
#define PRINT_DATA

// Declaration of local prototypes.
void matrix_generate( FLA_Obj A );
int Mvmult_n_unb( FLA_Obj A, FLA_Obj x, FLA_Obj y );
void DotProduct(double* at,double* x,int n,int ld,double* alpha);

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
  Mvmult_n_unb(A,x,y);
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
int Mvmult_n_unb( FLA_Obj A, FLA_Obj x, FLA_Obj y )
{
  FLA_Obj AT,              A0,
          AB,              a1t,
                           A2;

  FLA_Obj yT,              y0,
          yB,              psi1,
                           y2;
  double* buff_psi1;
  double* buff_a1t;
  double* buff_x;
  int n_a,ldim_A;
  n_a       = FLA_Obj_width(A);
  ldim_A    = FLA_Obj_col_stride(A);
  buff_x    = (double *) FLA_Obj_buffer_at_view(x);

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );
  
  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* *** */
                                              &a1t, 
                           AB,                &A2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* **** */
                                              &psi1, 
                           yB,                &y2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    buff_a1t  =  (double *) FLA_Obj_buffer_at_view(a1t);
    buff_psi1 =  (double *) FLA_Obj_buffer_at_view(psi1);
    DotProduct(buff_a1t,buff_x,n_a,ldim_A,buff_psi1);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  a1t, 
                            /* ** */           /* *** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                                                  psi1, 
                            /* ** */           /* **** */
                              &yB,                y2,     FLA_TOP );

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
void DotProduct(double* at,double* x,int n,int ld,double* alpha){
  int i;
  for(i=0;i<n;++i){
    *alpha+=at[i*ld]*x[i];
  }
}
