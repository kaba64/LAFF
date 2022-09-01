#include <stdio.h>
#include <FLAME.h>
#define PRINT_DATA

// Declaration of local prototypes.
void matrix_generate( FLA_Obj A );
void CopyTranspose(double* x,double* y,int m,int ldim);
int Transpose_unb( FLA_Obj A, FLA_Obj B );

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A,AT;
  
  printf("Please enter the size of the matrix m_A and n_A: \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  FLA_Obj_create(FLA_DOUBLE,n_A,m_A,0,0,&AT);
  matrix_generate(A);
  Transpose_unb(A,AT);
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " AT = [ ", AT, "%le", " ];" );
#endif
  printf("End of Program\n");
  FLA_Obj_free(&A);
  FLA_Obj_free(&AT);
  FLA_Finalize();
  return 0;
}
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Transpose_unb( FLA_Obj A, FLA_Obj B )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj BT,              B0,
          BB,              b1t,
                           B2;
  double* buff_a1;
  double* buff_b1t;
  int m_a1,ldim_B;

  ldim_B = FLA_Obj_col_stride(B);
  
  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_2x1( B,    &BT, 
                      &BB,            0, FLA_TOP );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( BT,                &B0, 
                        /* ** */            /* *** */
                                              &b1t, 
                           BB,                &B2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    buff_a1   =  (double *) FLA_Obj_buffer_at_view(a1);
    m_a1  = FLA_Obj_length(a1);
    buff_b1t   =  (double *) FLA_Obj_buffer_at_view(b1t);
    CopyTranspose(buff_a1,buff_b1t,m_a1,ldim_B);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0, 
                                                  b1t, 
                            /* ** */           /* *** */
                              &BB,                B2,     FLA_TOP );

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

void CopyTranspose(double* x,double* y,int m,int ldim){
  int i;
  for(i=0;i<m;++i){
    y[i*ldim]=x[i];
  }
}
