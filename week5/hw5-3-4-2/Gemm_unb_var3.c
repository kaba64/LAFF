#include <stdio.h>
# include <sys/time.h>
#include <FLAME.h>
#include "ger.h"

#define PRINT_DATA
// Declaration of local prototypes.
void matrix_generate( FLA_Obj A, char type);
double spendtimeinsecond();
int Gemm_unb_var3( FLA_Obj A, FLA_Obj B, FLA_Obj C );

int main( int argc, char *argv[] ) {
  int m_A,n_A,m_B,n_B;
  FLA_Obj A,B,C;
  double t1,t2;
  printf("Please enter the size of the matrix A : m_A and n_A \n");
  scanf("%d%d",&m_A,&n_A);
  printf("Please enter the size of the matrix B : m_B and n_B (m_B should be equal to n_A) \n");
  scanf("%d%d",&m_B,&n_B);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  FLA_Obj_create(FLA_DOUBLE,m_B,n_B,0,0,&B);
  FLA_Obj_create(FLA_DOUBLE,m_A,n_B,0,0,&C);
  matrix_generate(A,'F');
  matrix_generate(B,'F');
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " A = [ ", B, "%le", " ];" );
  FLA_Obj_show( " A = [ ", C, "%le", " ];" );
#endif
  t1=spendtimeinsecond();
  Gemm_unb_var3(A,B,C);
  t2=spendtimeinsecond();
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " A = [ ", B, "%le", " ];" );
  FLA_Obj_show( " A = [ ", C, "%le", " ];" );
#endif
  printf("The spent time for computing  A*x in second : %lf\n",t2-t1);
  FLA_Obj_free(&A);
  FLA_Obj_free(&B);
  FLA_Obj_free(&C);
  FLA_Finalize();
  return 0;
}
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int Gemm_unb_var3( FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj BT,              B0,
          BB,              b1t,
                           B2;

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
    ger_unb(a1,b1t,C);
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
//
void matrix_generate( FLA_Obj A, char type) {
  double  * buff_A;
  int  m_A, n_A, ldim_A;
  int  num;
  
  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );
  num = 1;
  if(type=='F'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = 0; i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
    }
  }else if(type=='S'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = 0; i < m_A; i++ ) {
	if(i>=j){
	  buff_A[ i + j * ldim_A ] = ( double ) num;
	  num++;
	}else{
	  buff_A[ i + j * ldim_A ]=buff_A[j+i*ldim_A];
	}
      }
    }
  }else if(type=='U'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = 0; i < j+1; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
    }
  }else if(type=='L'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = j; i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
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
