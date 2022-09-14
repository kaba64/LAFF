#include <stdio.h>
# include <sys/time.h>
#include <FLAME.h>
#include "axpy-gaussian.h"
#include "axpy-gaussian-rank1.h"
#include "ger.h"
#include "dot-gaussian.h"

#define PRINT_DATA
// Declaration of local prototypes.
void matrix_generate( FLA_Obj A, char type,char name);
double spendtimeinsecond();
int LU_unb_var5( FLA_Obj A );
int LTRSV_unb_var1( FLA_Obj L, FLA_Obj b );
int Utrsv_unb_var1( FLA_Obj U, FLA_Obj b );

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A,b;
  double t1,t2;
  printf("Please enter the size of the matrix A : m_A and n_A \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  FLA_Obj_create(FLA_DOUBLE,m_A,1,0,0,&b);
  matrix_generate(A,'F','A');
  matrix_generate(b,'F','b');
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " b = [ ", b, "%le", " ];" );
#endif
  t1=spendtimeinsecond();
  LU_unb_var5(A);
  LTRSV_unb_var1(A,b);
  Utrsv_unb_var1(A,b);
  t2=spendtimeinsecond();
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " x = [ ", b, "%le", " ];" );
#endif
  printf("The spent time for solving A*x=b in second : %lf\n",t2-t1);
  FLA_Obj_free(&A);
  FLA_Obj_free(&b);
  FLA_Finalize();
  return 0;
}
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int LU_unb_var5( FLA_Obj A )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/
    AXPY_Gaussian_unb(alpha11,a21);
    ger_unb(a21,a12t,A22);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}
//
int LTRSV_unb_var1( FLA_Obj L, FLA_Obj b )
{
  FLA_Obj LTL,   LTR,      L00,  l01,      L02, 
          LBL,   LBR,      l10t, lambda11, l12t,
                           L20,  l21,      L22;

  FLA_Obj bT,              b0,
          bB,              beta1,
                           b2;

  FLA_Part_2x2( L,    &LTL, &LTR,
                      &LBL, &LBR,     0, 0, FLA_TL );

  FLA_Part_2x1( b,    &bT, 
                      &bB,            0, FLA_TOP );

  while ( FLA_Obj_length( LTL ) < FLA_Obj_length( L ) ){

    FLA_Repart_2x2_to_3x3( LTL, /**/ LTR,       &L00,  /**/ &l01,      &L02,
                        /* ************* */   /* *************************** */
                                                &l10t, /**/ &lambda11, &l12t,
                           LBL, /**/ LBR,       &L20,  /**/ &l21,      &L22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( bT,                &b0, 
                        /* ** */            /* ***** */
                                              &beta1, 
                           bB,                &b2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    AXPY_Gaussian_Rank1_unb(beta1,l21,b2);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &LTL, /**/ &LTR,       L00,  l01,      /**/ L02,
                                                     l10t, lambda11, /**/ l12t,
                            /* ************** */  /* ************************* */
                              &LBL, /**/ &LBR,       L20,  l21,      /**/ L22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &bT,                b0, 
                                                  beta1, 
                            /* ** */           /* ***** */
                              &bB,                b2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}
//
int Utrsv_unb_var1( FLA_Obj U, FLA_Obj b )
{
  FLA_Obj UTL,   UTR,      U00,  u01,       U02, 
          UBL,   UBR,      u10t, upsilon11, u12t,
                           U20,  u21,       U22;

  FLA_Obj bT,              b0,
          bB,              beta1,
                           b2;

  FLA_Part_2x2( U,    &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_BR );

  FLA_Part_2x1( b,    &bT, 
                      &bB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( UBR ) < FLA_Obj_length( U ) ){

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00,  &u01,       /**/ &U02,
                                                &u10t, &upsilon11, /**/ &u12t,
                        /* ************* */   /* **************************** */
                           UBL, /**/ UBR,       &U20,  &u21,       /**/ &U22,
                           1, 1, FLA_TL );

    FLA_Repart_2x1_to_3x1( bT,                &b0, 
                                              &beta1, 
                        /* ** */            /* ***** */
                           bB,                &b2,        1, FLA_TOP );

    /*------------------------------------------------------------*/
    Dot_Gaussian_unb(beta1,u12t,b2);
    *(double*)FLA_Obj_buffer_at_view(beta1)/=*(double*)FLA_Obj_buffer_at_view(upsilon11);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00,  /**/ u01,       U02,
                            /* ************** */  /* ************************** */
                                                     u10t, /**/ upsilon11, u12t,
                              &UBL, /**/ &UBR,       U20,  /**/ u21,       U22,
                              FLA_BR );

    FLA_Cont_with_3x1_to_2x1( &bT,                b0, 
                            /* ** */           /* ***** */
                                                  beta1, 
                              &bB,                b2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}
//
void matrix_generate( FLA_Obj A, char type,char name) {
  double  * buff_A;
  int  m_A, n_A, ldim_A;
  int  num;
  double xin;
  
  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );
  num = 1;
  if(type=='F'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = 0; i < m_A; i++ ) {
	printf("Please enter %c(%d,%d) \t",name,i,j);
	scanf("%lf",&buff_A[ i + j * ldim_A ]);
        //buff_A[ i + j * ldim_A ] = ( double ) num;
        //num++;
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
