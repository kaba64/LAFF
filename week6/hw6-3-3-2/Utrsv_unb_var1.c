#include <stdio.h>
# include <sys/time.h>
#include <FLAME.h>
#include "dot-gaussian.h"

#define PRINT_DATA
// Declaration of local prototypes.
void matrix_generate( FLA_Obj A, char type,char name);
double spendtimeinsecond();
int Utrsv_unb_var1( FLA_Obj U, FLA_Obj b );

int main( int argc, char *argv[] ) {
  int m_U,n_U;
  FLA_Obj U,b;
  double t1,t2;
  printf("Please enter the size of the matrix U : m_U and n_U \n");
  scanf("%d%d",&m_U,&n_U);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_U,n_U,0,0,&U);
  FLA_Obj_create(FLA_DOUBLE,m_U,1,0,0,&b);
  matrix_generate(U,'F','U');
  matrix_generate(b,'F','b');
#ifdef PRINT_DATA
  FLA_Obj_show( " U = [ ", U, "%le", " ];" );
  FLA_Obj_show( " b = [ ", b, "%le", " ];" );
#endif
  t1=spendtimeinsecond();
  Utrsv_unb_var1(U,b);
  t2=spendtimeinsecond();
#ifdef PRINT_DATA
  FLA_Obj_show( " U = [ ", U, "%le", " ];" );
  FLA_Obj_show( " b = [ ", b, "%le", " ];" );
#endif
  printf("The spent time for computing Back substitution in second : %lf\n",t2-t1);
  FLA_Obj_free(&U);
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
