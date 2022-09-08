#include <stdio.h>
# include <sys/time.h>
#include <FLAME.h>
#include "dot.h"
#include "Dot_Both_column_vector.h"
#define PRINT_DATA

// Declaration of local prototypes.
void matrix_generate( FLA_Obj A, char type);
double spendtimeinsecond();
int Symv_u_unb_var1( FLA_Obj A, FLA_Obj x, FLA_Obj y );

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
  matrix_generate(A,'S');
  matrix_generate(x,'F');
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " X = [ ", x, "%le", " ];" );
  FLA_Obj_show( " Y = [ ", y, "%le", " ];" );
#endif
  t1=spendtimeinsecond();
  Symv_u_unb_var1(A,x,y);
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
int Symv_u_unb_var1( FLA_Obj A, FLA_Obj x, FLA_Obj y )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj xT,              x0,
          xB,              chi1,
                           x2;

  FLA_Obj yT,              y0,
          yB,              psi1,
                           y2;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_TOP );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                        /* ** */            /* **** */
                                              &chi1, 
                           xB,                &x2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* **** */
                                              &psi1, 
                           yB,                &y2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    Dot_Both_column_vector_unb(psi1,a01,x0);
    Dot_unb(psi1,alpha11,chi1);
    Dot_unb(psi1,a12t,x2);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &xT,                x0, 
                                                  chi1, 
                            /* ** */           /* **** */
                              &xB,                x2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                                                  psi1, 
                            /* ** */           /* **** */
                              &yB,                y2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}
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

double spendtimeinsecond()
{
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
