#include <stdlib.h>
#include <stdio.h>
#include <FLAME.h>

#define PRINT_DATA

// Declaration of local prototypes.                            

static void matrix_generate( FLA_Obj A );
int AXPY_unb(FLA_Obj alpha,FLA_Obj x,FLA_Obj y);

int main( int argc, char *argv[] ) {
  int      m_A, n_A;
  FLA_Obj  x, y,rho;
  double  * rho_ass;
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&rho);
  rho_ass = ( double * ) FLA_Obj_buffer_at_view(rho);
  FLA_Init();
  printf("Please enter alpha : \n");
  scanf("%lf", rho_ass);
  printf("Please the dimension of the matrix m_a and n_a : \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Obj_create(FLA_DOUBLE,m_A,1,0,0,&x); //FLA_INT   
  FLA_Obj_create(FLA_DOUBLE,m_A,1,0,0,&y);
  matrix_generate(x);
  matrix_generate(y);
  //matrix_generate(rho);
#ifdef PRINT_DATA
  FLA_Obj_show( " pi = [ ", x, "%le", " ];" );
  FLA_Obj_show( " taui = [ ", y, "%le", " ];" );
  FLA_Obj_show( " rho = [ ", rho, "%le", " ];" );
#endif
  AXPY_unb(rho, x,y);
#ifdef PRINT_DATA
  FLA_Obj_show( " pi = [ ", x, "%le", " ];" );
  FLA_Obj_show( " taui = [ ", y, "%le", " ];" );
#endif
  printf( "%% End of Program\n" );
  FLA_Finalize();
  return 0;
}
/* Copyright 2022 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: Kazem Bazesefidpar
                  kazemba@kth.se
*/
int AXPY_unb( FLA_Obj alpha, FLA_Obj x, FLA_Obj y )
{
  FLA_Obj xT,              x0,
          xB,              chi1,
                           x2;

  FLA_Obj yT,              y0,
          yB,              psi1,
                           y2;

  double* alphain, chi1in, psi1in;
  alphain = ( double * ) FLA_Obj_buffer_at_view(alpha);
  
  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_TOP );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_length( xT ) < FLA_Obj_length( x ) ){

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                        /* ** */            /* **** */
                                              &chi1, 
                           xB,                &x2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* **** */
                                              &psi1, 
                           yB,                &y2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    *(double*)FLA_Obj_buffer_at_view(psi1)+=*alphain*(*(double*)FLA_Obj_buffer_at_view(chi1));
    /*                       update line 1                        */
    /*                             :                              */
    /*                       update line n                        */

    /*------------------------------------------------------------*/

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
// https://github.com/flame/hqrrp/blob/master/libflame_sources/simple_test.c   
// ============================================================================  
static void matrix_generate( FLA_Obj A ) {
  double  * buff_A;
  int     m_A, n_A, ldim_A;
  int     i, j, num;

  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );
  // 
  // Matrix with integer values.   
  // ---------------------------  
  //  
  if( ( m_A > 0 )&&( n_A > 0 ) ) {
    num = 1;
    for ( j = 0; j < n_A; j++ ) {
      for ( i = ( j % m_A ); i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
      for ( i = 0; i < ( j % m_A ); i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
    }
    if( ( m_A > 0 )&&( n_A > 0 ) ) {
      buff_A[ 0 + 0 * ldim_A ] = 1.2;
    }
  }
}
