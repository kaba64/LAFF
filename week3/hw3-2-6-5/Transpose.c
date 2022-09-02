#include <stdio.h>
#include <FLAME.h>
#define PRINT_DATA

// Declaration of local prototypes.
void matrix_generate( FLA_Obj A );
void CopySymmetrize(double* xt,double* x,int n,int ldim,char type);
int Symmetrize_tri_unb( FLA_Obj A ,char type);

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A;
  
  printf("Please enter the size of the matrix m_A and n_A: \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  matrix_generate(A);
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
#endif
  Symmetrize_tri_unb(A,'U');
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
#endif
  printf("End of Program\n");
  FLA_Obj_free(&A);
  FLA_Finalize();
  return 0;
}
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int Symmetrize_tri_unb( FLA_Obj A, char type )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;
  double* buff_a01;
  double* buff_a10t;
  int n_a01t,ldim_A;

  ldim_A = FLA_Obj_col_stride(A);
  
  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/
    buff_a10t  =  (double *) FLA_Obj_buffer_at_view(a10t);
    n_a01t     =             FLA_Obj_width(a10t);
    buff_a01   =  (double *) FLA_Obj_buffer_at_view(a01);
    CopySymmetrize(buff_a10t,buff_a01,n_a01t,ldim_A,type);
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

void CopySymmetrize(double* xt,double* x,int n,int ldim,char type){
  int i;
  if(type=='U'){
    for(i=0;i<n;++i){
      xt[i*ldim]=x[i];
    }
  }else if(type=='L'){
    for(i=0;i<n;++i){
      x[i]=xt[i*ldim];
    }
  }
}
