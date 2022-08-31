#include <stdlib.h>
#include <stdio.h>
#include <FLAME.h>

#define PRINT_DATA

// Declaration of local prototypes.
int Set_to_diagonal_matrix_unb(FLA_Obj A,FLA_Obj X);
void SetZeroColumnVector(int* x,int n);
void InitializeVector(FLA_Obj x, char type);

int main( int argc, char *argv[] ) {
  int m_A, n_A;
  FLA_Obj  A,x;
  
  printf("Please enter the size of the matrix m_A and n_A: \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_INT,m_A,n_A,0,0,&A);
  FLA_Obj_create(FLA_INT,n_A,1,0,0,&x);
  InitializeVector(x,'I');
  Set_to_diagonal_matrix_unb(A,x);
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%d", " ];" );
#endif
  printf("End of Program\n" );
  FLA_Finalize();
  return 0;
}
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Set_to_diagonal_matrix_unb( FLA_Obj A, FLA_Obj X )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj XT,              X0,
          XB,              x1t,
                           X2;
  int* buff_a01;
  int* buff_alpha;
  int* buff_a21;
  int n_a01,n_a21;
  
  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( X,    &XT, 
                      &XB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( XT,                &X0, 
                        /* ** */            /* *** */
                                              &x1t, 
                           XB,                &X2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    buff_a01   = (int*) FLA_Obj_buffer_at_view(a01);
    n_a01      = FLA_Obj_length (a01);
    SetZeroColumnVector(buff_a01,n_a01);
    buff_alpha = ( int * ) FLA_Obj_buffer_at_view(alpha11);
    *buff_alpha = *(( int * ) FLA_Obj_buffer_at_view(x1t));
    buff_a21   = (int*) FLA_Obj_buffer_at_view(a21);
    n_a21      = FLA_Obj_length(a21);
    SetZeroColumnVector(buff_a21,n_a21);
    /*------------------------------------------------------------*/
    
    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &XT,                X0, 
                                                  x1t, 
                            /* ** */           /* *** */
                              &XB,                X2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}

void SetZeroColumnVector(int* x,int n){
  int i;
  for(i=0;i<n;++i){
    x[i]=0;
  }
}

void InitializeVector(FLA_Obj x,char type){
  //double* buff_x;
  int m_x, n_x;
  int i,num;
  num=0;
  if(type=='I'){
    int* buff_x;
    buff_x = ( int * ) FLA_Obj_buffer_at_view( x );
    m_x    = FLA_Obj_length( x );
    n_x    = FLA_Obj_width ( x );
    if(n_x==1){
      for(i=0;i<m_x;++i){
	buff_x[i]=(int) num;
	++num;
      }
    }else if(m_x==1){
      for(i=0;i<n_x;++i){
	buff_x[i]=(int) num;
        ++num;
      }
    }
  }else if(type=='D'){
    double* buff_x;
    buff_x = ( double * ) FLA_Obj_buffer_at_view( x );
    m_x    = FLA_Obj_length( x );
    n_x    = FLA_Obj_width ( x );
    if(n_x==1){
      for(i=0;i<m_x;++i){
        buff_x[i]=(double) num;
        ++num;
      }
    }else if(m_x==1){
      for(i=0;i<n_x;++i){
        buff_x[i]=(double) num;
        ++num;
      }
    }
  }
}
