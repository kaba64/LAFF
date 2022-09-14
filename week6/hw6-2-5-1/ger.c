#include "FLAME.h"
#include "axpy-gaussian-rank1.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int ger_unb( FLA_Obj x, FLA_Obj yt, FLA_Obj A )
{
  FLA_Obj yLt,    yRt,       y0t,  psi1,  y2t;

  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Part_1x2( yt,    &yLt,  &yRt,      0, FLA_LEFT );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  while ( FLA_Obj_width( yLt ) < FLA_Obj_width( yt ) ){

    FLA_Repart_1x2_to_1x3( yLt,  /**/ yRt,        &y0t, /**/ &psi1, &y2t,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    AXPY_Gaussian_Rank1_unb(psi1,x,a1);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &yLt,  /**/ &yRt,        y0t, psi1, /**/ y2t,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
