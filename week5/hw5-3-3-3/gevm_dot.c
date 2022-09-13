#include <FLAME.h>
#include "dot.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int gevm_dot_unb( FLA_Obj a, FLA_Obj B, FLA_Obj ct )
{
  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Obj cLt,    cRt,       c0t,  gamma1,  c2t;

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  FLA_Part_1x2( ct,    &cLt,  &cRt,      0, FLA_LEFT );

  while ( FLA_Obj_width( BL ) < FLA_Obj_width( B ) ){

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &b1, &B2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( cLt,  /**/ cRt,        &c0t, /**/ &gamma1, &c2t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    Dot_unb(gamma1,a,b1);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, b1, /**/ B2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &cLt,  /**/ &cRt,        c0t, gamma1, /**/ c2t,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
