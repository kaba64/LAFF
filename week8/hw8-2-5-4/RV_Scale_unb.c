#include <FLAME.h>
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int RV_Scale_unb( FLA_Obj alpha, FLA_Obj xt )
{
  FLA_Obj xLt,    xRt,       x0t,  chi1,  x2t;

  double alphain;
  alphain = *( double * ) FLA_Obj_buffer_at_view(alpha);
  
  FLA_Part_1x2( xt,    &xLt,  &xRt,      0, FLA_LEFT );

  while ( FLA_Obj_width( xLt ) < FLA_Obj_width( xt ) ){

    FLA_Repart_1x2_to_1x3( xLt,  /**/ xRt,        &x0t, /**/ &chi1, &x2t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    *(double*)FLA_Obj_buffer_at_view(chi1)=alphain*(*(double*)FLA_Obj_buffer_at_view(chi1));
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &xLt,  /**/ &xRt,        x0t, chi1, /**/ x2t,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
