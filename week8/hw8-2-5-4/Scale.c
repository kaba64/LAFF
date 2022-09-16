#include <FLAME.h>
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Scale_unb( FLA_Obj alpha, FLA_Obj x )
{
  FLA_Obj xT,              x0,
          xB,              chi1,
                           x2;

  double alphain;

  alphain = *( double * ) FLA_Obj_buffer_at_view(alpha);
  
  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_TOP );

  while ( FLA_Obj_length( xT ) < FLA_Obj_length( x ) ){

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                        /* ** */            /* **** */
                                              &chi1, 
                           xB,                &x2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    *(double*)FLA_Obj_buffer_at_view(chi1)=alphain*(*(double*)FLA_Obj_buffer_at_view(chi1));
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &xT,                x0, 
                                                  chi1, 
                            /* ** */           /* **** */
                              &xB,                x2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}
