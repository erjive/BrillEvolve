/*@@
  @file      brilldata.F
  @date
  @author    Carsten Gundlach, Miguel Alcubierre.
  @desc
             Construct Brill wave q function.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*#include "cctk_Functions.h"*/
#include "math.h"

  CCTK_REAL brillqshape(CCTK_REAL rho1, CCTK_REAL z1)
  /* Calculates the function q that appear in the conformal
     metric for Brill waves:

     ds^2  =  psi^4 ( e^(2q) (drho^2 + dz^2) + rho^2 dphi^2 )

    There are different choices for the form of q depending on
    the value of the parameter "brill_q":

    For the moment I will only consider the Gundlach (Holz) initial data:

    brill_q=2:
                                           2    2      2  c/2
                                    - [ ( r - r0 ) / sr  ]
                  q = a (rho/srho)  e
  */

  {
/*	    DECLARE_CCTK_ARGUMENTS;*/
	    DECLARE_CCTK_PARAMETERS;
/*	    DECLARE_CCTK_FUNCTIONS;*/

      CCTK_REAL q;

      if (CCTK_EQUALS(q_function,"gundlach"))
      {

        if (rho1 < 0)
        {
          rho1 = -rho1;
        }

         q = gundlacha*pow((rho1/gundlachsrho),2)*
                    expf(-(pow(rho1,2)+pow(z1,2)-pow(gundlachr0,2))/
                    pow(gundlachsigmar,2));


      }

      return q;
  }

