#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*#include "cctk_Functions.h"*/
#include "math.h"
#include "qshape.h"
/*#include "CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h"*/

/* Adapted from the IDBrillData thorn */


void brilldata(CCTK_ARGUMENTS);


void brilldata(CCTK_ARGUMENTS)
{

	DECLARE_CCTK_ARGUMENTS;
	DECLARE_CCTK_PARAMETERS;
/*	DECLARE_CCTK_FUNCTIONS;*/

	int i,j,k;
	int index;
	int istart, jstart, kstart, iend, jend, kend;
	/*int ierr=-1;*/
	/* char inputfile =""
	   char inputvar  = */

	CCTK_REAL x1,y1,z1,rho1,rho2;
	CCTK_REAL e2q,psi4;
	CCTK_REAL zero,one;

	/*CCTK_REAL brillq;*/

	zero = 0.0;
	one  = 1.0;

	istart = 0;
	jstart = 0;
	kstart = 0;

	iend = cctk_lsh[0];
	jend = cctk_lsh[1];
	kend = cctk_lsh[2];

	/* Read the initial data */

	/*ierr = IOUtil_RecoverVarsFromDatafiles (GH,IDfile,"BrillEvolve::brillpsi{ alias=’IDBrillMoL::phi’ }");*/

/*	if (ierr<0)
	{
		CCTK_WARN(0, "Brill wave initial data not found!");
	}
*/
#ifdef DEBUG_MOL
	printf("About to loop.\n");
#endif

	for (k=kstart; k<kend; k++)
	{
		for (j=jstart; j<jend; j++)
		{
			for (i=istart; i<iend; i++)
			{
				index = CCTK_GFINDEX3D(cctkGH,i,j,k);

				x1 = x[index];
				y1 = y[index];
				z1 = z[index];

				rho2 = x1*x1+y1*y1;
				rho1 = sqrt(rho2);

				e2q  = expf(2.0*brillqshape(rho1,z1));

				/*  Fudge division by rho^2 on axis. (Physically, y^/rho^2,
						x^2/rho^2 and xy/rho^2 are of course regular.)
						Transform Brills form of the physical metric to Cartesian
						coordinates via

						e^2q (drho^2 + dz^2) + rho^2 dphi^2 =
						e^2q (dx^2 + dy^2 + dz^2) + (1-e^2q) (xdy-ydx)^2/rho^2

						The individual coefficients can be read off as
						*/

				if (rho2>0.0)
				{
					gxx[index] = (e2q + (one - e2q)*y1*y1/rho2);
					gyy[index] = (e2q + (one - e2q)*x1*x1/rho2);
					gzz[index] = e2q;
					gxy[index] = - (one - e2q)*x1*y1/rho2;
				}
				else
				{
          //printf("x=%f,y=%f,z=%f \n",x1,y1,z1);
					gxx[index] = one;
					gyy[index] = one;
					gzz[index] = one;
					gxy[index] = zero;
				}

				gxz[index] = zero;
				gyz[index] = zero;

				kxx[index] = zero;
				kyy[index] = zero;
				kzz[index] = zero;
				kxy[index] = zero;
				kxz[index] = zero;
				kyz[index] = zero;
				/*printf("%f %f %f %f %s\n",x1,y1,z1,brillsource[index],qfunction);*/

			}
		}
	}

	if (CCTK_EQUALS(metric_type,"static conformal"))
	{
		*conformal_state = 1;

		for (k=kstart; k<kend; k++)
		{
			for (j=jstart; j<jend; j++)
			{
				for (i=istart; i<iend; i++)
				{
					index = CCTK_GFINDEX3D(cctkGH,i,j,k);

					psi[index] = brillpsi[index];

				}

			}
			printf("psi_val=%f\n",psi[index]);

		}

		printf("Conformal factor filled, static conformal\n");

	}
	else{

		*conformal_state = 0;

			for (k=kstart; k<kend; k++)
			{
				for (j=jstart; j<jend; j++)
				{
					for (i=istart; i<iend; i++)
					{
						index = CCTK_GFINDEX3D(cctkGH,i,j,k);

						psi4  = pow(brillpsi[index],4);

						gxx[index] = psi4*gxx[index];
						gyy[index] = psi4*gyy[index];
						gzz[index] = psi4*gzz[index];
						gxy[index] = psi4*gxy[index];

					}

				}

			}
		printf("Conformal factor filled, physical\n");
	}


}
