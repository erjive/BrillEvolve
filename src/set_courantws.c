#include "stdio.h"
#include "stdlib.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
//#include "cctk_Functions.h"
#include "math.h"

void setcourantwavespeed(CCTK_ARGUMENTS);

void setcourantwavespeed(CCTK_ARGUMENTS)
{
	DECLARE_CCTK_ARGUMENTS;
	DECLARE_CCTK_PARAMETERS;

  //Auxiliar variables
	int i,j,k;
	int index;
	int istart, jstart, kstart, iend, jend, kend;
  CCTK_REAL detg,guxx,guyy,guzz;
	CCTK_REAL vxmax,vymax,vzmax;

  //Reduction variables

  int vxindex,vyindex,vzindex;   /* grid variable index */
  int target_proc;               /* processor to hold the result */
  int reduction_handle;          /* handle for reduction operator */
  char *reduction_name;          /* reduction operator to use */
  int ierr = -1 ;                /* error flag*/

	istart = 0;
	jstart = 0;
	kstart = 0;

	iend = cctk_lsh[0];
	jend = cctk_lsh[1];
	kend = cctk_lsh[2];

  //printf("calculating gauge speed \n");
   #pragma omp parallel for collapse(3) private(i,j,k,index,detg,guxx,guyy,guzz)

  // #pragma omp private(i,j,k,index,detg,guxx,guyy,guzz,vxx,vyy,vzz)


	for (k=kstart; k<kend; k++)
	{
		for (j=jstart; j<jend; j++)
		{
			for (i=istart; i<iend; i++)
			{
       index = CCTK_GFINDEX3D(cctkGH,i,j,k);

       detg = -(pow(gxz[index],2)*gyy[index])+ 2.0*gxy[index]*gxz[index]*gyz[index]
              - gxx[index]*pow(gyz[index],2) - pow(gxy[index],2)*gzz[index]
              + gxx[index]*gyy[index]*gzz[index];

       //printf("detg=%f\n",detg);

       guxx = (-pow(gyz[index],2) + gyy[index]*gzz[index])/detg;

       if (CCTK_Equals (initial_data, "brilldata3D"))
       {
       guyy = (-pow(gxz[index],2) + gxx[index]*gzz[index])/detg;
       }
       guzz = (-pow(gxy[index],2) + gxx[index]*gyy[index])/detg;

       if (gauge_shock_avoid != 0)
       {
          vxx[index]=sqrtf(fabs((alp[index]*alp[index]+kappa_shock)*guxx));

          if (CCTK_Equals (initial_data, "brilldata3D")){
            vyy[index]=sqrtf(fabs((alp[index]*alp[index]+kappa_shock)*guyy));
          }

          vzz[index]=sqrtf(fabs((alp[index]*alp[index]+kappa_shock)*guzz));

/*          if (isnan(vxx[index])){
            vxx[index] = 0.0;
          }
          if (isnan(vyy[index])){
            vyy[index] = 0.0;
          }
          if (isnan(vzz[index])){
            vzz[index] = 0.0;
          }
*/
          //printf("vxx=%f,vyy=%f,vzz=%f \n",vxx[index],vyy[index],vzz[index]);
       }

			}
    }
  }
  //printf("all gauge speed calculus done \n");

  /* want to get the maximum for the wavetoy grid function */
  reduction_name = "maximum";

  vxindex = CCTK_VarIndex ("BrillEvolve::vxx");
  vyindex = CCTK_VarIndex ("BrillEvolve::vyy");
  vzindex = CCTK_VarIndex ("BrillEvolve::vzz");
  /* the reduction result will be obtained by processor 0 only */
  target_proc = 0;
  /* get the handle for the given reduction operator */
  reduction_handle = CCTK_ReductionHandle (reduction_name);

//  if (reduction_handle >= 0)
//  {
  ierr = CCTK_Reduce (cctkGH, target_proc, reduction_handle,1, CCTK_VARIABLE_REAL, &vxmax, 1, vxindex);
  ierr = CCTK_Reduce (cctkGH, target_proc, reduction_handle,1, CCTK_VARIABLE_REAL, &vzmax, 1, vzindex);

  if (CCTK_Equals (initial_data, "brilldata3D"))
  {
    ierr = CCTK_Reduce (cctkGH, target_proc, reduction_handle,1, CCTK_VARIABLE_REAL, &vymax, 1, vyindex);
  }

     if ( ierr == 0)
     {

        //printf("vxmax=%f,vymax%f,vzmax%f \n",vxmax,vymax,vzmax);
        if (CCTK_Equals (initial_data, "brilldata3D"))
        {
          *courant_wave_speed = fmax(fmax(vxmax,vymax),vzmax);
        }
        else if (CCTK_Equals (initial_data, "brilldata2D"))
        {
          //printf("vxmax=%f,vzmax%f \n",vxmax,vzmax);
          *courant_wave_speed = fmax(vxmax,vzmax);
        }

          //printf("COURANT WAVE SPEED = %f \n",*courant_wave_speed);
     }

     else
     {

         CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,"%s reduction failed", reduction_name);
     }

//  }

}
