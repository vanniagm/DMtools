/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      

#define OMEGA      
#define CDM_NUCLEON      
#define NEUTRINO
#define DECAYS
  /* Calculate relic density and display contribution of  individual channels */

#define CLEAN  to clean intermediate files

/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"





int main(int argc,char** argv)
{  int err,i,j,k,l,m,n,ll;
   char cdmName[10];
   char *mus="mus", *lsh="lsh";
   double ms,lam,ls;
   int spin2, charge3,cdim;
   int IN,JN,KNN,LN;
   double Omegah2Planck,Omegah2PlanckE;
   double XSMAX,XSMIN,LamMAX,pbcm2;
   int fast=1;
   double Beps=1.E-5, cut=0.01;
   double Omega,Xf;  
   double pA0[2],pA5[2],nA0[2],nA5[2],sigtot;
   double Nmass=0.939,Nmass2=122.,Nmassp=.938; /*nucleon mass*/
   double SCcoeff,SCcoeff2,SCcoeffp; 
   double nu[NZ], nu_bar[NZ],mu[NZ];
   double Ntot,Nutot;
   int forSun=1;
   int count,INDEX;
   double Emin=1;
   char fname[100];
   FILE * dat;


  
   dat=fopen ("Data1.dat","w");
   pbcm2 = 1.E-36;
   Omegah2Planck  = 0.1198;
   Omegah2PlanckE = 0.0026;
   IN=100;
   KNN=100;
   LN=0;

   ForceUG=0;  /* to Force Unitary Gauge assign 1 */

        for(i=0;i<IN;i++)
                {		 
                
	        ls=i*.01/(IN-1)+.0005; 
                assignValW(lsh,ls); 
                
		          for(k=0;k<KNN;k++)
		                {
		                ms=k*35./(KNN-1.)+45.;
		                assignValW(mus,ms);
		   
/*-----------------------------------MicrOmegasfunctions---------------------------------------------------*/
				  err=sortOddParticles(cdmName);
				  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}			  
				   qNumbers(cdmName, &spin2, &charge3, &cdim);
				 
				#ifdef OMEGA
				// to exclude processes with virtual W/Z in DM   annihilation      
				  VZdecay=0; VWdecay=0; cleanDecayTable();
				  Omega=darkOmega(&Xf,fast,Beps);
				   #endif				  
			      if (Omegah2Planck-3.*Omegah2PlanckE<=Omega&& Omega<= Omegah2Planck+3.*Omegah2PlanckE)
			      { 
			        
				  #ifdef CDM_NUCLEON        
				    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
				  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.)*pbcm2;
				   sigtot=SCcoeff*pow(54.*pA0[0]+(131.-54.)*nA0[0],2.)/pow(131.,2.);				  
				#endif
				
				   
				   fprintf(dat,"%.5E ",ms);
				   fprintf(dat,"%.5E ",ls);
				   fprintf(dat,"%.5E ",Omega);
				   fprintf(dat," %.5E  %.5E  %.5E %.5E %.5E \n",
				   SCcoeff*pow(nA0[0],2.),SCcoeff*pA0[0]*pA0[0],sigtot,
				   3.*SCcoeff*pow(nA5[0],2.),3.*SCcoeff*pA5[0]*pA5[0]); 
				#ifdef CLEAN
				#endif
	                            
			 	}
				   
				 /* killPlots();*/
									    
				
                        } /*ms*/
                }
 return 0;
}
