// Code from Chris Hans: http://www.stat.osu.edu/~hans/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

extern void rnorm_truncated(double*, int*, double*, double*, double*, double*);

extern "C"
{

#include "util_blasso.h"

// A Gibbs sampler for the Bayesian Lasso
void blassoGibbs(double* X, int* nn, int* pp, int* TT, int* BB, int* tthin,
   	            double* ttau, double* ssig2, double* pphi, 
			  double* sig2prior, int* fits2,
			  double* tauprior, int* fittau, int* modunc,
			  double* phiprior, int* fitphi, double* start, 
			  double* bdraws, double* sig2draws, 
			  double* taudraws, double* phidraws, double* marginc,
			  int* rrb, double* RB, double* YtY, double* YtX, 
			  int* NOISY, int* bprior, int* count)
{
   int		n = *nn; 			// no. observations
   int		p = *pp;			// no. variables
   int		T = *TT;			// no. Gibbs draws
   int		rb = *rrb;		// should we rao-blackwellize? 1(0)
   int		thin = *tthin;		// thinning the chain
   int		B = *BB;			// length of Burn in
   int		noisy = *NOISY;	// print progress
   double		tau = *ttau;		// penalty term
   double		sig2 = *ssig2;		// fixed error variance (or start value)
   double		phi = *pphi;		// prior variable incl. prob.
   double		post_shape=0.0, post_rate=0.0;	// posterior for sig2
   double		tau_rate=0.0, tau_shape=0.0;		// posterior for tau
   int		kvar=0;			// nvar in current model
   double		pscale;			// used to accomodate different
   							// priors on beta
   double		L1;				// l1 norm of betas
   int         i, j, k, tt, perc=2; // looping variables

   int		iter;
   int		nonzero;			// for var selection
   int		nvar;			// number of active predictors
   double		w;				// prob >< 0
   double		*beta;			// vector of current beta values
   double		mucon_pos, mucon_neg;  // conditional ``means'' for positive
   							// and negative components
   double		ldnpos, ldnneg, lpnpos, lpnneg; // normal dist stuff
   double		u;				// uniform RV
   int		tmpcnt;			// how many proposals for sig2 sampling
   int		badsamp;			// did sig2 sampler work?
   int		range;			// how far out should we go for sig2
   							// sampler (ie, how many cut points)?

   // for sampling from truncated normal
   int		one=1;
   double		sd, zero = 0.0;
   double		PosInf = R_PosInf, NegInf = R_NegInf;	// constants

   if(*fits2)
   {
	 post_shape = sig2prior[0] + 0.5*(double)n;

	 // add the extra bit for the different prior
	 if((*bprior) == 1) post_shape += 0.5*(double)p;
   }

   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////
   // Find the XtX matrix

   double* XtX = new double[p*p];

   for(i=0; i<p; i++)
   {
      for(j=0; j<p; j++)
      {
	    XtX[i*p + j] = 0.0;
	    for(k=0; k<n; k++) XtX[i*p + j] += X[i*n + k]*X[j*n + k];
      }

      RB[i] = 0.0; // the RB means...
   }

   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////

   // Set the starting values
   beta = new double[p];
   for(j=0; j<p; j++)
   {
      beta[j] = start[j];
	 if(*modunc && (beta[j]!=0)) kvar++;
   }

   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////
   //			START THE SAMPLER!!!			 //
   ////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////

   int progout = (thin*(T+B)*p) >= 40000;

   if(noisy && progout) // percent done
   {
      Rprintf("\nSampler Progress...\n| ");
	 for(i=0; i<100; i+=20) Rprintf("%d%%      ",i);
	 Rprintf("%d%% |\n|",i);
   }

   for(iter=0; iter < (T+B); iter++)
   {
      // Possibly thin the chain
      for(tt=0; tt<thin; tt++)
      {
         // allow for various priors on beta
	    if((*bprior) == 0) pscale = 1.0;
	    else pscale = sqrt(sig2);

	    L1 = 0.0;
	    nvar = 0; // reset at each iter

         // Loop over the full conditionals for beta
         for(j=0; j < p; j++)
         {
            // Conditional sd
		  sd = sqrt(sig2/XtX[j*p + j]);

            // Compute the conditional ``means''
		  mucon_pos = YtX[j]/XtX[j*p + j] - tau*sd*sd/pscale;
		  mucon_neg = YtX[j]/XtX[j*p + j] + tau*sd*sd/pscale;

            // Add the remaining piece
	       for(i=0; i<p; i++)
	       {
               if(i!=j)
	          {
	             mucon_pos -= beta[i]*XtX[i*p + j]/XtX[j*p + j];
	             mucon_neg -= beta[i]*XtX[i*p + j]/XtX[j*p + j];
	          }
	       }

		  ldnpos = dnorm(0.0, mucon_pos, sd, 1);
		  ldnneg = dnorm(0.0, mucon_neg, sd, 1);
		  lpnpos = pnorm(mucon_pos/sd, 0.0, 1.0, 1, 1);
		  lpnneg = pnorm(-mucon_neg/sd, 0.0, 1.0, 1,1);

            // nonzero may get set to 0 below if we're doing
		  // variable selection
		  nonzero = 1;
		  if(*modunc)
		  {
		     // weight in favor of = 0

               // If using the scaled prior (*bprior==1), 
			// the term is tau/(2*sigma). Otherwise the
			// term is just tau/2

			if((*bprior) == 1) w = 1.0/(1.0 + exp(log(phi) - log(1.0-phi) + log(tau) - log(2.0) - 0.5*log(sig2) + log(exp(lpnpos - ldnpos) + exp(lpnneg - ldnneg))));
			else w = 1.0/(1.0 + exp(log(phi) - log(1.0-phi) + log(tau) - log(2.0) + log(exp(lpnpos - ldnpos) + exp(lpnneg - ldnneg))));

               // Marginal inclusion probability
			//if((iter>=B) && (tt==(thin-1))) marginc[(iter-B)*p+j] = 1.0-w;
			if((iter>=B) && (tt==(thin-1))) marginc[j] += (1.0 - w) / (double)T;

			// CONTROL MODEL SIZE HERE!
			// If there are currently n variables in the model
			// and variable j is NOT already in the model,
			// force it to be zero
			//if((kvar==n) && (beta[j]==0.0)) w = 1.0;

			GetRNGstate();
			u = runif(0.0,1.0);
			PutRNGstate();

			if(u < w)
			{
			   // If beta[j] was currently nonzero, decrease
			   // the number of variables in the model
			   if(beta[j]!=0.0) kvar--;

                  beta[j] = 0.0;

			   nonzero = 0;
			}
			else
			{
			   // If beta[j] currently zero, increase
			   // the number of variables in the model
			   if(beta[j]==0.0) kvar++;
			}
		  }

            // If beta is nonzero for this draw
		  if(nonzero)
		  {
		     // Increment the number of active variables.
			// If fitting the full model, this will be
			// p at the end of the scan through beta
               nvar++; 

               // Compute the weight in favor of the nonnegative component
	          w = 1.0/(1.0 + exp(ldnpos - ldnneg + lpnneg - lpnpos));

               GetRNGstate();
	          u = runif(0.0,1.0);
	          PutRNGstate();

	          if(u < w)
		        rnorm_truncated(beta+j,&one,&mucon_pos,&sd,&zero,&PosInf);
	          else
		        rnorm_truncated(beta+j,&one,&mucon_neg,&sd,&NegInf,&zero);
		  }

	       // Update L1 norm
		  L1 += fabs(beta[j]);

            if((iter>=B) && (tt==(thin-1))) 
	       {
	          bdraws[(iter-B)*p + j] = beta[j];

	          // Compute rao-blackwellized mean if requested
	          if(rb==1) RB[j] += (w*(mucon_pos + sd*exp(dnorm(0.0,mucon_pos, sd, 1)-pnorm(mucon_pos/sd,0.0,1.0,1,1))) + (1-w)*(mucon_neg - sd*exp(dnorm(0.0, mucon_neg, sd, 1)-pnorm(-mucon_neg/sd, 0.0, 1.0, 1, 1))))/(double)T;
	       }
         }

	    // Update phi if so desired...
	    if(*modunc && *fitphi)
	    {
		  GetRNGstate();
		  phi = rbeta(phiprior[0] + kvar, phiprior[1] + p - kvar);
		  PutRNGstate();

            if((iter>=B) && (tt==(thin-1))) phidraws[iter-B] = phi;
	    }

	    // Sample sig2 if so desired...
	    if(*fits2)
	    {
            // This can all be made more efficient if we keep
		  // track of which variables are nonzero when doing
		  // variable selection...

            // compute the ``shape'' parameter
	       post_shape = sig2prior[0] + 0.5*(double)n;
            if((*bprior) == 1) post_shape += 0.5*(double)nvar;

            // compute the ``rate'' parameter
		  double tmp = *YtY;
		  if(nvar > 0)
		  {
		     for(j=0; j<p; j++) 
		     {
                  tmp += -2.0*YtX[j]*beta[j];
			   for(k=0; k<p; k++) tmp += beta[j]*XtX[j*p + k]*beta[k];
	          }
		  }
		  post_rate = sig2prior[1] + 0.5*tmp;

		  // If we are using the classic prior, just sample it...
		  if((*bprior)==0)
		  {
		     GetRNGstate();
		     sig2 = 1.0/rgamma(post_shape, 1.0/post_rate);
			PutRNGstate();
	       }
		  else // if using the rescaled prior
		  {
		     // rejection sampling here!
			// start with 3 sd's, increment if necessary
			range = 3;
			do
			{
			   badsamp = sig2_rej_samp(&sig2, post_shape, post_rate, 
			   					  tau, L1, range++, &tmpcnt);
			} while(badsamp);

			//if((iter>=B) && (tt==(thin-1))) count[iter-B] = tmpcnt;
		  }

		  if((iter>=B) && (tt==(thin-1))) sig2draws[iter-B] = sig2;
	    }

         // Sample tau if so desired
	    if(*fittau)
	    {
	       if((*bprior)==0) tau_rate = L1 + tauprior[1];
		  else tau_rate = L1/sqrt(sig2) + tauprior[1];

		  tau_shape = tauprior[0] + (double)nvar;

		  GetRNGstate();
		  tau = rgamma(tau_shape, 1.0/tau_rate);
		  PutRNGstate();

		  if((iter>=B) && (tt==(thin-1))) taudraws[iter-B] = tau;
	    }
      }

	 if(noisy && progout && ( (thin*(iter+1))/(double)(thin*(T+B)) >= perc/100.0))
	 {
         Rprintf("*");
	    perc += 2;
	 }
   }

   if(noisy && progout) Rprintf("|\n\n");

   delete[] beta; beta = NULL;
   delete[] XtX; XtX = NULL;

   return;
}

} // end extern "C"

