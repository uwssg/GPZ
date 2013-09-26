#include "gaussian_process_driver.h"
#include "goto_tools.h"
#include "eigen_wrapper.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>

enum{ib0,ik0,ib1,ik1,im0,im1};

double lnchoose(int i1, int i2){
///returns the natural logarithm of i1 choose i2
   
   int i,j,k,l;
   double lnum,ldenom,xx;
   
   if(i1<i2){
     return 0.0;
   }
   
   if(i2<i1/2){
     i2=i1-i2;
   }
   
   lnum=0.0;
   for(j=i1;j>i1-i2;j--){
     lnum+=log(double(j));
   }
   
   ldenom=0.0;
   for(j=1;j<=i2;j++){
     ldenom+=log(double(j));
   }
   
   xx=lnum-ldenom;
  
  return xx;
   
}

gpnoisy::gpnoisy(){
  initialized=0;
  dim=-1;
  kk=15;
  room=10000;
  roomstep=10000;

  lambda=0.001;
  
  srchtime=0.0;
  ittime=0.0;
  dettime=0.0;
  covtime=0.0;
  invtime=0.0;
  splittime=0.0;
  
  chose2=0;
  dofswitch=0;
  priorswitch=0;
  only1=0;
  prior_width=1.0;
  
  sprintf(statname,"generic_status_file.sav");
  
}

gpnoisy::~gpnoisy(){

  int i;
  
  if(initialized==1){
    delete [] fn;
    
    delete kptr;
    
    for(i=0;i<room;i++)delete [] noise[i];
    delete [] noise;
  }
  else printf("this is odd; you are deleting gp without initializing it\n");
  
}

void gpnoisy::set_lambda(double ll){
  lambda=ll;
}

double gpnoisy::p_lambda(){
  return lambda;
}

void gpnoisy::initialize(int pin,double **seed, double *seedfn,
                         double **seednoise){
			 
  int i,j,k,l;
  double *mins,*maxs;
  FILE *output;
  
  initialized=1;
  
  if(dim<0){
      printf("WARNING you cannot initialize gpnoisy until you set dim\n");
      exit(1);
  }
  
  mins=new double[dim];
  maxs=new double[dim];
  for(i=0;i<dim;i++){
      mins[i]=0.0;
      maxs[i]=1.0;
  }
  //these are not actual maximum and minimum values;
  //they are set to 0 and 1 so that the distance function
  //in kd_tree actually returns the Euclidean distance without
  //any special normalizations
  
  room=pin;
  fn=new double[pin];
  noise=new double*[pin];
  for(i=0;i<pin;i++)noise[i]=new double[dim];
  
  for(i=0;i<pin;i++){
    fn[i]=seedfn[i];
    for(j=0;j<dim;j++)noise[i][j]=seednoise[i][j];
  }
  
  //printf("statname is %s\n",statname);
  
  output=fopen(statname,"a");
  fprintf(output,"about to build tree\n");
  fclose(output);
  kptr=new kd_tree(dim,pin,seed,mins,maxs);
  kptr->check_tree(-1);
  
  output=fopen(statname,"a");
  fprintf(output,"tree diagnostic %d\n",kptr->diagnostic);
  fclose(output);
  
  if(kptr->diagnostic!=1){
      printf("WARNING kd_tree improperly built\n");
      exit(1);
  }
  
  pts=kptr->pts;
  
  delete [] mins;
  delete [] maxs;
  
}

void gpnoisy::add_pt(double *newpt, double newfn, double *newnoise){
  
  int i,j,k,l;
  double *buff,**nbuff;
  
  if(pts<room){
    fn[pts]=newfn;
    for(i=0;i<dim;i++){
      noise[pts][i]=newnoise[i];
    }
  }
  else{
    //printf("need to make room\n");
    buff=new double[pts];
    nbuff=new double*[pts];
    for(i=0;i<pts;i++){
      buff[i]=fn[i];
      nbuff[i]=new double[dim];
      for(j=0;j<dim;j++)nbuff[i][j]=noise[i][j];
      delete [] noise[i];
    }
    delete [] noise;
    delete [] fn;
    room+=roomstep;
    fn=new double[room];
    noise=new double*[room];
    for(i=0;i<room;i++)noise[i]=new double[dim];
    for(i=0;i<pts;i++){
      fn[i]=buff[i];
      for(j=0;j<dim;j++)noise[i][j]=nbuff[i][j];
      delete [] nbuff[i];
    }
    fn[pts]=newfn;
    for(i=0;i<dim;i++)noise[pts][i]=newnoise[i];
    delete [] buff;
    delete [] nbuff;
  }
  kptr->add(newpt);
  pts=kptr->pts;

}

void gpnoisy::get_pdf(
    double *pt,double *qnoise,
    double *xpdf, double *pdf,double *pdf_1, double *pdf_2,int npdf,
    double *chi1out, double *mu_1out, double *sig_1out,  
    double *chi2out, double *mu_2aout,double *sig_2aout, double *naout, 
                     double *mu_2bout, double *sig_2bout,double *nbout, 
		     double *ntotout, int skip){

  int i,j,k,l,n[2],splitdex,s[2],ii,jj;
  
  double mu[2],sig2[2],before,after,nn,xx1,xx2,choosefactor;
  double zbar[2];
  double sigbest,norm,sp1,sp2;
  double chi1,chi2,mu_1,sig2_1,chimin,llav;
  double zbar_1,sp_1,det1,det2;
  double deltab,prior,ddeltab1,ddeltab0,deltaz2;
  int ideltaz2;

  double mm_det[16],xx,yy;
  int rows[4],cols[4],ip,jp;
  
  double gg1det,gg2det,gg_1det;
  double prechi2,prechi1;
  
  double rawchi2,rawchi1,rawchi_1;
  double dx1db,dx2db,ddx1db,ddx2db;

   
  int oneisvalid,twoisvalid;

  
  FILE *output;
  
  
  int *neigh;
  neigh=new int[kk];
  
  double *dd,**ggq,**gg1,**gg1in;
  dd=new double[kk];
  ggq=new double*[2];
  ggq[0]=new double[kk];
  ggq[1]=new double[kk];
  gg1=new double*[kk];
  gg1in=new double*[kk];
  for(i=0;i<kk;i++){
      gg1[i]=new double[kk];
      gg1in[i]=new double[kk];
  }

  
  double **gg2,**gg2in;
  gg2=new double*[kk];
  gg2in=new double*[kk];
  for(i=0;i<kk;i++){
      gg2[i]=new double[kk];
      gg2in[i]=new double[kk];
  }
  

  double *tosort;
  tosort=new double[kk];
  
  double *ggq_1;
  ggq_1=new double[kk];
  
  double **gg_1,**gg_1in;
  gg_1=new double*[kk];
  gg_1in=new double*[kk];
  for(i=0;i<kk;i++){
      gg_1[i]=new double[kk];
      gg_1in[i]=new double[kk];
  }


    for(i=0;i<npdf;i++){
      pdf[i]=0.0;//the combined pdf
      
      pdf_1[i]=0.0;//the unimodal pdf
      
      pdf_2[i]=0.0;//the bimodal pdf
    }
    
    before=double(time(NULL));
   
    
    //find the neighbors;
    //skip>=0 denotes the index of a training galaxy
    //that is to be discarded from the model
    //(for use when optimizing on validation galaxies)
    
    if(skip<0){
        kptr->nn_srch(pt,kk,neigh,dd);
    }
    else{
        int *neighbuff;
	double *ddbuff;
	
	neighbuff=new int[kk+1];
	ddbuff=new double[kk+1];
        kptr->nn_srch(pt,kk+1,neighbuff,ddbuff);
	for(i=0,j=0;j<kk+1 && i<kk;j++){
	    neigh[i]=neighbuff[j];
	    dd[i]=ddbuff[j];
	    if(neighbuff[j]!=skip)i++;
        }
	
	if(i!=kk){
	    printf("WARNING did not get i to kk \n");
	    exit(1);
	}
	
	xx1=0.0;
	for(i=0;i<dim;i++){
	    xx1+=pt[i]*pt[i];
	}
	xx1=sqrt(xx1);
	
	for(i=0;i<kk;i++){
	    if(dd[i]/xx1<1.0e-10 || isnan(dd[i]/xx1)){
	        output=fopen(statname,"a");
	        fprintf(output,"BE AWARE\n");
		fprintf(output,"you seem to be validating using galaxies as their own neigbors\n");
		fprintf(output,"dd%d = %e -- %e\n\n",i,dd[i],xx1);
	        fclose(output);
	    }
	}
	delete [] neighbuff;
	delete [] ddbuff;
    }
    after=double(time(NULL));
    srchtime+=after-before;
    before=after;
    
    //calculate the characteristic length scale (llav)
    //which is used as a normalizing factor in the
    //covariogram
    ii=0;
    llav=0.0;
    for(i=0;i<kk;i++){
       for(j=i+1;j<kk;j++){
         if(ii>kk*(kk-1)/2-1){
            printf("WARNING ii overstepped in userpredict\n");
            exit(1);
         }
         xx1=0.0;
         for(k=0;k<dim;k++){
	     xx1+=power(kptr->data[neigh[i]][k]-kptr->data[neigh[j]][k],2);
	 }
         llav+=sqrt(xx1);
         ii++;
       }
    }   
    llav=llav/double(kk*(kk-1)/2);
     
 
    //divide the nearest neighbor galaxies along z;
    //split them so that the sum of the variances
    //in the two sub-populations is minimized
    
    double *sorted;
    sorted=new double[kk];
    for(i=0;i<kk;i++){
      tosort[i]=fn[neigh[i]];
    }
    //sort(tosort,neigh,kk);
    sort_and_check(tosort,sorted,neigh,kk);
    delete [] sorted;
    
    double var1,var2;
    k=-1;
    for(i=4;i<kk-4;i++){
      
      //testing possible divisions in the sorted population;
      //'i' represents the number of galaxies in the low-z
      //sub-population
      
      xx1=0.0;
      var1=0.0;
      xx2=0.0;
      var2=0.0;
      
      for(j=0;j<i;j++){
        xx1+=fn[neigh[j]];
        var1+=fn[neigh[j]]*fn[neigh[j]];
      }
      for(;j<kk;j++){
        xx2+=fn[neigh[j]];
	var2+=fn[neigh[j]]*fn[neigh[j]];
      }
      
      xx1=xx1/double(i);
      var1=var1/double(i-1)-double(i)*xx1*xx1/double(i-1);
      
      xx2=xx2/double(kk-i);
      var2=var2/double(kk-i-1)-double(kk-i)*xx2*xx2/double(kk-i-1);
    
      nn=var1+var2;
      
      if(k<0 || nn<sigbest){
        sigbest=nn;
	splitdex=i;
	k=1;
      }
      
    }//find the best splitting index
    
    //n[i] is the number of galaxies in the ith sub-population
    //s[i] is the starting index of the ith sub-population  
    n[0]=splitdex;
    n[1]=kk-splitdex;
    s[0]=0;
    s[1]=splitdex;
    
    after=double(time(NULL));
    splittime+=after-before;
    before=after;
    
    //calculate the covariance matrix elements relating the query point
    //to the galaxies in the 1st sub-population
    for(i=0;i<n[0];i++){
      j=s[0]+i;
      ggq[0][i]=covariogram(pt,qnoise,kptr->data[neigh[j]],noise[neigh[j]],llav,dim);
    }
    
    //ditto for the second subpopulation
    for(i=0;i<n[1];i++){
      j=s[1]+i;
      ggq[1][i]=covariogram(pt,qnoise,kptr->data[neigh[j]],noise[neigh[j]],llav,dim);
    }
    after=double(time(NULL));
    covtime+=after-before;
    before=after;
    
    //calculate the algebraic mean of z in each sub-population
    for(ii=0;ii<2;ii++){
      nn=0.0;
      zbar[ii]=0.0;
      for(i=s[ii];i<s[ii]+n[ii] && i<kk;i++)zbar[ii]+=fn[neigh[i]];
      zbar[ii]=zbar[ii]/double(n[ii]); 
    }
    
    after=double(time(NULL));
    ittime+=after-before;
    before=after;
    
    //calculate the covariance matrix relating the 1st sub-population galaxies
    //to each other
    for(i=0;i<n[0];i++){
      gg1[i][i]=1.0+lambda;
      k=s[0]+i;
      for(j=i+1;j<n[0];j++){
        l=s[0]+j;
        gg1[i][j]=covariogram(kptr->data[neigh[k]],noise[neigh[k]],\
	kptr->data[neigh[l]],noise[neigh[l]],llav,dim);
	
	gg1[j][i]=gg1[i][j];
      }
    }
    after=double(time(NULL));
    covtime+=after-before;
    before=after;
  
    //invert the gg1 covariance matrix;
    //gg1det contains the determinant of the matrix
    gg1det=fabs(invert_lapack(gg1,gg1in,n[0],1));
    after=double(time(NULL));
    invtime+=after-before;
    before=after;
    
    //repeat the process for the 2nd sub-population
    for(i=0;i<n[1];i++){
      gg2[i][i]=1.0+lambda;
      k=s[1]+i;
      for(j=i+1;j<n[1];j++){
        l=s[1]+j;
	
        gg2[i][j]=covariogram(kptr->data[neigh[k]],noise[neigh[k]],\
	kptr->data[neigh[l]],noise[neigh[l]],llav,dim);
	
	gg2[j][i]=gg2[i][j];
      }
    }
    
    after=double(time(NULL));
    covtime+=after-before;
    before=after;
    
    gg2det=fabs(invert_lapack(gg2,gg2in,n[1],1));    
    
    after=double(time(NULL));
    invtime+=after-before;
    before=after;
    
    //calculate sp1, the normalizing factor for gg1
    xx1=0.0;
    for(i=0;i<n[0];i++){
      xx1+=
      (fn[neigh[s[0]+i]]-zbar[0])*(fn[neigh[s[0]+i]]-zbar[0])*gg1in[i][i];
      
      for(j=i+1;j<n[0];j++){
        xx1+=2.0*(fn[neigh[s[0]+i]]-zbar[0])*(fn[neigh[s[0]+j]]-zbar[0])*\
        gg1in[i][j];
       }
    }
    
    sp1=xx1/double(n[0]);
    rawchi1=xx1;
    
    //calculate sp2, the normalizing factor for gg2
    xx2=0.0;
    for(i=0;i<n[1];i++){
      xx2+=(fn[neigh[s[1]+i]]-zbar[1])*(fn[neigh[s[1]+i]]-zbar[1])*\
      gg2in[i][i];
      for(j=i+1;j<n[1];j++){
        xx2+=2.0*(fn[neigh[s[1]+i]]-zbar[1])*(fn[neigh[s[1]+j]]-zbar[1])*\
        gg2in[i][j];
       }
    }
    rawchi2=xx2;
    
    sp2=xx2/double(n[1]);
    
    //calculate gg_1, the covariance matrix for the unimodal model
    for(i=0;i<kk;i++){
      gg_1[i][i]=1.0+lambda;
      for(j=i+1;j<kk;j++){
        if(i<n[0] && j<n[0]){
	  gg_1[i][j]=gg1[i][j];
	}
	else if(i>=n[0] && j>=n[0]){
	  gg_1[i][j]=gg2[i-n[0]][j-n[0]];
	}
	else{
	  gg_1[i][j]=covariogram(kptr->data[neigh[i]],noise[neigh[i]],\
      	  kptr->data[neigh[j]],noise[neigh[j]],llav,dim);
	}
     	gg_1[j][i]=gg_1[i][j];
	
      }
    }
    
    gg_1det=fabs(invert_lapack(gg_1,gg_1in,kk,1));
   
    
    //calculate the uncertainty^2 associated with the 1st mode of the bi-modal model
    sig2[0]=1.0+lambda;
    for(i=0;i<n[0];i++){
      sig2[0]-=ggq[0][i]*ggq[0][i]*gg1in[i][i];
      for(j=i+1;j<n[0];j++){
        sig2[0]-=2.0*ggq[0][i]*ggq[0][j]*gg1in[i][j];
      }
    }
    sig2[0]=sig2[0]*sp1;
    
    //calculate the uncertainty^2 associated with the 2nd mode of the bi-modal model
    sig2[1]=1.0+lambda;
    for(i=0;i<n[1];i++){
      sig2[1]-=ggq[1][i]*ggq[1][i]*gg2in[i][i];
      for(j=i+1;j<n[1];j++){
        sig2[1]-=2.0*ggq[1][i]*ggq[1][j]*gg2in[i][j];
      }
    }
    sig2[1]=sig2[1]*sp2;
   
    //this is a correction factor that takes into account the normalization of
    //the two modes as well as the degeneracy involved in dividing N_k galaxies
    //into two sub-populations of size n[0] and n[1]
    choosefactor
    =double(n[0])*log(double(n[0])/double(kk))
    +double(n[1])*log(double(n[1])/double(kk))
    +lnchoose(kk,n[0]);
    
    //this is the chi^2 associated with the likelihood of the
    //bi-modal model given the data (see equation 22 of our paper)
    chi2=double(kk)+
    double(n[0])*log(sp1)+double(n[1])*log(sp2)+
    log(gg1det)+log(gg2det)+
    double(kk)*log(2.0*pi)-2.0*choosefactor;
    
    //are we imposing a prior based on the separation of the two modes?   
    if(priorswitch==1){
      
      //deltab is the separation of the two modes
      deltab=fabs(zbar[1]-zbar[0]);
      
      //deltaz2 is the squared width of the prior
      //prior_width is a user-set normalization factor, if you want
      if(sp1>sp2){
        deltaz2=sp1*prior_width;
	ideltaz2=ik0;
      }
      else{
        deltaz2=sp2*prior_width;
	ideltaz2=ik1;
      }
      
      prior=double(kk)*(deltaz2/(deltab*deltab)-log(deltaz2))+double(kk)*log(2.0*pi);
      
      chi2+=prior;
    }
    
    prechi2=chi2;
        
    after=double(time(NULL));
    ittime+=after-before;
    before=after;
    
   //now we need to calculate the second derivative of chi2 with respect
   //to the model parameters (zbar and the normalization of the covariance matrix)
   //this will give us the Hessian we use for model comparison
   
   //the Hessian is stored in mm_det[]
   
   //I neglect minus signs in the calculation of dx1db,ddx1db, etc.
   //because they all end up cancelling when it comes time to assemble
   //the Hessian.
   
      dx1db=0.0;
      ddx1db=0.0;
      
      for(i=0;i<n[0];i++){
        for(j=0;j<n[0];j++){
	  dx1db+=gg1in[i][j]*(fn[neigh[j]]-zbar[0]);
	  ddx1db+=gg1in[i][j];
	}
      } 
    
    dx2db=0.0;
    ddx2db=0.0;
    
    for(i=s[1];i<kk;i++){
      k=i-s[1];
      for(j=s[1];j<kk;j++){
        l=j-s[1];
	dx2db+=gg2in[k][l]*(fn[neigh[j]]-zbar[1]);
	ddx2db+=gg2in[k][l];
      }
    }
    
    for(i=0;i<16;i++)mm_det[i]=0.0;
    
    mm_det[ib0*4+ib0]=ddx1db/sp1;
    
    mm_det[ib0*4+ik0]=dx1db/(sp1*sp1);
    mm_det[ik0*4+ib0]=mm_det[ib0*4+ik0];
    
    mm_det[ik0*4+ik0]=double(n[0])*0.5/(sp1*sp1);
    
    //////////////////
    
    mm_det[ib1*4+ib1]=ddx2db/sp2;
    
    mm_det[ib1*4+ik1]=dx2db/(sp2*sp2);
    mm_det[ik1*4+ib1]=mm_det[ib1*4+ik1];
    
    mm_det[ik1*4+ik1]=double(n[1])*0.5/(sp2*sp2);
    
    //include the parameters of the prior in the Hessian
    if(priorswitch==1){
      if(zbar[1]>zbar[0]){
        ddeltab1=1.0;
	ddeltab0=-1.0;
      }
      else{
        ddeltab1=-1.0;
	ddeltab0=1.0;
      }
      
      mm_det[ib0*4+ib0]+=3.0*double(kk)*deltaz2/power(deltab,4);
      mm_det[ib1*4+ib1]+=3.0*double(kk)*deltaz2/power(deltab,4);
      
      mm_det[ib0*4+ib1]+=(-3.0)*double(kk)*deltaz2/power(deltab,4);
      mm_det[ib1*4+ib0]+=(-3.0)*double(kk)*deltaz2/power(deltab,4);
      
      mm_det[ideltaz2*4+ideltaz2]+=0.5*double(kk)*prior_width*prior_width/(deltaz2*deltaz2);
      
      mm_det[ideltaz2*4+ib0]+=(-1.0)*double(kk)*ddeltab0*prior_width/(deltab*deltab*deltab);
      mm_det[ib0*4+ideltaz2]+=(-1.0)*double(kk)*ddeltab0*prior_width/(deltab*deltab*deltab);
      
      mm_det[ideltaz2*4+ib1]+=(-1.0)*double(kk)*ddeltab1*prior_width/(deltab*deltab*deltab);
      mm_det[ib1*4+ideltaz2]+=(-1.0)*double(kk)*ddeltab1*prior_width/(deltab*deltab*deltab);
      
      
      
      //remember: want -d^2lnP/dtheta dphi
      //but ddeltab1*ddeltab0 = -1
      //ddeltab_i^2=1
      
    }//if priorswitch==1
   
    
	for(i=0;i<4;i++){
	  rows[i]=1;
	  cols[i]=1;
	}
	
	det2=fabs(get_determinant(mm_det,rows,cols,4));
	chi2+=log(det2)-4.0*log(2.0*pi);
	
	
	after=double(time(NULL));
	dettime+=after-before;
	before=after;
  
  /////////////////////now to construct the single-mode model
      
    
    zbar_1=0.0;
    norm=0.0;
    
    for(i=0;i<kk;i++){
      zbar_1+=fn[neigh[i]];
      norm+=1.0;
    }
    zbar_1=zbar_1/norm;
   
    for(i=0;i<n[0];i++){
      ggq_1[i]=ggq[0][i];
    }
    for(;i<kk;i++){
      j=i-n[0];
      ggq_1[i]=ggq[1][j];
    }
    
    after=double(time(NULL));
    covtime+=after-before;
    before=after;
    
    xx1=0.0;
    for(i=0;i<kk;i++){
      xx1+=(fn[neigh[i]]-zbar_1)*(fn[neigh[i]]-zbar_1)*gg_1in[i][i];
      for(j=i+1;j<kk;j++){
        xx1+=2.0*(fn[neigh[i]]-zbar_1)*(fn[neigh[j]]-zbar_1)*\
	gg_1in[i][j];
      }
    }
    
    rawchi_1=xx1;
    sp_1=xx1/double(kk);
    

    //chi1 is the chi^2 associated with the likelihood of the single-mode
    //model given the data (again, see equation 22)
    chi1=double(kk)+double(kk)*log(sp_1)+log(gg_1det)+double(kk)*log(2.0*pi);
   
    //flat prior on the 1-mode model 
    if(priorswitch==1){
      chi1+=2.0*double(kk)*log(xpdf[npdf-1]);
    }
    
    prechi1=chi1;
    
    after=double(time(NULL));
    ittime+=after-before;
    before=after;
   
    
    /////////////calculate the Hessian for the 1-mode model
      
      dx1db=0.0;
      
      ddx1db=0.0;
  
      for(i=0;i<kk;i++){
        for(j=0;j<kk;j++){
	
	  dx1db+=gg_1in[i][j]*(fn[neigh[j]]-zbar_1);
	 
	  ddx1db+=gg_1in[i][j];
	
	}
      } 
      

      
	
	mm_det[ib0*2+ib0]=ddx1db/sp_1;
	
	mm_det[ib0*2+ik0]=dx1db/(sp_1*sp_1);
	mm_det[ik0*2+ib0]=mm_det[ib0*2+ik0];
	
	mm_det[ik0*2+ik0]=double(kk)*0.5/(sp_1*sp_1);
	
	
	for(i=0;i<2;i++){
	  rows[i]=1;
	  cols[i]=1;
	}
	
	det1=fabs(get_determinant(mm_det,rows,cols,2));
    
	chi1+=log(det1)-2.0*log(2.0*pi);
	
	if(det1==0.0 || isnan(det1)){
	  printf("WARNING chi1 det %e -- chi1 %e sp1 %e gg_1det %e\n",det1,chi1,sp1,gg_1det);
	  printf("kk %d n0 %d n1 %d\n",kk,n[0],n[1]);
	  printf("gg1det %e gg2det %e\n",gg1det,gg2det);
	
	  
	  for(i=0;i<2;i++){
	    for(j=0;j<2;j++)printf("%.3e ",mm_det[i*3+j]);
	    printf("\n");
	  }
	  printf("inversion err %e\n",check_inversion(gg_1,gg_1in,kk));
	  exit(1);
	}
	
	
	
       after=double(time(NULL));
       dettime+=after-before;
       before=after;
       

  
  //now calculate the uncertainty squared associated with the one-mode model
    
    sig2_1=1.0+lambda;
    for(i=0;i<kk;i++){
      sig2_1-=ggq_1[i]*ggq_1[i]*gg_1in[i][i];
      for(j=i+1;j<kk;j++){
        sig2_1-=2.0*ggq_1[i]*ggq_1[j]*gg_1in[i][j];
      }
    }
    sig2_1=sig2_1*sp_1;
  
  //if we are dividing chi1 and chi2 by the number of nearest neighbors
  //when choosing between models
  if(dofswitch==1){
   chi1=chi1/double(kk);
   chi2=chi2/double(kk);
  }
  
  
  //if twoisvalid==1, it means that the bi-modal model is
  //an acceptable alternative
  twoisvalid=1;
  if(!(det2>0.0))twoisvalid=0;
  if(!(sig2[0]>0.0))twoisvalid=0;
  if(!(sig2[1]>0.0))twoisvalid=0;
  if(isnan(det2))twoisvalid=0;
  if(isnan(sig2[0]))twoisvalid=0;
  if(isnan(sig2[1]))twoisvalid=0;
  if(isnan(chi2))twoisvalid=0;
  if(zbar[1]==zbar[0])twoisvalid=0;
  if(only1==1)twoisvalid=0;
  
  //if oneisvalid==1, it means that the unimodal model is
  //an acceptable alternative
  oneisvalid=1;
  if(!(det1>0.0))oneisvalid=0;
  if(isnan(det1))oneisvalid=0;
  if(!(sig2_1>0.0))oneisvalid=0;
  if(isnan(sig2_1))oneisvalid=0;
  if(isnan(chi1))oneisvalid=0;
 
  //find the minimum of the two chi^2 values for purposes
  //of doing the weighted sum of models
  if(oneisvalid==1)chimin=chi1;
  else if(twoisvalid==1)chimin=chi2;
  else chimin=0.0;
  
  if(oneisvalid==1 && twoisvalid==1){
    if(chi1<chi2)chimin=chi1;
    else chimin=chi2;
  }
  
  //count the number of times the bi-modal model
  //was a better fit
  if(chi2<chi1 && twoisvalid==1)chose2++;
  chi1out[0]=chi1;
  chi2out[0]=chi2;
  
  //build pdf_2 and add the bi-modal model to the
  //total pdf
  if(twoisvalid==1){
    
    mu[0]=zbar[0];
    for(i=0;i<n[0];i++){
      for(j=0;j<n[0];j++){
        mu[0]+=ggq[0][i]*gg1in[i][j]*(fn[neigh[s[0]+j]]-zbar[0]);
      } 
    }
   
    mu[1]=zbar[1];
    for(i=0;i<n[1];i++){
     for(j=0;j<n[1];j++){
       mu[1]+=ggq[1][i]*gg2in[i][j]*(fn[neigh[s[1]+j]]-zbar[1]);
     }
    }
    
  
    if(sig2[0]<=0.0 || sig2[1]<=0.0){
      printf("WARNING sig %e %e \n",sig2[0],sig2[1]);
      printf("n %d %d\n",n[0],n[1]);
      printf("mu %e %e\n",mu[0],mu[1]);
      printf("kk %d s %d %d\n",kk,s[0],s[1]);
      printf("sig1 %e\n",sig2_1);
      exit(1);
    }
    
    
    nn=chi2-chimin;
    for(i=0;i<npdf;i++){
      xx1=power(mu[0]-xpdf[i],2)/sig2[0];
      xx2=power(mu[1]-xpdf[i],2)/sig2[1];
      
      
      yy=(double(n[0])*exp(-0.5*xx1)/(double(kk)*sqrt(sig2[0]))+\
      double(n[1])*exp(-0.5*xx2)/(double(kk)*sqrt(sig2[1])));
      
      pdf[i]+=yy*exp(-0.5*nn);
      pdf_2[i]+=yy;
      
    
    }
  
  }//if the 2 GP is acceptable
  
  
  //build pdf_1 and add the unimodal model to the total pdf
  if(oneisvalid==1){

    mu_1=zbar_1;
    for(i=0;i<kk;i++){
      for(j=0;j<kk;j++){
        mu_1+=ggq_1[i]*gg_1in[i][j]*(fn[neigh[j]]-zbar_1);
      }
    }

    if(sig2_1<=0.0){
      printf("WARNING sig2_1 %e sp_1 %e\n",sig2_1,sp_1);
      exit(1);
    }
    
    nn=chi1-chimin;
    for(i=0;i<npdf;i++){
      xx1=(mu_1-xpdf[i])*(mu_1-xpdf[i])/sig2_1;
      yy=(exp(-0.5*xx1)/sqrt(sig2_1));
      
      pdf[i]+=yy*exp(-0.5*nn);
      pdf_1[i]+=yy;
      
    }
    
  
  }//if the 1 GP model was acceptable
  
  if(oneisvalid==1){
    mu_1out[0]=mu_1;
    sig_1out[0]=sig2_1;
  }
  else{
      mu_1out[0]=2.5;
      sig_1out[0]=1.0e20;
      chi1out[0]=1.0e20;
  }
  
  if(twoisvalid==1){
      mu_2aout[0]=mu[0];
      sig_2aout[0]=sig2[0];
      naout[0]=double(n[0]);
      mu_2bout[0]=mu[1];
      sig_2bout[0]=sig2[1];
      nbout[0]=double(n[1]);
      ntotout[0]=double(kk);
  }
  else{
      mu_2aout[0]=2.5;
      mu_2bout[0]=2.5;
      sig_2aout[0]=1.0e20;
      sig_2bout[0]=1.0e20;
      naout[0]=1.0;
      nbout[0]=1.0;
      ntotout[0]=2.0;
      chi2out[0]=1.0e20;
  
  }
  
  //if both models are invalid, return a flat
  //pdf
  if(oneisvalid==0 && twoisvalid==0){
  
    failures++;
    nn=0.0;
    for(i=0;i<npdf;i++){
      pdf[i]=1.0;
      nn+=pdf[i]*(xpdf[1]-xpdf[0]);
    }
  }
    
    //normalize the total pdf
    nn=0.0;
    for(i=0;i<npdf;i++){
      nn+=pdf[i]*(xpdf[1]-xpdf[0]);
    }
    if(nn==0.0){
        //printf("WARNING norm %e\n",norm);
	//for(i=0;i<npdf;i++)printf("%e %e %e\n",xpdf[i],pdf[i],pdf_1[i]);
	nn=1.0e-10;
    }
    for(i=0;i<npdf;i++){
      pdf[i]=pdf[i]/nn;
    }
    norm=nn;
    
    //normalize pdf_1
    nn=0.0;
    for(i=0;i<npdf;i++){
      nn+=pdf_1[i]*(xpdf[1]-xpdf[0]);
    }
    for(i=0;i<npdf;i++)pdf_1[i]=pdf_1[i]/nn;
    
    //normalize pdf_2
    nn=0.0;
    for(i=0;i<npdf;i++){
      nn+=pdf_2[i]*(xpdf[1]-xpdf[0]);
    }
    for(i=0;i<npdf;i++){
      pdf_2[i]=pdf_2[i]/nn;
    }
    
    after=double(time(NULL));
    ittime+=after-before;
    before=after;

    for(i=0;i<npdf;i++){
      if(isnan(pdf[i])){
      
      double mean1,mean[2];
      mean1=0.0;
      for(i=0;i<kk;i++)mean1+=fn[neigh[i]];
      mean1=mean1/double(kk);
      
      mean[0]=0.0;
      for(i=0;i<n[0];i++)mean[0]+=fn[neigh[i]];
      mean[0]=mean[0]/double(n[0]);
      
      mean[1]=0.0;
      for(i=s[1];i<s[1]+n[1];i++)mean[1]+=fn[neigh[i]];
      mean[1]=mean[1]/double(n[1]);
      
      printf("WARNING pdf isnan\n");
      printf("chimin %e \n",chimin);
      printf("chi2 %e chi 1 %e prior %e\n",chi2,chi1,prior);
      printf("prechi2 %e prechi1 %e\n",prechi2,prechi1);
      printf("dets %e %e %e\n",gg1det,gg2det,gg_1det);
      printf("sp %e %e %e\n",sp1,sp2,sp_1);
      printf("choose %e n %d %d\n",choosefactor,n[0],n[1]);
      printf("det %e %e\n",det1,det2);
      printf("choose test %e\n",lnchoose(kk,n[0]));
      printf("oneisvalid %d twoisvalid %d\n",oneisvalid,twoisvalid);
      printf("%e %e %e\n",sig2_1,sig2[0],sig2[1]);
      printf("%e %e %e\n",mu_1,mu[0],mu[1]);
      printf("%e %e %e\n",mean1,mean[0],mean[1]);
      printf("pdf1 %e norm %e\n",pdf_1[i],norm);
      exit(1);
      }
    }
    

  
 delete [] neigh;
 delete [] dd;
 delete [] ggq[0];
 delete [] ggq[1];
 delete [] ggq;
 for(i=0;i<kk;i++){
     delete [] gg1[i];
     delete [] gg2[i];
     delete [] gg_1[i];
     delete [] gg1in[i];
     delete [] gg2in[i];
     delete [] gg_1in[i];
 }
 delete [] gg1;
 delete [] gg2;
 delete [] gg_1;
 delete [] gg1in;
 delete [] gg2in;
 delete [] gg_1in;
 delete [] tosort;
 delete [] ggq_1;
 

}


void gpnoisy::get_time(int denom){
  printf("\n");
  printf("srchtime %e\n",srchtime/double(denom));
  printf("invtime %e\n",invtime/double(denom));
  printf("dettime %e\n",dettime/double(denom));
  printf("ittime %e\n",ittime/double(denom));
  printf("splittime %e\n",splittime/double(denom));
  printf("covartime %e\n",covtime/double(denom));
  printf("\n");
}


void gpnoisy::write_data(char *name){
  int i,j,k,l;
  
  FILE *output;
  output=fopen(name,"w");
  for(i=0;i<pts;i++){
    fprintf(output,"%e ",fn[i]);
    for(j=0;j<dim;j++)fprintf(output,"%e ",kptr->data[i][j]);
    fprintf(output,"\n");
  }
  fclose(output);
}

