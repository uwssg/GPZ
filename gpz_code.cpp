#include "goto_tools.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gaussian_process_driver.h"

double lsquared,ltotheone;

double covar(double *p1, double *n1, double *p2, double *n2, double llav,\
int dim){
  int i;
  double dd,core,ans,term,dx,ddx,denom,ll,num,aa,bb;

  ans=1.0;
  for(i=0;i<dim;i++){
    dd=power(p1[i]-p2[i],2);
    dd=dd/(n1[i]+n2[i]+lsquared*llav*llav);
    core=exp(-0.5*dd);
    
    denom=sqrt(n1[i]+n2[i]+lsquared*llav*llav);
    
    term=core*ltotheone*llav/denom;    
    ans=ans*term;
  }
  
  ///simple squared exponential covariogram
  /*
  dd=0.0;
  for(i=0;i<dim;i++){
      dd+=(p1[i]-p2[i])*(p1[i]-p2[i]);
  }
  ans=exp(-0.5*dd/lsquared);
  */
  
  ///Rasmussen and Wiliams (2006) eqn 4.29
  /*aa=1.0+lsquared;
  bb=1.0+lsquared;
  num=lsquared;
  for(i=0;i<dim;i++){
      aa+=p1[i]*p1[i]*lsquared;
      bb+=p2[i]*p2[i]*lsquared;
      num+=p1[i]*p2[i]*lsquared;
  }
  denom=sqrt(aa*bb);
  ans=asin(num/denom);*/
  
  if(ans<0.0){
    printf("WARNING covar is giving a neg %e\n",ans);
    for(i=0;i<dim;i++)printf("%e %e\n",p1[i],n1[i]);
    printf("\n");
    for(i=0;i<dim;i++)printf("%e %e\n",p2[i],n2[i]);
    printf("\n");
    printf("core %e\n",core);
    exit(1);
  }
  
  return ans;
}



main(int iargc, char *argv[]){

int i,j,k,l,ntest,ntrain,nvalid,dim,ncovar,ii,*bindex;
int kset,oswitch,jj,pdfsteps,normbyvar;
int eeset,c2best;

int kk,kkbest,kkmin,kkmax,kkstep,failures_best;

double floor,lnfloor;
double eebest,ee,biasbest;
double llin;

double ndll,chi1,chi2;

double **test,**tssig,**valid,**vsig,**train,**trsig;
double *ftest,*ftrain,*fvalid,*v,*min,*max,**dmu,nn,rr,mm;
double worst,bias,ll,llbest;

double before,after,pdftime,entropy;
double mu,sig,fbar,mode,pmode;

double *xpdf,*pdf,*pdf_1,*pdf_2,pdfsig;
double *mean,*var,dx;

double llmin,llmax,llstep;

double mu1,sig1,mu2a,sig2a,mu2b,sig2b,na,nb,ntot;

char testname[100],trainname[100],validname[100],keyword[100],paramname[100];
char outname[100],pdfname[100],statname[100];
int neighdex;
double ddneigh,zmax;
FILE *output,*input,*status;

gpnoisy gg;
oswitch=1;
normbyvar=0;

ndll=10.0;

///floor and lnfloor set the minimum possible
///value of ln[P(truth)] with optimizing
///hyper parameters
lnfloor=-60.0;
floor=exp(lnfloor);

dim=5;
kset=200;
eeset=0;



llmin=-1.0;
llmax=10.0;

pdfsteps=500;
kkmin=50;
kkmax=150;
kkstep=10;
zmax=5.0;

for(j=0;argv[1][j]!=0;j++)paramname[j]=argv[1][j];
paramname[j]=0;

printf("paramname %s\n",paramname);

input=fopen(paramname,"r");
while(fscanf(input,"%s",keyword)>0){
  if(compare_char(keyword,"#dim")==1){
    fscanf(input,"%d",&dim);
    
  }
  else if(compare_char(keyword,"#train")==1){
    fscanf(input,"%s",trainname);
   
  }
  else if(compare_char(keyword,"#test")==1){
    fscanf(input,"%s",testname);
   
  }
  else if(compare_char(keyword,"#valid")==1){
    fscanf(input,"%s",validname);
    
  }
  else if(compare_char(keyword,"#dof")==1){
    gg.dofswitch=1;
  }
  else if(compare_char(keyword,"#prior")==1){
    gg.priorswitch=1;
  }
  else if(compare_char(keyword,"#output")==1){
    fscanf(input,"%s",outname);
  }
  else if(compare_char(keyword,"#status")==1){
    fscanf(input,"%s",statname);
  }
  else if(compare_char(keyword,"#hyperparams")==1){
   oswitch=0;
   fscanf(input,"%le %d",&llin,&kset);
   printf("kset %d\n",kset);
 }
 else if(compare_char(keyword,"#llrange")==1){
   fscanf(input,"%le %le %le",&llmin,&llmax,&ndll);
 }
 else if(compare_char(keyword,"#kkrange")==1){
   fscanf(input,"%d %d %d",&kkmin,&kkmax,&kkstep);
 }
 else if(compare_char(keyword,"#normbyvar")==1){
   normbyvar=1;
 }
 else if(compare_char(keyword,"#only1")==1){
     gg.only1=1;
 }
 else if(compare_char(keyword,"#prior_width")==1){
     fscanf(input,"%le",&gg.prior_width);
 }
 else if(compare_char(keyword,"#zmax")==1){
     fscanf(input,"%le",&zmax);
 }

}
fclose(input);

dx=zmax/double(pdfsteps);

if(llmax<llmin){
  nn=llmax;
  llmax=llmin;
  llmin=nn;
}
llstep=(llmax-llmin)/ndll;

if(kkmax<kkmin){
    i=kkmax;
    kkmax=kkmin;
    kkmin=i;
}

////figure out how many data points are in the
////test, training, and validation data sets

input=fopen(trainname,"r");
for(ntrain=0;fscanf(input,"%le",&ee)>0;ntrain++){
  for(i=1;i<2*(dim)+1;i++)fscanf(input,"%e",&ee);
}
fclose(input);


input=fopen(testname,"r");
for(ntest=0;fscanf(input,"%le",&ee)>0;ntest++){
  for(i=1;i<2*(dim)+1;i++)fscanf(input,"%e",&ee);
}
fclose(input);


input=fopen(validname,"r");
for(nvalid=0;fscanf(input,"%le",&ee)>0;nvalid++){
  for(i=1;i<2*(dim)+1;i++)fscanf(input,"%e",&ee);
}
fclose(input);


/////read in the test, training, and validation sets

v=new double[dim*2+1];
test=new double*[ntest];
tssig=new double*[ntest];
ftest=new double[ntest];
input=fopen(testname,"r");
i=0;
while(fscanf(input,"%le",&v[0])>0){
  
  for(j=1;j<2*(dim)+1;j++)fscanf(input,"%le",&v[j]);

    test[i]=new double[dim];
    tssig[i]=new double[dim];
    for(j=0;j<dim;j++)test[i][j]=v[j];
    for(j=0;j<dim;j++)tssig[i][j]=v[dim+j];
    ftest[i]=v[2*dim];

    i++;
  
  
}
fclose(input);
if(i!=ntest){
    printf("WARNING miscounted ntest\n");
    exit(1);
}


train=new double*[ntrain];
trsig=new double*[ntrain];
ftrain=new double[ntrain];
input=fopen(trainname,"r");
i=0;
while(fscanf(input,"%le",&v[0])>0){
  
  for(j=1;j<2*(dim)+1;j++)fscanf(input,"%le",&v[j]);

    train[i]=new double[dim];
    trsig[i]=new double[dim];
    for(j=0;j<dim;j++)train[i][j]=v[j];
    for(j=0;j<dim;j++)trsig[i][j]=v[dim+j];
    ftrain[i]=v[2*dim];
   
    i++;
   
  
}
fclose(input);
if(i!=ntrain){
    printf("WARNING miscounted ntrain\n");
    exit(1);
}

valid=new double*[nvalid];
vsig=new double*[nvalid];
fvalid=new double[nvalid];
input=fopen(validname,"r");
i=0;
while(fscanf(input,"%le",&v[0])>0){
  
  for(j=1;j<2*(dim)+1;j++)fscanf(input,"%le",&v[j]);

    valid[i]=new double[dim];
    vsig[i]=new double[dim];
    for(j=0;j<dim;j++)valid[i][j]=v[j];
    for(j=0;j<dim;j++)vsig[i][j]=v[j+dim];
    fvalid[i]=v[2*dim];
    i++;
  
  
}
fclose(input);
if(i!=nvalid){
    printf("WARNING miscounted nvalid\n");
    exit(1);
}

status=fopen(statname,"w");
fprintf(status,"finally ntrain %d ntest %d nvalid %d\n",ntrain,ntest,nvalid);
fclose(status);


///option to normalize the flux and noise data so that
///each filter is on equal footing
if(normbyvar==1){
  status=fopen(statname,"a");
  fprintf(status,"norming by variance\n");
  fclose(status);
  
  mean=new double[dim+1];
  var=new double[dim+1];
  for(i=0;i<dim+1;i++){
    mean[i]=0.0;
    var[i]=0.0;
  }
  
  for(i=0;i<ntrain;i++){
    for(j=0;j<dim;j++){
       mean[j]+=train[i][j];
       var[j]+=train[i][j]*train[i][j];
    }
    mean[dim]+=ftrain[i];
    var[dim]+=ftrain[i]*ftrain[i];
  }
  
  status=fopen(statname,"a");
  for(i=0;i<dim+1;i++){
    mean[i]=mean[i]/double(ntrain);
    var[i]=var[i]/double(ntrain-1)-\
    double(ntrain)*mean[i]*mean[i]/double(ntrain-1);
    
    printf("i %d var %e\n",i,var[i]);
    
    var[i]=sqrt(var[i]);
    fprintf(status,"   mean%d %e pm %e\n",i,mean[i],var[i]);
    
  }
  fclose(status);
   
  for(i=0;i<ntrain;i++){
    for(j=1;j<dim;j++){
      train[i][j]=train[i][j]*var[0]/var[j];
      trsig[i][j]=trsig[i][j]*var[0]/var[j];
    }
  }
  
  for(i=0;i<nvalid;i++){
    for(j=1;j<dim;j++){
      valid[i][j]=valid[i][j]*var[0]/var[j];
      vsig[i][j]=vsig[i][j]*var[0]/var[j];
    }
  }
  
  for(i=0;i<ntest;i++){
    for(j=1;j<dim;j++){
      test[i][j]=test[i][j]*var[0]/var[j];
      tssig[i][j]=tssig[i][j]*var[0]/var[j];
    }
  }
  
  delete [] mean;
  delete [] var;
}  

///convert the noise arrays to noise^2 arrays so that
///we do not have to waste time squaring the noise in the
///covariogram every time we call it
for(i=0;i<ntrain;i++){
  for(j=0;j<dim;j++)trsig[i][j]=trsig[i][j]*trsig[i][j];
}
for(i=0;i<ntest;i++){
  for(j=0;j<dim;j++)tssig[i][j]=tssig[i][j]*tssig[i][j];
}
for(i=0;i<nvalid;i++){
  for(j=0;j<dim;j++)vsig[i][j]=vsig[i][j]*vsig[i][j];
}

gg.covariogram=covar;
gg.dim=dim;
gg.kk=kset;
gg.initialize(ntrain,train,ftrain,trsig);

for(i=0;statname[i]!=0;i++){
    gg.statname[i]=statname[i];
}
gg.statname[i]=0;

//this is the "nugget" in the Gaussian Process' covariogram;
//we are just setting it to a small number here;
//it helps ensure that the covariogram is non-singular
gg.set_lambda(1.0e-6);


bias=0.0;
ee=0.0;
before=double(time(NULL));
eebest=-1.0;
  bindex=new int[nvalid];
  xpdf=new double[pdfsteps];
  for(i=0;i<pdfsteps;i++){
    xpdf[i]=(double(i)+0.5)*dx;
  }
  
  ///find the index in xpdf[] of the redshift value
  ///of each validation galaxy
  for(i=0;i<nvalid;i++){
    for(j=0;j<pdfsteps-1 && xpdf[j]<fvalid[i];j++);
    if(j>0 && fvalid[i]-xpdf[j-1]<xpdf[j]-fvalid[i])j--;
    bindex[i]=j;
  }
  

pdf=new double[pdfsteps];
pdf_1=new double[pdfsteps];
pdf_2=new double[pdfsteps];

///add the validation galaxies to the training set, so that
///they can be used for infering the pdfs
for(i=0;i<nvalid;i++){
    gg.add_pt(valid[i],fvalid[i],vsig[i]);
}


///option to optimize ell and N_k
lsquared=1.0;
ltotheone=1.0;
if(oswitch==1){
 k=0;
 for(kk=kkmin;kk<=kkmax;kk+=kkstep){
   gg.kk=kk;
   for(ll=llmin;ll<=llmax;ll+=llstep){   
     lsquared=exp(ll*log(10.0));
     ltotheone=sqrt(lsquared);
   
     gg.chose2=0;
     gg.failures=0;
     ee=0.0;
     bias=0.0;
   
     for(i=0;i<nvalid;i++){
     
       ///get the pdf for each validation galaxy;
       ///the last argument of get_pdf tells the code not to use
       ///the galaxy in question as a nearest neighbor of itself
     
       gg.get_pdf(valid[i],vsig[i],xpdf,pdf,pdf_1,pdf_2,pdfsteps,&nn,&nn,
       &nn,&nn,&nn,&nn,&nn,&nn,&nn,&nn,&nn,ntrain+i);
    
       //optimize on the sum of ln[P(truth)]
       //with a minimum possible value set by lnfloor
       if(pdf[bindex[i]]>floor)nn=log(pdf[bindex[i]]);
       else nn=lnfloor;    
       ee+=nn;

       if(isnan(ee)){
          printf("ee %e kk %d ll %e nn %e\n",ee,kk,ll,nn);
          exit(1); 
       }
    }//loop over validation galaxies
  
    if(eeset==0 || ee>eebest){
     eebest=ee;
     llbest=ll;
     kkbest=kk;
     failures_best=gg.failures;
     c2best=gg.chose2;
     
     status=fopen(statname,"a");
     fprintf(status,
     "     eebest %e -- log10(ll) %e kk %d failures %d chose2 %d\n",\
     eebest,llbest,kkbest,failures_best,c2best);
     fclose(status);
     eeset=1;
     
   }
   
   k++;
   if(k%10==0 && k!=0){
     status=fopen(statname,"a");
     after=double(time(NULL));
     fprintf(status,\
     "per optimization step %e eebest %e --  log10(ll) %e kk %d bs %d\n",\
     (after-before)/double(k),eebest,llbest,kkbest,failures_best);
     
     printf(\
     "per optimization step %e eebest %e --  log10(ll) %e kkbest %d bs %d \n",\
     (after-before)/double(k),eebest,\
     llbest,kkbest,failures_best);     
     fclose(status);
   }
 
 
  }//loop over ll
  }//loop over kk

  lsquared=exp(llbest*log(10.0));
  ltotheone=sqrt(lsquared);
  gg.kk=kkbest;

  status=fopen(statname,"a");
  after=double(time(NULL));
  fprintf(status,\
        "total optimization time %e per galaxy %e\n",
	after-before,(after-before)/double(k));
        fprintf(status,"log10(ll) %e kk %d bs %d ee %e \n",\
        llbest,kkbest,failures_best,eebest);
  fclose(status);
}//if oswitch==1
else{
 //or you can just assert the values of ll and N_k

 lsquared=exp(llin*log(10.0));
 ltotheone=sqrt(lsquared);
 gg.kk=kset;
}
////////////////

bias=0.0;
ee=0.0;
gg.chose2=0;
gg.failures=0;

k=0;
pdftime=0.0;
before=double(time(NULL));
output=fopen(outname,"w");
fprintf(output,"# chi1 mu1 sig1 chi2 mu2a sig2a na mu2b sig2 nb ntot");
fprintf(output,"\n");
fclose(output);

printf("lsquared %e kk %d\n",lsquared,gg.kk);
status=fopen(statname,"a");
fprintf(status,"lsquared %e kk %d\n\n",lsquared,gg.kk);
fclose(status);

for(i=0;i<ntest;i++){
 
 gg.get_pdf(test[i],tssig[i],xpdf,pdf,pdf_1,pdf_2,pdfsteps,
            &chi1,&mu1,&sig1,
            &chi2,&mu2a,&sig2a,&na,
            &mu2b,&sig2b,&nb,&ntot,-1);


 output=fopen(outname,"a");
 //fprintf(output,"# chi1 mu1 sig1 chi2 mu2a sig2a na mu2b sig2 nb ntot");
 fprintf(output,"%e %e %e %e %e %e %e %e %e %e %e\n",
 chi1,mu1,sig1,
 chi2,mu2a,sig2a,na,
 mu2b,sig2b,nb,ntot);
 fclose(output);
 
 if(i%1000==0){
   output=fopen(statname,"a");
   after=double(time(NULL));
   fprintf(output,"got %d in %e av %e chose2 %d\n",
   i,after-before,(after-before)/double(i+1),gg.chose2);
   fclose(output);
   
   printf("got %d\n",i);
   
 }
 

}
after=double(time(NULL));

status=fopen(statname,"a");
fprintf(status,"bias %e ee %e chose2 %d\n",bias/double(ntest),ee/double(ntest),gg.chose2);
fclose(status);

}
