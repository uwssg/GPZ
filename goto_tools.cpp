#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "goto_tools.h"


int compare_char(char *s1, char *s2){
 //are two character strings the same?
 //if so, return 1
 //if not, return 0

 int i;
 //printf("\ncomparing %s %s\n",s1,s2);
 for(i=0;i<letters && (s1[i]!=0 || s2[i]!=0);i++){
  if(s1[i]!=s2[i])return 0;
 }
 return 1;
 
}

double power(double arg,int raised){

  //return arg raised to the integer power 'raised'

  int n;
  double ans;

  if(raised==0)return 1.0;
  else{ans=1.0;
  for(n=0;n<raised;n++){
    ans=ans*arg;
  }

  return ans;
  }

}

void kill(char *words){
 double junk;
 FILE *output;
 
 //write the character string words[] to the terminal
 //then hang, waiting for an input double
 
 printf("%s\n",words);
 scanf("%lf",&junk);
}

/*just going to lift this out of Numerical Recipes p 109*/
void polint(double *xa, double *ya, int n, double x, double *y, double *dy){
/*recall: because this is from NR, it has arrays start at element unity*/
	int i,m,ns=1,isfd;
	double den,dif,dift,ho,hp,w;
	double c[n+1],d[n+1];
	char scream[100];
	dif=fabs(x-xa[1]);
	for(i=1;i<=n;i++){
		if((dift=fabs(x-xa[i]))<dif){
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for(m=1;m<n;m++){
		for(i=1;i<=n-m;i++){
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if((den=ho-hp)==0.0){printf("Error in routine polint ");
				for(isfd=1;isfd<=n;isfd++)printf(" (%e, %e) ",xa[isfd],ya[isfd]);
				printf(" want %e \n",x); 
				sprintf(scream,"stop");
				kill(scream);
			
			
			}
			/*This error can occur only if two input xas re (to within roundoff) identical*/
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y+=(*dy=(2*ns<(n-m)?c[ns+1]:d[ns--]));
	}
}


double interpolate(double *x, double *y, double target, int el){

  //this is a wrapper for the Numerical Recipes polynomial interpolation
  //routine polint()
  
  //x[] is the array of x values
  
  //y[] is the array of y values
  
  //'target' is the x value for which you want to interpolate a y value
  
  //'el' is the number of elements in x[] and y[]
  
  //x[] must be monotonic, but it can be either ascending or descending

 double *xt,*yt;
 double xint[5],yint[5],err,ans;
 int i,n,min;
 
 
 if(x[0]>x[1]){
  xt=new double[el];
  yt=new double[el];
  for(i=0;i<el;i++){
   xt[i]=x[el-1-i];
   yt[i]=y[el-1-i];
  
  }
 }
 else{
  xt=x;
  yt=y;
 }

 
 for(i=0;xt[i]<target && i<el;i++);
 if(xt[i]==target){ans=yt[i];}
 else{
  if(i<2)min=0;
  else if(i>el-2)min=el-4;
  else min=i-2;
  for(n=1;n<5;n++){
   xint[n]=xt[min+n-1];
   yint[n]=yt[min+n-1];
  }
  //printf("min was %d %e i %d %e target %e lim %e %e\n",min,xt[min],i,xt[i],target,xt[0],xt[el-1]);
  polint(xint,yint,4,target,&ans,&err);
 }
 
 if(x[0]>x[1]){
  delete [] xt;
  delete [] yt;
 }
 
 //if(!(ans>0.0))scanf("%lf",&junk);
 return ans;

}


//the routines below are just the merge sort algorithm from Numerical Recipes;
//(2nd edition); Press, Teukolsky, Vetterling, Flannery 1992

int scanner(double *m, int *indices, int index, int el){
/*this will take the matrix m and put everything in it with value
greater than element m[index] to the right of that element
and everything less than m[index] to the left; it then returns
the new index of that anchored element (which you now *know* is in
the right spot*/

	double toobig,toosmall;
	int i,j,n,k,pp;
	int itoobig,itoosmall;
        FILE *output;

	itoobig=0;
	itoosmall=0;
	for(i=0;i<el;i++){
	  if(m[i]<=m[index])itoosmall++;
	  
	} 
	
	toobig=m[index];
	itoobig=indices[index];
	
	indices[index]=indices[itoosmall-1];
	m[index]=m[itoosmall-1];
	
	index=itoosmall-1;
	m[index]=toobig;
	indices[index]=itoobig;

	i=0;
	j=el-1;
	n=0;
	pp=1;
	while(index<el-1 && n<el-1){
	 for(;m[i]<=m[index] && i<index;i++,n++);
	 itoobig=indices[i];
	 toobig=m[i];
	 
	 for(;m[j]>m[index] && j>index;j--,n++){
	  
	 }
	 itoosmall=indices[j];
	 toosmall=m[j];
	
	 
	 
	 if(toosmall<toobig){
	
	 //in case it ran to the end and i found m[index]
	  m[i]=toosmall;
	  indices[i]=itoosmall;
	 
	  m[j]=toobig;
	  indices[j]=itoobig;
	 }
	 
	}

	return index;
}

void sort(double *m, int *indices, int el){
	double *newm,junk;
	int k,i,*newi,ii,diff;
	
	if(el==2){
	 if(m[0]>m[1]){
	  
	   junk=m[1];
	   m[1]=m[0];
	   m[0]=junk;
	   
	   k=indices[1];
	   indices[1]=indices[0];
	   indices[0]=k;
	  
	 }
	}
	else if(el>2){
	
	diff=0; //just in case all elements are identical
	for(ii=1;ii<el;ii++)if(m[ii]!=m[0])diff=1;
	  
	if(diff==1){

	   i=scanner(m,indices,el/2,el);

	   if(i+1<el){
	     newm=&m[i+1];
	    sort(m,indices,i+1);
	     newi=&indices[i+1];
	    sort(newm,newi,el-i-1);
	  }
	  else{
	    sort(m,indices,el-1);
	  }
	
	 }
	
	}
}

void check_sort(double *sorted, double *unsorted, int *inn, int *inn_0, int n){
    int i,j;
    double diff;
    
    for(i=0;i<n;i++){
        for(j=0;j<n && inn_0[j]!=inn[i];j++);
	
        diff=fabs(sorted[i]-unsorted[j]);
	if(unsorted[inn[i]]!=0.0){
	    diff=diff/fabs(unsorted[inn[i]]);
	}
	if(diff>1.0e-10){
	    printf("WARNING sort failed to associate\n");
	    exit(1);
	}
    }
    
    for(i=1;i<n;i++){
        if(sorted[i]<sorted[i-1]){
	    printf("WARNING sort failed to sort\n");
	    exit(1);
	}
    }
    
    
}


void sort_and_check(double *list, double *sorted, int *inn, int n){
    int i,*inn_0;
   
    inn_0=new int[n];
    for(i=0;i<n;i++){
        sorted[i]=list[i];
	inn_0[i]=inn[i];
    }
    sort(sorted,inn,n);
   // printf("sorted\n");
    check_sort(sorted,list,inn,inn_0,n);
   // printf("checked\n");
    
    delete [] inn_0;
}

