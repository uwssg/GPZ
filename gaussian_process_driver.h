#include "kd.h"

class gpnoisy{
  
  private:
    
    double *fn,**noise,lambda;
    int initialized,room,roomstep;

  public:
    kd_tree *kptr;
    int dim,kk,pts,failures,only1;
    int dofswitch,priorswitch,chose2;
    double srchtime,invtime,ittime,splittime;
    double dettime,covtime,prior_width;
    double dav,dsig;
    double maxgg,maxggin;
    char statname[100];
    
    gpnoisy();
    ~gpnoisy();
   
    double(*covariogram)(double*,double*,double*,double*,double,int);
   
    void initialize(int,double**,double*,double**);
   
    void get_pdf(double*,double*,double*,double*,
                      double*,double*,int,
		      double*,double*,double*,double*,double*,
		      double*,double*,double*,
		      double*,double*,double*,int);
		      
    void add_pt(double*,double,double*);
   
    void write_data(char*);
   
    void set_lambda(double);
   
    double p_lambda();

    void get_time(int);

  
};
