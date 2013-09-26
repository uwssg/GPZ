//this is the 27 February 2012 rewrite

class kd_tree{
 private:
  int **tree;
  int masterparent;
  int room,roomstep;
  double tol;
  int ncube,callcubes;
  
  int nkernel;
  
  void confirm(int,int,int,int);
  void organize(int*,int,int,int,int);
  int find_node(double*);
  void neigh_check(double*,int,int*,double*,int,int);
  void kernel_check(double*,double*,int*,int,int);
   void reassign(int);
   void descend(int);
  
 public:
  int dim,pts,diagnostic,xplr,ktests;
  
  int **cdex;
  double **cubecenter,*cubevol,**cubemax,**cubemin;
  int **cube_mindex,**cube_maxdex;
  
  double **data,*maxs,*mins;

  kd_tree(int,int,double**,double*,double*);
  ~kd_tree();
 
 void check_tree(int);
 double distance(double*,double*);
 void find_cubes();

 void write_tree(char*);
 void add(double*);
 void remove(int);
 void count(int,int*);
 void nn_srch(double*,int,int*,double*);
 int kernel_srch(double*,double*,int*);

};

class kd_cube{
  
  private:
    int dim,pts,*inn,*icalled,ocalled,ocalled1,nhack,zeroct;
    int *priority;
    double **minbound,**maxbound,*tosort;
    double avict;
    double *mins,*maxs;
    
    double *center,*range;
    
    void organize(int,int,double**,int,double*);
    
  public:
     int varswitch;
  
    ~kd_cube();
    kd_cube(int,int,double**,double*,double*,double*,int,int*);
    double pmin(int,int);
    double pmax(int,int);
    void picalled();
};
