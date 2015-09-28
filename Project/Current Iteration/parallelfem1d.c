# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include<sys/time.h>
# include <omp.h>

#define NUM_THREADS 4

int main ( void );
void assemble ( double adiag[], double aleft[], double arite[], double f[], 
  double h[], int indx[], int nl, int node[], int nu, int nquad, int nsub, 
  double ul, double ur, double xn[], double xquad[] );
double ff ( double x );
void geometry ( double h[], int ibc, int indx[], int nl, int node[], int nsub, 
  int *nu, double xl, double xn[], double xquad[], double xr );
void init ( int *ibc, int *nquad, double *ul, double *ur, double *xl, 
  double *xr );
void output ( double f[], int ibc, int indx[], int nsub, int nu, double ul, 
  double ur, double xn[] );
void phi ( int il, double x, double *phii, double *phiix, double xleft, 
  double xrite );
double pp ( double x );
void prsys ( double adiag[], double aleft[], double arite[], double f[], 
  int nu );
double qq ( double x );
void solve ( double adiag[], double aleft[], double arite[], double f[], 
  int nu );
void timestamp ( void );

/******************************************************************************/

int main ( void )
{
# define NSUB 80000
# define NL 4



  //double adiag[NSUB+1];
  double *adiag;
  adiag = (double *)malloc(sizeof(double)*(NSUB+1));
  //double aleft[NSUB+1];
  double *aleft;
  aleft = (double *)malloc(sizeof(double)*(NSUB+1));
  //double arite[NSUB+1];
  double *arite;
  arite = (double *)malloc(sizeof(double)*(NSUB+1));

  double *f;
  f = (double *)malloc(sizeof(double)*(NSUB+1));
 // double h[NSUB];
  double *h;
  h = (double *)malloc(sizeof(double)*(NSUB));
  int ibc;
  //int indx[NSUB+1];
  int *indx;
  indx = (int *)malloc(sizeof(double)*(NSUB+1));
  //int node[NL*NSUB];
  int *node;
  node = (int*)malloc(sizeof(double)*(NL*NSUB));
  int nquad;
  int nu;
  double ul;
  double ur;
  double xl;
  //double xn[NSUB+1];
  double *xn;
  xn = (double *)malloc(sizeof(double)*(NSUB+1));
  //double xquad[NSUB];
  double *xquad;
  xquad = (double *)malloc(sizeof(double)*(NSUB));
  double xr;

 
  timestamp ( );


  FILE *fp, *fopen();
  fp =fopen("paralleloutput.txt","a");
  fprintf (fp, "\n" );
  fprintf (fp, "FEM1D\n" );
  fprintf (fp, "  C version\n" );
  fprintf (fp, "\n" );
  fprintf (fp, "  Solve the two-point boundary value problem\n" );
  fprintf (fp, "\n" );
  fprintf (fp, "  - d/dX (P dU/dX) + Q U  =  F\n" );
  fprintf (fp, "\n" );
  fprintf (fp, "  on the interval [XL,XR], specifying\n" );
  fprintf (fp, "  the value of U or U' at each end.\n" );
  fprintf (fp, "\n" );
  fprintf (fp, "  The interval [XL,XR] is broken into NSUB = %d subintervals\n", NSUB );
  fprintf (fp, "  Number of basis functions per element is NL = %d\n", NL );

  fclose(fp);


/***********************************************
  Start timing 
************************************************/

struct timeval start, end;
gettimeofday(&start, NULL);

/***********************************************

************************************************/

/*
  Initialize the data.

*/
  init ( &ibc, &nquad, &ul, &ur, &xl, &xr );
/*
  Compute the geometric quantities.
*/
  geometry ( h, ibc, indx, NL, node, NSUB, &nu, xl, xn, xquad, xr );
/*
  Assemble the linear system.
*/
  assemble ( adiag, aleft, arite, f, h, indx, NL, node, nu, nquad, 
    NSUB, ul, ur, xn, xquad );
/*
  Print out the linear system.
*/
  prsys ( adiag, aleft, arite, f, nu );
/*
  Solve the linear system.
*/
  solve ( adiag, aleft, arite, f, nu );
/*
  Print out the solution.
*/
  output ( f, ibc, indx, NSUB, nu, ul, ur, xn );
/*
  Terminate.
*/

  gettimeofday(&end, NULL);
  double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
  
  
  printf("Parallel Time = %12.10f seconds\n",delta);

  FILE *dp, *fopen();
  dp =fopen("paralleloutput.txt","a");
  fprintf (dp, "\n" );
  fprintf (dp, "FEM1D:\n" );
  fprintf (dp, "  Normal end of execution.\n\n" );
  //timestamp ( );
  fprintf(dp,"Parallel Time = %12.10f Seconds\n",delta);
 printf ( "*****************************************\n" );
  fclose(dp);
  
  printf ( "\n" );

  

  return 0;
# undef NL
# undef NSUB
}
/******************************************************************************/

void init ( int *ibc, int *nquad, double *ul, double *ur, double *xl, 
  double *xr )
{
/*
  IBC declares what the boundary conditions are.
*/
  *ibc = 1;
/*
  NQUAD is the number of quadrature points per subinterval.
  The program as currently written cannot handle any value for
  NQUAD except 1.
*/
  *nquad = 1;
/*
  Set the values of U or U' at the endpoints.
*/
  *ul = 0.0;
  *ur = 1.0;
/*
  Define the location of the endpoints of the interval.
*/
  *xl = 0.0;
  *xr = 1.0;
/*
  Print out the values that have been set.
*/
  FILE *fp, *fopen();
  fp =fopen("paralleloutput.txt","a");
  
  fprintf (fp, "\n" );
  fprintf (fp, "  The equation is to be solved for\n" );
  fprintf (fp, "  X greater than XL = %f\n", *xl );
  fprintf (fp, "  and less than XR = %f\n", *xr );
  fprintf (fp, "\n" );
  fprintf (fp, "  The boundary conditions are:\n" );
  fprintf (fp, "\n" );

  if ( *ibc == 1 || *ibc == 3 )
  {
    fprintf (fp, "  At X = XL, U = %f\n", *ul );
  }
  else
  {
    fprintf (fp, "  At X = XL, U' = %f\n", *ul );
  }

  if ( *ibc == 2 || *ibc == 3 )
  {
    fprintf (fp, "  At X = XR, U = %f\n", *ur );
  }
  else
  {
    fprintf (fp, "  At X = XR, U' = %f\n", *ur );
  }

  fprintf (fp, "\n" );
  fprintf (fp, "  Number of quadrature points per element is %d\n", *nquad );

  fclose(fp);
  
  return;
}
/******************************************************************************/

void geometry ( double h[], int ibc, int indx[], int nl, int node[], int nsub, 
  int *nu, double xl, double xn[], double xquad[], double xr )

{

  double gstart = omp_get_wtime();

  //Open text file for printing

  FILE *fp, *fopen();
  fp =fopen("paralleloutput.txt","a");
  

  int i;
/*
  Set the value of XN, the locations of the nodes.
*/
  fprintf (fp, "\n" );
  fprintf (fp, "  Node      Location\n" );
  fprintf (fp, "\n" );

//Loop G1




omp_set_num_threads(4);
#pragma omp parallel 
  {
    #pragma omp for
    for ( i = 0; i <= nsub; i++ )
    {
      xn[i]  =  ( ( double ) ( nsub - i ) * xl 
                + ( double )          i   * xr ) 
                / ( double ) ( nsub );
    }

//Print statement withing Loop G1 caused excessive slow down when parallelising
//hence, it has been taken out of loop
    #pragma omp single    
      for ( i = 0; i <= nsub; i++ )
        {
          fprintf (fp, "  %8d  %14f \n", i, xn[i] );
        }
      
  }

fprintf (fp, "\n" );
fprintf (fp, "Subint    Length\n" );
fprintf (fp, "\n" );
 

omp_set_num_threads(4);
#pragma omp parallel 
{
//Loop G2
  #pragma omp for 
    for ( int i = 0; i < nsub; i++ )
    {
      h[i] = xn[i+1] - xn[i];
      //fprintf (fp, "  %8d  %14f\n", i+1, h[i] );
    }
  #pragma omp single
    {
      for ( i = 0; i < nsub; i++ )
      {
        fprintf (fp, "  %8d  %14f\n", i+1, h[i] );
      }
    }
  }


/*  
  Set the quadrature points, each of which is the midpoint
  of its subinterval.
*/

      fprintf (fp, "\n" );
      fprintf (fp, "Subint    Quadrature point\n" );
      fprintf (fp, "\n" );

 /*
  Similar to the previous parallel block Loops G3 & G4 intially contained 
  print statements within the loops. While the mundane computation and
  filling on arrays was a task well suited to parallelisation, the printing
  of the statements caused slowdown, hence the printing operation was placed
  into a separate loop. The implicit barrier at the end of the for loop
  may cause some parallel overhead but subsequent tests found that this
  section had an average time of 0.0867s compared to 0.1025s with serial
  execution.
 */ 


#pragma omp parallel num_threads(NUM_THREADS)
{
 //Loop G3
  #pragma omp for
      for ( i = 0; i < nsub; i++ )
      {
        xquad[i] = 0.5 * ( xn[i] + xn[i+1] );
       
      }
  #pragma omp single nowait
      {
        for ( i = 0; i < nsub; i++ )
        {
          fprintf (fp, "  %8d  %14f\n", i+1, xquad[i] );
        }
      /*
        Set the value of NODE, which records, for each interval,
        the node numbers at the left and right.
      */
        fprintf (fp, "\n" );
        fprintf (fp, "Subint  Left Node  Right Node\n" );
        fprintf (fp, "\n" );
      }
     //Loop G4
      #pragma omp for
        for ( i = 0; i < nsub; i++ )
        {
          node[0+i*2] = i;
          node[1+i*2] = i + 1;   
        }
    #pragma omp single 
      {
        for ( i = 0; i < nsub; i++ )
        {
          fprintf (fp, "  %8d  %8d  %8d\n", i+1, node[0+i*2], node[1+i*2] );
        }
      }
  }





/*
  Starting with node 0, see if an unknown is associated with
  the node.  If so, give it an index.
*/
  *nu = 0;
/*
  Handle first node.
*/
  i = 0;
  if ( ibc == 1 || ibc == 3 )
  {
    indx[i] = -1;
  }
  else
  {
    *nu = *nu + 1;
    indx[i] = *nu;
  }
/*
  Handle nodes 1 through nsub-1
*/


/*
The original program had this section written very poorly
had an unnecessary reduction, the *nu value being
written to the indx[] array could simply be replaced
with the increment i, since this loop deals with 
nodes 1 - nsub-1, *ns will simply equal nsub-1.

No wait used since the *nu assignment is not related to 
the operations in the loop hence the first available
thread can execute the statement in the single block.

After parallelisation, this block was timed to be

0.000262022s versus 0.000402927s hence there is 
significant speed up.

*/
double start = omp_get_wtime(); 
#pragma omp parallel
  {
    //int nthrds = omp_get_num_threads; 
    #pragma omp for nowait
      for ( i = 1; i < nsub; i++ )
      {
        indx[i] = i;
      }
    #pragma omp single
      {
        *nu = nsub-1;
      }
  }
double end = omp_get_wtime();
printf("Specific block time: %.16g\n",end - start);

/*
  Handle the last node.
/*/
  i = nsub;

  if ( ibc == 2 || ibc == 3 )
  {
    indx[i] = -1;
  }
  else
  {
    *nu = *nu + 1;
    indx[i] = *nu;
  }

  fprintf (fp, "\n" );
  fprintf (fp, "  Number of unknowns NU = %8d\n", *nu );
  fprintf (fp, "\n" );
  fprintf (fp, "  Node  Unknown\n" );
  fprintf (fp, "\n" );


  for ( i = 0; i <= nsub; i++ )
  {
    fprintf (fp, "  %8d  %8d\n", i, indx[i] );
  }

  fclose(fp);
double gend = omp_get_wtime();
printf("geometry() block time: %.16g\n", gend - gstart);

  return;
}
/******************************************************************************/

void assemble ( double adiag[], double aleft[], double arite[], double f[], 
  double h[], int indx[], int nl, int node[], int nu, int nquad, int nsub, 
  double ul, double ur, double xn[], double xquad[] )
{
  double aij;
  double he;
  int i;
  int ie;
  int ig;
  int il;
  int iq;
  int iu;
  int jg;
  int jl;
  int ju;
  double phii;
  double phiix;
  double phij;
  double phijx;
  double x;
  double xleft;
  double xquade;
  double xrite;
/*

Loop A1: Complete independence between loops

The following code shows a previous attempt to parallelise Loop A1, a description
of this method was shown in the report. 
*/
/*
#pragma omp parallel num_threads(NUM_THREADS)
  {
    #pragma omp sections 
    {
      #pragma omp section
        for ( i = 0; i < nu; i++ )
        {
          f[i] = 0.0;
        }
      #pragma omp section
        for ( i = 0; i < nu; i++ )
        {
          adiag[i] = 0.0;
        }
      #pragma omp section
        for ( i = 0; i < nu; i++ )
        {
          aleft[i] = 0.0;
        }
      #pragma omp section
        for ( i = 0; i < nu; i++ )
        {
          arite[i] = 0.0;
        }
    }
  }
*/

/*
It was realised that the previous author had unneccessarily split up the 
filling of matrix into four separate loops, to simplify the code, the
computation was put into one loop as shown, this loop was then parallelised.
The reasoning for rewriting the loop this way was to simplify the parallelisation
of the four loops so as to not have four separate loop work sharing constructs.

This provides a great deal of speed up for this section, the serial version gave 
a time of 0.003906, 0.003751, 0.003867 giving an average time of 0.003841s.

The parallel version give a time of 0.001520, 0.001578, 0.001610, giving an
average time of 0.001569s. This more than halves the time.
*/

double astart = omp_get_wtime();

#pragma omp parallel num_threads(NUM_THREADS)
  { 
    #pragma omp for
      for ( i = 0; i < nu; i++ )
      {
        f[i] = 0.0;
        adiag[i] = 0.0;
        aleft[i] = 0.0;
        arite[i] = 0.0;
      }
  }

/*
  For interval number IE,
*/
#pragma omp parallel private(he,xleft,xrite,xquade,ig,iu,jg,ju,aij,phii,phiix) 
  {
    #pragma omp for 
      for ( ie = 0; ie < nsub; ie++ )
      {
        int id = omp_get_thread_num();
        double he = h[ie];
        double xleft = xn[node[0+ie*2]];
        double xrite = xn[node[1+ie*2]];
       // printf("%f %f %f ie: %i id:%i\n",xrite,xleft,he,ie,id );
    /*
      consider each quadrature point IQ,
    */
        for ( iq = 0; iq < nquad; iq++ )
        {
          double xquade = xquad[ie];
         printf("%f ie: %i iq: %i id:%i\n", xquade,ie,iq ,id);
    /*
      and evaluate the integrals associated with the basis functions
      for the left, and for the right nodes.
    */
          for ( il = 1; il <= nl; il++ )
          {
            int ig = node[il-1+ie*2];
            int iu = indx[ig] - 1;
       // printf("%i %i %i\n", ig,iu,id);

            if ( 0 <= iu )
            {
              phi ( il, xquade, &phii, &phiix, xleft, xrite );//////?
              f[iu] = f[iu] + he * ff ( xquade ) * phii;//////?
             //  printf("%f\n",f[iu] );
    /*
      Take care of boundary nodes at which U' was specified.
    */
              if ( ig == 0 )
              {
                int x = 0.0;
                f[iu] = f[iu] - pp ( x ) * ul; 
                //printf("%f\n",f[iu] );
              }
              else if ( ig == nsub )
              {
                int x = 1.0;
                f[iu] = f[iu] + pp ( x ) * ur;
                //printf("%f\n",f[iu] );
              }
    /*
      Evaluate the integrals that take a product of the basis
      function times itself, or times the other basis function
      that is nonzero in this interval.
    */
              for ( jl = 1; jl <= nl; jl++ )
              {
                int jg = node[jl-1+ie*2];

                int ju = indx[jg] - 1;
//printf("jg=%i ju=%i\n",jg ,ju);
                phi ( jl, xquade, &phij, &phijx, xleft, xrite );

               double aij = he * ( pp ( xquade ) * phiix * phijx 
                           + qq ( xquade ) * phii  * phij   );
               //printf("%f\n",aij );
    /*
      If there is no variable associated with the node, then it's
      a specified boundary value, so we multiply the coefficient
      times the specified boundary value and subtract it from the
      right hand side.
    */
      //Race Conditions here
                if ( ju < 0 )
                {
                  if ( jg == 0 )
                  {
                    f[iu] = f[iu] - aij * ul;
                    //printf("%f\n",f[iu] );
                  }
                  else if ( jg == nsub )
                  {               
                    f[iu] = f[iu] - aij * ur;
                    //printf("%f\n",f[iu] );
                  }
                }
    /*
      Otherwise, we add the coefficient we've just computed to the
      diagonal, or left or right entries of row IU of the matrix.
    */
                else
                {
                  if ( iu == ju )
                  {
                    adiag[iu] = adiag[iu] + aij;
                  }
                  else if ( ju < iu )
                  {
                    aleft[iu] = aleft[iu] + aij;
                  }
                  else
                  {
                    arite[iu] = arite[iu] + aij;
                  }
                }
              }
            }
          }
        }
      }
}
  double aend = omp_get_wtime();
  printf("assemble() block time: %.16g\n",aend - astart);
  return;
 
}
/******************************************************************************/

void prsys ( double adiag[], double aleft[], double arite[], double f[], 
  int nu )
{
  double pstart = omp_get_wtime();
  int i;

  FILE *fp, *fopen();
  fp =fopen("paralleloutput.txt","a");
  

  fprintf (fp, "\n" );
  fprintf (fp, "Printout of tridiagonal linear system:\n" );
  fprintf (fp, "\n" );
  fprintf (fp, "Equation  ALEFT  ADIAG  ARITE  RHS\n" );
  fprintf (fp, "\n" );

  for ( i = 0; i < nu; i++ )
  {
    fprintf (fp, "  %8d  %14f  %14f  %14f  %14f\n",
      i + 1, aleft[i], adiag[i], arite[i], f[i] );
  }

  fclose(fp);
  double pend = omp_get_wtime();
  printf("prsys() block time: %.16g\n",pend - pstart);

  return;
}

/******************************************************************************/

void solve ( double adiag[], double aleft[], double arite[], double f[], 
  int nu )
{
  
  double sstart = omp_get_wtime();

  int i;
/*
  Carry out Gauss elimination on the matrix, saving information
  needed for the backsolve.
*/
  arite[0] = arite[0] / adiag[0];

/*
  Loop S1: Serial Section - 0.004854s
  This section has loops that have far too many dependencies between each 
  iterations, since the iterations are not independent as well as
  the relative shot duration of the serial section, this section
  will not be parallelised.
*/


  for ( i = 1; i < nu - 1; i++ )
  {
    adiag[i] = adiag[i] - aleft[i] * arite[i-1];
    arite[i] = arite[i] / adiag[i];
  }
  adiag[nu-1] = adiag[nu-1] - aleft[nu-1] * arite[nu-2];
/*
  Carry out the same elimination steps on F that were done to the
  matrix.
*/
  f[0] = f[0] / adiag[0];
  for ( i = 1; i < nu; i++ )
  {
    f[i] = ( f[i] - aleft[i] * f[i-1] ) / adiag[i];
  }
/*
  And now carry out the steps of "back substitution".
*/
  for ( i = nu - 2; 0 <= i; i-- )
  {
    f[i] = f[i] - arite[i] * f[i+1];
  }

  double send = omp_get_wtime();
  printf("solve() block time: %.16g\n",send - sstart);


  return;
}
/******************************************************************************/

void output ( double f[], int ibc, int indx[], int nsub, int nu, double ul, 
  double ur, double xn[] )

{
double ostart = omp_get_wtime();
FILE *fp, *fopen();
  fp =fopen("paralleloutput.txt","a");
  int i;
  double *u;
  u = (double *)malloc(sizeof(double)*(nsub));

  fprintf (fp, "\n" );
  fprintf (fp, "  Computed solution coefficients:\n" );
  fprintf (fp, "\n" );
  fprintf (fp, "  Node    X(I)        U(X(I))\n" );
  fprintf (fp, "\n" );
/*
Loop O1: Time - 0.07260s (Serial)
*/
#pragma omp parallel num_threads(NUM_THREADS)
{
  #pragma omp for
    for ( i = 0; i <= nsub; i++ )
    {
  /*
    If we're at the first node, check the boundary condition.
  */
      if ( i == 0 )
      {
        if ( ibc == 1 || ibc == 3 )
        {
          u[i] = ul;
        }
        else
        {
          u[i] = f[indx[i]-1];
        }
      }
  /*
    If we're at the last node, check the boundary condition.
  */
      else if ( i == nsub )
      {
        if ( ibc == 2 || ibc == 3 )
        {
          u[i] = ur;
        }
        else
        {
          u[i] = f[indx[i]-1];
        }
      }
  /*
    Any other node, we're sure the value is stored in F.
  */
      else
      {
        u[i] = f[indx[i]-1];
      }

    }

 // print statement is taken out of for loop
#pragma omp single
    {
      for ( i = 0; i <= nsub; i++ )
        {
          fprintf (fp, "  %8d  %8f  %14f\n", i, xn[i], u[i] );
        }
      fclose(fp);
    }
}
  
  double oend = omp_get_wtime();
  printf("ouput() block time: %.16g\n",oend - ostart);
  
  
  return;
}

/******************************************************************************
All functions from this point on are called in the above functions


******************************************************************************/

double ff ( double x )

{
  double value;

  value = 0.0;

  return value;
}

/******************************************************************************/

void phi ( int il, double x, double *phii, double *phiix, double xleft, 
  double xrite )
{
  if ( xleft <= x && x <= xrite )
  {
    if ( il == 1 )
    {
      *phii = ( xrite - x ) / ( xrite - xleft );
      *phiix =         -1.0 / ( xrite - xleft );
    }
    else
    {
      *phii = ( x - xleft ) / ( xrite - xleft );
      *phiix = 1.0          / ( xrite - xleft );
    }
  }
/*
  If X is outside of the interval, just set everything to 0.
*/
  else
  {
    *phii  = 0.0;
    *phiix = 0.0;
  }

  return;
}
/******************************************************************************/

double pp ( double x )
{
  double value;

  value = 1.0;

  return value;
}

/******************************************************************************/

double qq ( double x )
{
  double value;

  value = 0.0;

  return value;
}

/******************************************************************************/

void timestamp ( void )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );
  printf ( "*****************************************\n" );
  printf ( "Parallel Code: from %s\n", time_buffer );
  printf ( "*****************************************\n" );
  FILE *fp, *fopen();
  fp =fopen("paralleloutput.txt","w+");

  fprintf (fp, "%s\n", time_buffer );

  fclose(fp);

  return;
# undef TIME_SIZE
}


