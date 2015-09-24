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
  
#pragma omp parallel for num_threads(NUM_THREADS) 
  for ( i = 0; i <= nsub; i++ )
  {
    xn[i]  =  ( ( double ) ( nsub - i ) * xl 
              + ( double )          i   * xr ) 
              / ( double ) ( nsub );
  }

  for ( i = 0; i <= nsub; i++ )
    {
      fprintf (fp, "  %8d  %14f \n", i, xn[i] );
    }
    
/*
  Set the lengths of each subinterval.
  First thread that gets here, print out statements, 
  implied barrier at the end of block
*/

          fprintf (fp, "\n" );
          fprintf (fp, "Subint    Length\n" );
          fprintf (fp, "\n" );


//Loop G2

    for ( i = 0; i < nsub; i++ )
    {
      h[i] = xn[i+1] - xn[i];
      fprintf (fp, "  %8d  %14f\n", i+1, h[i] );
    }
/*  
  Set the quadrature points, each of which is the midpoint
  of its subinterval.
*/

      fprintf (fp, "\n" );
      fprintf (fp, "Subint    Quadrature point\n" );
      fprintf (fp, "\n" );

  
 //Loop G3
  for ( i = 0; i < nsub; i++ )
  {
    xquad[i] = 0.5 * ( xn[i] + xn[i+1] );
   fprintf (fp, "  %8d  %14f\n", i+1, xquad[i] );
  }
/*
  Set the value of NODE, which records, for each interval,
  the node numbers at the left and right.
*/
  fprintf (fp, "\n" );
  fprintf (fp, "Subint  Left Node  Right Node\n" );
  fprintf (fp, "\n" );
 //Loop G4
  for ( i = 0; i < nsub; i++ )
  {
    node[0+i*2] = i;
    node[1+i*2] = i + 1;
    fprintf (fp, "  %8d  %8d  %8d\n", i+1, node[0+i*2], node[1+i*2] );
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
  for ( i = 1; i < nsub; i++ )
  {
    *nu = *nu + 1;
    indx[i] = *nu;
  }
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

#pragma omp parallel num_threads(NUM_THREADS)
  { 
    #pragma omp for nowait
      for ( i = 0; i < nu; i++ )
      {
        f[i] = 0.0;
      }
    #pragma omp for nowait
      for ( i = 0; i < nu; i++ )
      {
        adiag[i] = 0.0;
      }
    #pragma omp for nowait
      for ( i = 0; i < nu; i++ )
      {
        aleft[i] = 0.0;
      }
     #pragma omp for nowait 
      for ( i = 0; i < nu; i++ )
      {
        arite[i] = 0.0;
      }
    }

/*
  For interval number IE,
*/

  for ( ie = 0; ie < nsub; ie++ )
  {
    he = h[ie];
    xleft = xn[node[0+ie*2]];
    xrite = xn[node[1+ie*2]];
/*
  consider each quadrature point IQ,
*/
    for ( iq = 0; iq < nquad; iq++ )
    {
      xquade = xquad[ie];
/*
  and evaluate the integrals associated with the basis functions
  for the left, and for the right nodes.
*/
      for ( il = 1; il <= nl; il++ )
      {
        ig = node[il-1+ie*2];
        iu = indx[ig] - 1;

        if ( 0 <= iu )
        {
          phi ( il, xquade, &phii, &phiix, xleft, xrite );
          f[iu] = f[iu] + he * ff ( xquade ) * phii;
/*
  Take care of boundary nodes at which U' was specified.
*/
          if ( ig == 0 )
          {
            x = 0.0;
            f[iu] = f[iu] - pp ( x ) * ul;
          }
          else if ( ig == nsub )
          {
            x = 1.0;
            f[iu] = f[iu] + pp ( x ) * ur;
          }
/*
  Evaluate the integrals that take a product of the basis
  function times itself, or times the other basis function
  that is nonzero in this interval.
*/
          for ( jl = 1; jl <= nl; jl++ )
          {
            jg = node[jl-1+ie*2];
            ju = indx[jg] - 1;

            phi ( jl, xquade, &phij, &phijx, xleft, xrite );

            aij = he * ( pp ( xquade ) * phiix * phijx 
                       + qq ( xquade ) * phii  * phij   );
/*
  If there is no variable associated with the node, then it's
  a specified boundary value, so we multiply the coefficient
  times the specified boundary value and subtract it from the
  right hand side.
*/
            if ( ju < 0 )
            {
              if ( jg == 0 )
              {
                f[iu] = f[iu] - aij * ul;
              }
              else if ( jg == nsub )
              {               
                f[iu] = f[iu] - aij * ur;
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
  return;
}
/******************************************************************************/

void prsys ( double adiag[], double aleft[], double arite[], double f[], 
  int nu )
{
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

  return;
}

/******************************************************************************/

void solve ( double adiag[], double aleft[], double arite[], double f[], 
  int nu )
{
  int i;
/*
  Carry out Gauss elimination on the matrix, saving information
  needed for the backsolve.
*/
  arite[0] = arite[0] / adiag[0];

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

  return;
}
/******************************************************************************/

void output ( double f[], int ibc, int indx[], int nsub, int nu, double ul, 
  double ur, double xn[] )

{
FILE *fp, *fopen();
  fp =fopen("paralleloutput.txt","a");
  int i;
  double u;

  fprintf (fp, "\n" );
  fprintf (fp, "  Computed solution coefficients:\n" );
  fprintf (fp, "\n" );
  fprintf (fp, "  Node    X(I)        U(X(I))\n" );
  fprintf (fp, "\n" );

  for ( i = 0; i <= nsub; i++ )
  {
/*
  If we're at the first node, check the boundary condition.
*/
    if ( i == 0 )
    {
      if ( ibc == 1 || ibc == 3 )
      {
        u = ul;
      }
      else
      {
        u = f[indx[i]-1];
      }
    }
/*
  If we're at the last node, check the boundary condition.
*/
    else if ( i == nsub )
    {
      if ( ibc == 2 || ibc == 3 )
      {
        u = ur;
      }
      else
      {
        u = f[indx[i]-1];
      }
    }
/*
  Any other node, we're sure the value is stored in F.
*/
    else
    {
      u = f[indx[i]-1];
    }

    fprintf (fp, "  %8d  %8f  %14f\n", i, xn[i], u );
  }

  fclose(fp);
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


