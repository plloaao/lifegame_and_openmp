# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <omp.h>

#define DDefaultOutputFilename "./Life_%04d.txt"
#define DNumIterForPartialResults 25


//This code is useful to check execution times data of the code
// it prints up the ex times for init, update, write functions and the total wtime

int main(int argc, char *argv[]);
int *life_init ( char *filename, double prob, int m, int n, int *seed );
void life_update ( int m, int n, int grid[] );
void life_write ( char *output_filename, int m, int n, int grid[] );
void life_read ( char *input_filename, int m, int n, int grid[] );
double r8_uniform_01 ( int *seed );
int s_len_trim ( char *s );
void timestamp ( void );


int main(int argc, char *argv[])
{
  char *filename=NULL, output_filename[100];
  int it;
  int it_max=10;
  int m=10;
  int n=10;
  int *grid;
  double prob;
  int seed;

  double start_time, run_time;
  double start_time_init, start_time_upd, start_time_write;
  double run_time_init = 0, run_time_upd = 0, run_time_write = 0;

  start_time = omp_get_wtime();
  //timestamp ( );
//  printf ( "\n" );
//  printf ( "LIFE GAME SERIAL\n" );
//  printf ( "  C version\n" );
//  printf ( "  Carry out a few steps of John Conway's\n" );
//  printf ( "  Game of Life.\n" );
//  printf ( "  Parameters: lifegame [FicheroEstadoInicial] [TamañoX] [TamañoY] [NumIteraciones].\n" );
//  printf ( "\n" );

  if (argc>4)
    it_max = atoi(argv[4]);
  if (argc>3)
    n = atoi(argv[3]);
  if (argc>2)
    m = atoi(argv[2]);
  if (argc>1)
    filename = argv[1];
  else
    filename = NULL;

  prob = 0.20;
  seed = 123456789;

  for ( it = 0; it <= it_max; it++ )
  {
    if ( it == 0 )
    {
      start_time_init = omp_get_wtime();

      grid = life_init ( filename, prob, m, n, &seed );

      run_time_init += omp_get_wtime() - start_time_init;



    }
    else
    {
      start_time_upd = omp_get_wtime();

      life_update ( m, n, grid );

      run_time_upd += omp_get_wtime() - start_time_upd;


    }
//    if ((it%DNumIterForPartialResults)==0)
//    {
//      sprintf(output_filename,DDefaultOutputFilename,it);
//      life_write ( output_filename, m, n, grid );
//      printf ( "  %s\n", output_filename );
//    }
  }

  sprintf(output_filename,DDefaultOutputFilename,it-1); // a fine ciclo l'avrà comunque aggiornata

        start_time_write = omp_get_wtime();
  life_write ( output_filename, m, n, grid );
      run_time_write = omp_get_wtime() - start_time_write;

/*
  Free memory.
*/
  free ( grid );

/*
  Terminate.
*/
//  printf ( "\n" );
//  printf ( "LIFE GAME SERIAL\n" );
//  printf ( "  Normal end of execution.\n" );
//  printf ( "\n" );
  //timestamp ( );
        printf( "Initialization time: %fs\n", run_time_init);
              printf( "Update time: %fs\n", run_time_upd);
                           printf( "Writing time: %fs\n\n", run_time_write);
  run_time = omp_get_wtime() - start_time;
  printf("mult nxm = %dx%d for %d iter in %f seconds \n\n", n,m, it-1, run_time);

  return 0;
}


/******************************************************************************/

int *life_init ( char *filename, double prob, int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    LIFE_INIT initializes the life grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, double PROB, the probability that a grid cell
    should be alive.

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, int LIFE_INIT[(1+M+1)*(1+N+1)], the initial grid.
*/
{
  int *grid;
  int i;
  int j;
  double r;


  grid = ( int * ) malloc ( ( m + 2 ) * ( n + 2 ) * sizeof ( int ) );
  if (grid==NULL)
    perror("Error malloc grid:");

  if (filename!=NULL)
  {
      /* Read input file */
      printf("Reading Input filename %s\n",filename);
      life_read (filename, m, n, grid);
  }
  else
  {
    for ( j = 0; j <= n + 1; j++ )
    {
      for ( i = 0; i <= m + 1; i++ )
      {
        grid[i+j*(m+2)] = 0;
      }
    }

    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= m; i++ )
      {
        r = r8_uniform_01 ( seed );
        if ( r <= prob )
        {
          grid[i+j*(m+2)] = 1;
        }
      }
    }
  }

  return grid;
}
/******************************************************************************/

void life_update ( int m, int n, int grid[] )

/******************************************************************************/
/*
  Purpose:

    LIFE_UPDATE updates a Life grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input/output, int GRID[(1+M+1)*(1+N+1)], the data.
*/
{
  int i;
  int j;
  int *s;
  int i_prev, i_next, j_prev, j_next;

  s = ( int * ) malloc ( m * n * sizeof ( int ) );
#pragma omp parallel for private(i, i_prev,i_next,j_prev,j_next)
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      i_prev = (1 < i) ? i - 1 : m; // grazie a queste info posso determinare i vicini posto i-j
      i_next = (i < m ) ? i + 1 : 1;
      j_prev = (1 < j) ? j - 1 : n;
      j_next = (j < n) ? j + 1 : 1;
      s[i-1+(j-1)*m] = // in questo punto risalgo ai valori di grid al punto giusto
      //no problem poiché ricavo valori dalla shared ma senza modificarla ma sto scrivendo in una s shared
      // è un problema? 1) prova a fare private s; potrebbe funzionare fino a un certo punto per via di cache
          grid[i_prev+(j_prev)*(m+2)] + grid[i_prev+j*(m+2)] + grid[i_prev+(j_next)*(m+2)]
        + grid[i  +(j_prev)*(m+2)]                     + grid[i  +(j_next)*(m+2)]
        + grid[i_next+(j_prev)*(m+2)] + grid[i_next+j*(m+2)] + grid[i_next+(j_next)*(m+2)];
    }
  }
/*
  Any dead cell with 3 live neighbors becomes alive.
  Any living cell with less than 2 or more than 3 neighbors dies.
*/

//qui valuta delle condizioni su s quindi penso che non richieda molto ma poi gli tocca aggiornare la grid che è una shared
// 2) private grid -> si può privatizzare solo fino a un certo punto a mio avviso
// il verificare richiede tanto costo? cosa vuol dire? accedere a una variabile e verificarla??
//l'accedere a una variabile interrompe il lavoro degli altri?
//con il dynamic privatizzatimo i chunck?

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      if ( grid[i+j*(m+2)] == 0 )
      {
        if ( s[i-1+(j-1)*m] == 3 )
        {
          grid[i+j*(m+2)] = 1;
        }
      }
      else if ( grid[i+j*(m+2)] == 1 )
      {
        if ( s[i-1+(j-1)*m] < 2 || 3 < s[i-1+(j-1)*m] )
        {
          grid[i+j*(m+2)] = 0;
        }
      }
    }
  }

  free ( s );

  return;
}
/******************************************************************************/

void life_write ( char *output_filename, int m, int n, int grid[] )

/******************************************************************************/
/*
  Purpose:

    LIFE_WRITE writes a grid to a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output file name.

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input, int GRID[(1+M+1)*(1+N+1)], the data.
*/
{
  int i;
  int j;
  FILE *output_unit;
/*
  Open the file.
*/
  output_unit = fopen ( output_filename, "wt" );
/*
  Write the data.
*/
  for ( j = 1; j <= n ; j++ )
  {
    for ( i = 1; i <= m ; i++ )
    {
      fprintf ( output_unit, " %d", grid[i+j*(m+2)] );
    }
    fprintf ( output_unit, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output_unit );

  return;
}
/******************************************************************************/

void life_read ( char *filename, int m, int n, int grid[] )

/******************************************************************************/
/*
  Purpose:

    LIFE_READ reads a file to a grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Parameters:

    Input, char *OUTPUT_FILENAME, the output file name.

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input, int GRID[(1+M+1)*(1+N+1)], the data.
*/
{
  int i;
  int j;
  FILE *input_unit;
/*
  input the file.
*/
  input_unit = fopen ( filename, "rt" );
  if (input_unit==NULL)
    perror("Reading input file:");
/*
  Read the data.
*/
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      fscanf ( input_unit, "%d", &(grid[i+j*(m+2)]) );
    }
  }
  /* Set the grid borderline to 0's */
  for ( j = 0; j <= n +1; j++ )
  {
    grid[0+j*(m+2)] = 0;
    grid[(m+1)+j*(m+2)] = 0;

  }
  for ( i = 1; i <= m; i++ )
  {
    grid[i+0*(m+2)] = 0;
    grid[i+(n+1)*(m+2)] = 0;
  }
/*
  Close the file.
*/
  fclose ( input_unit );

  return;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a pseudorandom R8 scaled to [0,1].

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    P A Lewis, A S Goodman, J M Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int i4_huge = 2147483647;
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

int s_len_trim ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LEN_TRIM returns the length of a string to the last nonblank.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 April 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
    If S_LEN_TRIM is 0, then the string is entirely blank.
*/
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
//# define TIME_SIZE 40
//
//  static char time_buffer[TIME_SIZE];
//  const struct tm *tm;
//  size_t len;
//  time_t now;
//
//  now = time ( NULL );
//  tm = localtime ( &now );
//
//  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );
//
//  fprintf ( stdout, "%s\n", time_buffer );
//
//  return;
//# undef TIME_SIZE
}

