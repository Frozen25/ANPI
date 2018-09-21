


#include <iostream>
#include <fstream>
#include <cctype>
#include <cmath>
#include <cfloat>
#include "stdafx.h"
#include "math.h"
#include <float.h>



template <class T, class U>
void jenkinsTraub(const bmt::polynomial<T>& poly,
                    std::vector<U>& roots) {
    
    static_assert(std::is_floating_point<T>::value ||
                  boost::is_complex<T>::value,
                  "T must be floating point or complex");
    static_assert(std::is_floating_point<U>::value ||
                  boost::is_complex<U>::value,
                  "U must be floating point or complex");

   


}




static T sr, si, tr, ti, pvr, pvi, are, mre, eta, infin;
static int nn;
static T *pr, *pi, *hr, *hi, *qpr, *qpi, *qhr, *qhi, *shr, *shi; 

static void noshft( const int l1 );
static void fxshft( const int l2, T *zr, T *zi, int *conv );
static void vrshft( const int l3, T *zr, T *zi, int *conv );
static void calct( int *bol );
static void nexth( const int bol );
static void polyev( const int nn, const T sr, const T si, const T pr[], const T pi[], T qr[], T qi[], T *pvr, T *pvi );
static T errev( const int nn, const T qr[], const T qi[], const T ms, const T mp, const T are, const T mre );
static void cauchy( const int nn, T pt[], T q[], T *fn_val );
static T scale( const int nn, const T pt[], const T eta, const T infin, const T smalno, const T base );
static void cdivid( const T ar, const T ai, const T br, const T bi, T *cr, T *ci );
static T cmod( const T r, const T i );
static void mcon( T *eta, T *infiny, T *smalno, T *base );

template <class T>
int cpoly( const T *opr, const T *opi, int degree, T *zeror, T *zeroi ) 
   {
   int cnt1, cnt2, idnn2, i, conv;
   T xx, yy, cosr, sinr, smalno, base, xxx, zr, zi, bnd;

   mcon( &eta, &infin, &smalno, &base );
   are = eta;
   mre = 2.0 * sqrt( 2.0 ) * eta;
   xx = 0.70710678;
   yy = -xx;
   cosr = -0.060756474;
   sinr = -0.99756405;
   nn = degree;  

   // Algorithm fails if the leading coefficient is zero
   if( opr[ 0 ] == 0 && opi[ 0 ] == 0 )
      return -1;

   // Allocate arrays
   pr = new T [ degree+1 ];
   pi = new T [ degree+1 ];
   hr = new T [ degree+1 ];
   hi = new T [ degree+1 ];
   qpr= new T [ degree+1 ];
   qpi= new T [ degree+1 ];
   qhr= new T [ degree+1 ];
   qhi= new T [ degree+1 ];
   shr= new T [ degree+1 ];
   shi= new T [ degree+1 ];

   // Remove the zeros at the origin if any
   while( opr[ nn ] == 0 && opi[ nn ] == 0 )
      {
      idnn2 = degree - nn;
      zeror[ idnn2 ] = 0;
      zeroi[ idnn2 ] = 0;
      nn--;
      }

   // Make a copy of the coefficients
   for( i = 0; i <= nn; i++ )
      {
      pr[ i ] = opr[ i ];
      pi[ i ] = opi[ i ];
      shr[ i ] = cmod( pr[ i ], pi[ i ] );
      }

   // Scale the polynomial
   bnd = scale( nn, shr, eta, infin, smalno, base );
   if( bnd != 1 )
      for( i = 0; i <= nn; i++ )
         {
         pr[ i ] *= bnd;
         pi[ i ] *= bnd;
         }

search: 
   if( nn <= 1 )
      {
      cdivid( -pr[ 1 ], -pi[ 1 ], pr[ 0 ], pi[ 0 ], &zeror[ degree-1 ], &zeroi[ degree-1 ] );
      goto finish;
      }

   // Calculate bnd, alower bound on the modulus of the zeros
   for( i = 0; i<= nn; i++ )
      shr[ i ] = cmod( pr[ i ], pi[ i ] );

   cauchy( nn, shr, shi, &bnd );
   
   // Outer loop to control 2 Major passes with different sequences of shifts
   for( cnt1 = 1; cnt1 <= 2; cnt1++ )
      {
      // First stage  calculation , no shift
      noshft( 5 );

      // Inner loop to select a shift
      for( cnt2 = 1; cnt2 <= 9; cnt2++ )
         {
         // Shift is chosen with modulus bnd and amplitude rotated by 94 degree from the previous shif
         xxx = cosr * xx - sinr * yy;
         yy = sinr * xx + cosr * yy;
         xx = xxx;
         sr = bnd * xx;
         si = bnd * yy;

         // Second stage calculation, fixed shift
         fxshft( 10 * cnt2, &zr, &zi, &conv );
         if( conv )
            {
            // The second stage jumps directly to the third stage ieration
            // If successful the zero is stored and the polynomial deflated
            idnn2 = degree - nn;
            zeror[ idnn2 ] = zr;
            zeroi[ idnn2 ] = zi;
            nn--;
            for( i = 0; i <= nn; i++ )
               {
               pr[ i ] = qpr[ i ];
               pi[ i ] = qpi[ i ];
               }
            goto search;
            }
         // If the iteration is unsuccessful another shift is chosen
         }
      // if 9 shifts fail, the outer loop is repeated with another sequence of shifts
      }

   // The zerofinder has failed on two major passes
   // return empty handed with the number of roots found (less than the original degree)
   degree -= nn;

finish:
   // Deallocate arrays
   delete [] pr;
   delete [] pi;
   delete [] hr;
   delete [] hi;
   delete [] qpr;
   delete [] qpi;
   delete [] qhr;
   delete [] qhi;
   delete [] shr;
   delete [] shi;

   return degree;       
   }


// COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
// POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.
//
template <class T>
static void noshft( const int l1 )
   {
   int i, j, jj, n, nm1;
   T xni, t1, t2;

   n = nn;
   nm1 = n - 1;
   for( i = 0; i < n; i++ )
      {
      xni = nn - i;
      hr[ i ] = xni * pr[ i ] / n;
      hi[ i ] = xni * pi[ i ] / n;
      }
   for( jj = 1; jj <= l1; jj++ )
      {
      if( cmod( hr[ n - 1 ], hi[ n - 1 ] ) > eta * 10 * cmod( pr[ n - 1 ], pi[ n - 1 ] ) )
         {
         cdivid( -pr[ nn ], -pi[ nn ], hr[ n - 1 ], hi[ n - 1 ], &tr, &ti );
         for( i = 0; i < nm1; i++ )
            {
            j = nn - i - 1;
            t1 = hr[ j - 1 ];
            t2 = hi[ j - 1 ];
            hr[ j ] = tr * t1 - ti * t2 + pr[ j ];
            hi[ j ] = tr * t2 + ti * t1 + pi[ j ];
            }
         hr[ 0 ] = pr[ 0 ];
         hi[ 0 ] = pi[ 0 ];
         }
      else
         {
         // If the constant term is essentially zero, shift H coefficients
         for( i = 0; i < nm1; i++ )
            {
            j = nn - i - 1;
            hr[ j ] = hr[ j - 1 ];
            hi[ j ] = hi[ j - 1 ];
            }
         hr[ 0 ] = 0;
         hi[ 0 ] = 0;
         }
      }
   }

// COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR CONVERGENCE.
// INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
// APPROXIMATE ZERO IF SUCCESSFUL.
// L2 - LIMIT OF FIXED SHIFT STEPS
// ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
// CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
//
template <class T>
static void fxshft( const int l2, T *zr, T *zi, int *conv )
   {
   int i, j, n;
   int test, pasd, bol;
   T otr, oti, svsr, svsi;

   n = nn;
   polyev( nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi );
   test = 1;
   pasd = 0;

   // Calculate first T = -P(S)/H(S)
   calct( &bol );

   // Main loop for second stage
   for( j = 1; j <= l2; j++ )
      {
      otr = tr;
      oti = ti;

      // Compute the next H Polynomial and new t
      nexth( bol );
      calct( &bol );
      *zr = sr + tr;
      *zi = si + ti;

      // Test for convergence unless stage 3 has failed once or this
      // is the last H Polynomial
      if( !( bol || !test || j == 12 ) )
         if( cmod( tr - otr, ti - oti ) < 0.5 * cmod( *zr, *zi ) )
            {
            if( pasd )
               {
               // The weak convergence test has been passwed twice, start the third stage
               // Iteration, after saving the current H polynomial and shift
               for( i = 0; i < n; i++ )
                  {
                  shr[ i ] = hr[ i ];
                  shi[ i ] = hi[ i ];
                  }
               svsr = sr;
               svsi = si;
               vrshft( 10, zr, zi, conv );
               if( *conv ) return;

               //The iteration failed to converge. Turn off testing and restore h,s,pv and T
               test = 0;
               for( i = 0; i < n; i++ )
                  {
                  hr[ i ] = shr[ i ];
                  hi[ i ] = shi[ i ];
                  }
               sr = svsr;
               si = svsi;
               polyev( nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi );
               calct( &bol );
               continue;
               }
            pasd = 1;
            }
         else
            pasd = 0;
      }

   // Attempt an iteration with final H polynomial from second stage
   vrshft( 10, zr, zi, conv );
   }

// CARRIES OUT THE THIRD STAGE ITERATION.
// L3 - LIMIT OF STEPS IN STAGE 3.
// ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
//           ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE ON EXIT.
// CONV    -  .TRUE. IF ITERATION CONVERGES
//
template <class T>
static void vrshft( const int l3, T *zr, T *zi, int *conv )
   {
   int b, bol;
   int i, j;
   T mp, ms, omp, relstp, r1, r2, tp;

   *conv = 0;
   b = 0;
   sr = *zr;
   si = *zi;

   // Main loop for stage three
   for( i = 1; i <= l3; i++ )
      {
      // Evaluate P at S and test for convergence
      polyev( nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi );
      mp = cmod( pvr, pvi );
      ms = cmod( sr, si );
      if( mp <= 20 * errev( nn, qpr, qpi, ms, mp, are, mre ) )
         {
         // Polynomial value is smaller in value than a bound onthe error
         // in evaluationg P, terminate the ietartion
         *conv = 1;
         *zr = sr;
         *zi = si;
         return;
         }
      if( i != 1 )
         {
         if( !( b || mp < omp || relstp >= 0.05 ) )
            {
            // Iteration has stalled. Probably a cluster of zeros. Do 5 fixed 
            // shift steps into the cluster to force one zero to dominate
            tp = relstp;
            b = 1;
            if( relstp < eta ) tp = eta;
            r1 = sqrt( tp );
            r2 = sr * ( 1 + r1 ) - si * r1;
            si = sr * r1 + si * ( 1 + r1 );
            sr = r2;
            polyev( nn, sr, si, pr, pi, qpr, qpi, &pvr, &pvi );
            for( j = 1; j <= 5; j++ )
               {
               calct( &bol );
               nexth( bol );
               }
            omp = infin;
            goto _20;
            }
         
         // Exit if polynomial value increase significantly
         if( mp *0.1 > omp ) return;
         }

      omp = mp;

      // Calculate next iterate
_20:  calct( &bol );
      nexth( bol );
      calct( &bol );
      if( !bol )
         {
         relstp = cmod( tr, ti ) / cmod( sr, si );
         sr += tr;
         si += ti;
         }
      }
   }

// COMPUTES  T = -P(S)/H(S).
// BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
template <class T>
static void calct( int *bol )
   {
   int n;
   T hvr, hvi;

   n = nn;

   // evaluate h(s)
   polyev( n - 1, sr, si, hr, hi, qhr, qhi, &hvr, &hvi );
   *bol = cmod( hvr, hvi ) <= are * 10 * cmod( hr[ n - 1 ], hi[ n - 1 ] ) ? 1 : 0;
   if( !*bol )
      {
      cdivid( -pvr, -pvi, hvr, hvi, &tr, &ti );
      return;
      }

   tr = 0;
   ti = 0;
   }

// CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
// BOOL   -  LOGICAL, IF .TRUE. H(S) IS ESSENTIALLY ZERO
//
template <class T>
static void nexth( const int bol )
   {
   int j, n;
   T t1, t2;

   n = nn;
   if( !bol )
      {
      for( j = 1; j < n; j++ )
         {
         t1 = qhr[ j - 1 ];
         t2 = qhi[ j - 1 ];
         hr[ j ] = tr * t1 - ti * t2 + qpr[ j ];
         hi[ j ] = tr * t2 + ti * t1 + qpi[ j ];
         }
      hr[ 0 ] = qpr[ 0 ];
      hi[ 0 ] = qpi[ 0 ];
      return;
      }

   // If h[s] is zero replace H with qh
   for( j = 1; j < n; j++ )
      {
      hr[ j ] = qhr[ j - 1 ];
      hi[ j ] = qhi[ j - 1 ];
      }
   hr[ 0 ] = 0;
   hi[ 0 ] = 0;
   }

// EVALUATES A POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
// PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
//  
template <class T>
static void polyev( const int nn, const T sr, const T si, const T pr[], const T pi[], T qr[], T qi[], T *pvr, T *pvi )  
   {
   int i;
   T t;

   qr[ 0 ] = pr[ 0 ];
   qi[ 0 ] = pi[ 0 ];
   *pvr = qr[ 0 ];
   *pvi = qi[ 0 ];

   for( i = 1; i <= nn; i++ )
      {
      t = ( *pvr ) * sr - ( *pvi ) * si + pr[ i ];
      *pvi = ( *pvr ) * si + ( *pvi ) * sr + pi[ i ];
      *pvr = t;
      qr[ i ] = *pvr;
      qi[ i ] = *pvi;
      }
   }

// BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER RECURRENCE.
// QR,QI - THE PARTIAL SUMS
// MS    -MODULUS OF THE POINT
// MP    -MODULUS OF POLYNOMIAL VALUE
// ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION
//
template <class T>
static T errev( const int nn, const T qr[], const T qi[], const T ms, const T mp, const T are, const T mre )
   {
   int i;
   T e;

   e = cmod( qr[ 0 ], qi[ 0 ] ) * mre / ( are + mre );
   for( i = 0; i <= nn; i++ )
      e = e * ms + cmod( qr[ i ], qi[ i ] );

   return e * ( are + mre ) - mp * mre;
   }

// CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
// POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.
//
template <class T>
static void cauchy( const int nn, T pt[], T q[], T *fn_val )
   {
   int i, n;
   T x, xm, f, dx, df;

   pt[ nn ] = -pt[ nn ];

   // Compute upper estimate bound
   n = nn;
   x = exp( log( -pt[ nn ] ) - log( pt[ 0 ] ) ) / n;
   if( pt[ n - 1 ] != 0 )
      {
      // Newton step at the origin is better, use it
      xm = -pt[ nn ] / pt[ n - 1 ];
      if( xm < x ) x = xm;
      }

   // Chop the interval (0,x) until f < 0
   while(1)
      {
      xm = x * 0.1;
      f = pt[ 0 ];
      for( i = 1; i <= nn; i++ )
         f = f * xm + pt[ i ];
      if( f <= 0 )
         break;
      x = xm;
      }
   dx = x;
   
   // Do Newton iteration until x converges to two decimal places
   while( fabs( dx / x ) > 0.005 )
      {
      q[ 0 ] = pt[ 0 ];
      for( i = 1; i <= nn; i++ )
         q[ i ] = q[ i - 1 ] * x + pt[ i ];
      f = q[ nn ];
      df = q[ 0 ];
      for( i = 1; i < n; i++ )
         df = df * x + q[ i ];
      dx = f / df;
      x -= dx;
      }

   *fn_val = x;
   }

// RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL.
// THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW
// INTERFERING WITH THE CONVERGENCE CRITERION.  THE FACTOR IS A POWER OF THE
// BASE.
// PT - MODULUS OF COEFFICIENTS OF P
// ETA, INFIN, SMALNO, BASE - CONSTANTS DESCRIBING THE FLOATING POINT ARITHMETIC.
//
template <class T>   
static T scale( const int nn, const T pt[], const T eta, const T infin, const T smalno, const T base )
   {
   int i, l;
   T hi, lo, max, min, x, sc;
   T fn_val;

   // Find largest and smallest moduli of coefficients
   hi = sqrt( infin );
   lo = smalno / eta;
   max = 0;
   min = infin;

   for( i = 0; i <= nn; i++ )
      {
      x = pt[ i ];
      if( x > max ) max = x;
      if( x != 0 && x < min ) min = x;
      }

   // Scale only if there are very large or very small components
   fn_val = 1;
   if( min >= lo && max <= hi ) return fn_val;
   x = lo / min;
   if( x <= 1 )
      sc = 1 / ( sqrt( max )* sqrt( min ) );
   else
      {
      sc = x;
      if( infin / sc > max ) sc = 1;
      }
   l = (int)( log( sc ) / log(base ) + 0.5 );
   fn_val = pow( base , l );
   return fn_val;
   }

// COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.
//
template <class T>
static void cdivid( const T ar, const T ai, const T br, const T bi, T *cr, T *ci )
   {
   T r, d, t, infin;

   if( br == 0 && bi == 0 )
      {
      // Division by zero, c = infinity
      mcon( &t, &infin, &t, &t );
      *cr = infin;
      *ci = infin;
      return;
      }

   if( fabs( br ) < fabs( bi ) )
      {
      r = br/ bi;
      d = bi + r * br;
      *cr = ( ar * r + ai ) / d;
      *ci = ( ai * r - ar ) / d;
      return;
      }

   r = bi / br;
   d = br + r * bi;
   *cr = ( ar + ai * r ) / d;
   *ci = ( ai - ar * r ) / d;
   }

// MODULUS OF A COMPLEX NUMBER AVOIDING OVERFLOW.
//
template <class T>
static T cmod( const T r, const T i )
   {
   T ar, ai;

   ar = fabs( r );
   ai = fabs( i );
   if( ar < ai )
      return ai * sqrt( 1.0 + pow( ( ar / ai ), 2.0 ) );

   if( ar > ai )
      return ar * sqrt( 1.0 + pow( ( ai / ar ), 2.0 ) );

   return ar * sqrt( 2.0 );
   }

// MCON PROVIDES MACHINE CONSTANTS USED IN VARIOUS PARTS OF THE PROGRAM.
// THE USER MAY EITHER SET THEM DIRECTLY OR USE THE STATEMENTS BELOW TO
// COMPUTE THEM. THE MEANING OF THE FOUR CONSTANTS ARE -
// ETA       THE MAXIMUM RELATIVE REPRESENTATION ERROR WHICH CAN BE DESCRIBED
//           AS THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
//           1.0_dp + ETA &gt; 1.0.
// INFINY    THE LARGEST FLOATING-POINT NUMBER
// SMALNO    THE SMALLEST POSITIVE FLOATING-POINT NUMBER
// BASE      THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED
//
template <class T>
static void mcon( T *eta, T *infiny, T *smalno, T *base )
   {
   *base = DBL_RADIX;
   *eta = DBL_EPSILON;
   *infiny = DBL_MAX;
   *smalno = DBL_MIN;
   }

///////////////////////////////////////////////////////////
////                    REAL 
///////////////////////////////////////////////////////////
using namespace std{

    #define MAXDEGREE 100
    #define MDP1 MAXDEGREE+1

    void rpoly_ak1(T op[MDP1], int* Degree, T zeror[MAXDEGREE], T zeroi[MAXDEGREE]);
    void Fxshfr_ak1(int L2, int* NZ, T sr, T bnd, T K[MDP1], int N, T p[MDP1], int NN, T qp[MDP1], T* lzi, T* lzr, T* szi, T* szr);
    void QuadSD_ak1(int NN, T u, T v, T p[MDP1], T q[MDP1], T* a, T* b);
    int calcSC_ak1(int N, T a, T b, T* a1, T* a3, T* a7, T* c, T* d, T* e, T* f, T* g, T* h, T K[MDP1], T u, T v, T qk[MDP1]);
    void nextK_ak1(int N, int tFlag, T a, T b, T a1, T* a3, T* a7, T K[MDP1], T qk[MDP1], T qp[MDP1]);
    void newest_ak1(int tFlag, T* uu, T* vv, T a, T a1, T a3, T a7, T b, T c, T d, T f, T g, T h, T u, T v, T K[MDP1], int N, T p[MDP1]);
    void QuadIT_ak1(int N, int* NZ, T uu, T vv, T* szr, T* szi, T* lzr, T* lzi, T qp[MDP1], int NN, T* a, T* b, T p[MDP1], T qk[MDP1], T* a1, T* a3, T* a7, T* d, T* e, T* f, T* g, T* h, T K[MDP1]);
    void RealIT_ak1(int* iFlag, int* NZ, T* sss, int N, T p[MDP1], int NN, T qp[MDP1], T* szr, T* szi, T K[MDP1], T qk[MDP1]);
    void Quad_ak1(T a, T b1, T c, T* sr, T* si, T* lr, T* li);

    template <class T>
    void rpoly_ak1(T op[MDP1], int* Degree, T zeror[MAXDEGREE], T zeroi[MAXDEGREE]){

        int i, j, jj, l, N, NM1, NN, NZ, zerok;

        T K[MDP1], p[MDP1], pt[MDP1], qp[MDP1], temp[MDP1];
        T bnd, df, dx, factor, ff, moduli_max, moduli_min, sc, x, xm;
        T aa, bb, cc, lzi, lzr, sr, szi, szr, t, xx, xxx, yy;

        const T RADFAC = 3.14159265358979323846/180; // Degrees-to-radians conversion factor = pi/180
        const T lb2 = log(2.0); // Dummy variable to avoid re-calculating this value in loop below
        const T lo = FLT_MIN/DBL_EPSILON;
        const T cosr = cos(94.0*RADFAC); // = -0.069756474
        const T sinr = sin(94.0*RADFAC); // = 0.99756405

        if ((*Degree) > MAXDEGREE){
            cout << "\nThe entered Degree is greater than MAXDEGREE. Exiting rpoly. No further action taken.\n";
            *Degree = -1;
            return;
        } // End ((*Degree) > MAXDEGREE)

        //Do a quick check to see if leading coefficient is 0
        if (op[0] != 0){

        N = *Degree;
        xx = sqrt(0.5); // = 0.70710678
        yy = -xx;

        // Remove zeros at the origin, if any
        j = 0;
        while (op[N] == 0){
            zeror[j] = zeroi[j] = 0.0;
            N--;
            j++;
        } // End while (op[N] == 0)

        NN = N + 1;

        // Make a copy of the coefficients
        for (i = 0; i < NN; i++)   p[i] = op[i];

        while (N >= 1){ // Main loop
            // Start the algorithm for one zero
            if (N <= 2){
            // Calculate the final zero or pair of zeros
                if (N < 2){
                    zeror[(*Degree) - 1] = -(p[1]/p[0]);
                    zeroi[(*Degree) - 1] = 0.0;
                } // End if (N < 2)
                else { // else N == 2
                    Quad_ak1(p[0], p[1], p[2], &zeror[(*Degree) - 2], &zeroi[(*Degree) - 2], &zeror[(*Degree) - 1], &zeroi[(*Degree) - 1]);
                } // End else N == 2
                break;
            } // End if (N <= 2)

            // Find the largest and smallest moduli of the coefficients

            moduli_max = 0.0;
            moduli_min = FLT_MAX;

            for (i = 0; i < NN; i++){
                x = fabs(p[i]);
                if (x > moduli_max)   moduli_max = x;
                if ((x != 0) && (x < moduli_min))   moduli_min = x;
            } // End for i

            // Scale if there are large or very small coefficients
            // Computes a scale factor to multiply the coefficients of the polynomial. The scaling
            // is done to avoid overflow and to avoid undetected underflow interfering with the
            // convergence criterion.
            // The factor is a power of the base.

            sc = lo/moduli_min;

            if (((sc <= 1.0) && (moduli_max >= 10)) || ((sc > 1.0) && (FLT_MAX/sc >= moduli_max))){
                sc = ((sc == 0) ? FLT_MIN : sc);
                l = (int)(log(sc)/lb2 + 0.5);
                factor = pow(2.0, l);
                if (factor != 1.0){
                    for (i = 0; i < NN; i++)   p[i] *= factor;
                } // End if (factor != 1.0)
            } // End if (((sc <= 1.0) && (moduli_max >= 10)) || ((sc > 1.0) && (FLT_MAX/sc >= moduli_max)))

            // Compute lower bound on moduli of zeros

            for (i = 0; i < NN; i++)   pt[i] = fabs(p[i]);
            pt[N] = -(pt[N]);

            NM1 = N - 1;

            // Compute upper estimate of bound

            x = exp((log(-pt[N]) - log(pt[0]))/(T)N);

            if (pt[NM1] != 0) {
                // If Newton step at the origin is better, use it
                xm = -pt[N]/pt[NM1];
                x = ((xm < x) ? xm : x);
            } // End if (pt[NM1] != 0)

            // Chop the interval (0, x) until ff <= 0

            xm = x;
            do {
                x = xm;
                xm = 0.1*x;
                ff = pt[0];
                for (i = 1; i < NN; i++)   ff = ff *xm + pt[i];
            } while (ff > 0); // End do-while loop

            dx = x;

            // Do Newton iteration until x converges to two decimal places

            while (fabs(dx/x) > 0.005) {
                df = ff = pt[0];
                for (i = 1; i < N; i++){
                    ff = x*ff + pt[i];
                    df = x*df + ff;
                } // End for i
                ff = x*ff + pt[N];
                dx = ff/df;
                x -= dx;
            } // End while loop

            bnd = x;

            // Compute the derivative as the initial K polynomial and do 5 steps with no shift

            for (i = 1; i < N; i++)   K[i] = (T)(N - i)*p[i]/((T)N);
            K[0] = p[0];

            aa = p[N];
            bb = p[NM1];
            zerok = ((K[NM1] == 0) ? 1 : 0);

            for (jj = 0; jj < 5; jj++) {
                cc = K[NM1];
                if (zerok){
                    // Use unscaled form of recurrence
                    for (i = 0; i < NM1; i++){
                        j = NM1 - i;
                        K[j] = K[j - 1];
                    } // End for i
                    K[0] = 0;
                   zerok = ((K[NM1] == 0) ? 1 : 0);
                } // End if (zerok)

                else { // else !zerok
                    // Used scaled form of recurrence if value of K at 0 is nonzero
                    t = -aa/cc;
                    for (i = 0; i < NM1; i++){
                        j = NM1 - i;
                        K[j] = t*K[j - 1] + p[j];
                    } // End for i
                    K[0] = p[0];
                    zerok = ((fabs(K[NM1]) <= fabs(bb)*DBL_EPSILON*10.0) ? 1 : 0);
                } // End else !zerok

            } // End for jj

            // Save K for restarts with new shifts
            for (i = 0; i < N; i++)   temp[i] = K[i];

            // Loop to select the quadratic corresponding to each new shift

            for (jj = 1; jj <= 20; jj++){

                // Quadratic corresponds to a T shift to a non-real point and its
                // complex conjugate. The point has modulus BND and amplitude rotated
                // by 94 degrees from the previous shift.

                xxx = -(sinr*yy) + cosr*xx;
                yy = sinr*xx + cosr*yy;
                xx = xxx;
                sr = bnd*xx;

                // Second stage calculation, fixed quadratic

                Fxshfr_ak1(20*jj, &NZ, sr, bnd, K, N, p, NN, qp, &lzi, &lzr, &szi, &szr);

                if (NZ != 0){

                    // The second stage jumps directly to one of the third stage iterations and
                    // returns here if successful. Deflate the polynomial, store the zero or
                    // zeros, and return to the main algorithm.

                    j = (*Degree) - N;
                    zeror[j] = szr;
                    zeroi[j] = szi;
                    NN = NN - NZ;
                    N = NN - 1;
                    for (i = 0; i < NN; i++)   p[i] = qp[i];
                    if (NZ != 1){
                        zeror[j + 1] = lzr;
                        zeroi[j + 1] = lzi;
                    } // End if (NZ != 1)
                    break;
                } // End if (NZ != 0)
                else { // Else (NZ == 0)

                    // If the iteration is unsuccessful, another quadratic is chosen after restoring K
                    for (i = 0; i < N; i++)   K[i] = temp[i];
                } // End else (NZ == 0)

            } // End for jj

            // Return with failure if no convergence with 20 shifts

            if (jj > 20) {
                cout << "\nFailure. No convergence after 20 shifts. Program terminated.\n";
                *Degree -= N;
                break;
            } // End if (jj > 20)

        } // End while (N >= 1)

        } // End if op[0] != 0
        else { // else op[0] == 0
            cout << "\nThe leading coefficient is zero. No further action taken. Program terminated.\n";
            *Degree = 0;
        } // End else op[0] == 0

        return;
    } // End rpoly_ak1

    template <class T>
    void Fxshfr_ak1(int L2, int* NZ, T sr, T bnd, T K[MDP1], int N, T p[MDP1], int NN, T qp[MDP1], T* lzi, T* lzr, T* szi, T* szr){

    // Computes up to L2 fixed shift K-polynomials, testing for convergence in the linear or
    // quadratic case. Initiates one of the variable shift iterations and returns with the
    // number of zeros found.

    // L2 limit of fixed shift steps
    // NZ number of zeros found

    int fflag, i, iFlag, j, spass, stry, tFlag, vpass, vtry;
    T a, a1, a3, a7, b, betas, betav, c, d, e, f, g, h, oss, ots, otv, ovv, s, ss, ts, tss, tv, tvv, u, ui, v, vi, vv;
    T qk[MDP1], svk[MDP1];

    *NZ = 0;
    betav = betas = 0.25;
    u = -(2.0*sr);
    oss = sr;
    ovv = v = bnd;

    //Evaluate polynomial by synthetic division
    QuadSD_ak1(NN, u, v, p, qp, &a, &b);

    tFlag = calcSC_ak1(N, a, b, &a1, &a3, &a7, &c, &d, &e, &f, &g, &h, K, u, v, qk);

    for (j = 0; j < L2; j++){

        //Calculate next K polynomial and estimate v
        nextK_ak1(N, tFlag, a, b, a1, &a3, &a7, K, qk, qp);
        tFlag = calcSC_ak1(N, a, b, &a1, &a3, &a7, &c, &d, &e, &f, &g, &h, K, u, v, qk);
        newest_ak1(tFlag, &ui, &vi, a, a1, a3, a7, b, c, d, f, g, h, u, v, K, N, p);

        vv = vi;

        // Estimate s

        ss = ((K[N - 1] != 0.0) ? -(p[N]/K[N - 1]) : 0.0);

        ts = tv = 1.0;

        if ((j != 0) && (tFlag != 3)){

           // Compute relative measures of convergence of s and v sequences

            tv = ((vv != 0.0) ? fabs((vv - ovv)/vv) : tv);
            ts = ((ss != 0.0) ? fabs((ss - oss)/ss) : ts);

            // If decreasing, multiply the two most recent convergence measures

            tvv = ((tv < otv) ? tv*otv : 1.0);
            tss = ((ts < ots) ? ts*ots : 1.0);

            // Compare with convergence criteria

            vpass = ((tvv < betav) ? 1 : 0);
            spass = ((tss < betas) ? 1 : 0);

            if ((spass) || (vpass)){

                // At least one sequence has passed the convergence test.
                // Store variables before iterating

                for (i = 0; i < N; i++)   svk[i] = K[i];

                s = ss;

                // Choose iteration according to the fastest converging sequence

                stry = vtry = 0;
                fflag = 1;

                do {

                    iFlag = 1; // Begin each loop by assuming RealIT will be called UNLESS iFlag changed below

                    if ((fflag && ((fflag = 0) == 0)) && ((spass) && (!vpass || (tss < tvv)))){
                        ; // Do nothing. Provides a quick "short circuit".
                    } // End if (fflag)

                    else { // else !fflag
                        QuadIT_ak1(N, NZ, ui, vi, szr, szi, lzr, lzi, qp, NN, &a, &b, p, qk, &a1, &a3, &a7, &d, &e, &f, &g, &h, K);

                        if ((*NZ) > 0)   return;

                        // Quadratic iteration has failed. Flag that it has been tried and decrease the
                        // convergence criterion

                        vtry = 1;
                        betav *= 0.25;

                        // Try linear iteration if it has not been tried and the s sequence is converging
                        if (stry || (!spass)){
                            iFlag = 0;
                        } // End if (stry || (!spass))
                        else {
                            for (i = 0; i < N; i++)   K[i] = svk[i];
                        } // End if (stry || !spass)

                    } // End else !fflag

                    if (iFlag != 0){
                        RealIT_ak1(&iFlag, NZ, &s, N, p, NN, qp, szr, szi, K, qk);

                        if ((*NZ) > 0)   return;

                        // Linear iteration has failed. Flag that it has been tried and decrease the
                        // convergence criterion

                        stry = 1;
                        betas *= 0.25;

                        if (iFlag != 0){

                            // If linear iteration signals an almost T real zero, attempt quadratic iteration

                            ui = -(s + s);
                            vi = s*s;
                            continue;

                        } // End if (iFlag != 0)
                    } // End if (iFlag != 0)

                    // Restore variables
                    for (i = 0; i < N; i++)   K[i] = svk[i];

                    // Try quadratic iteration if it has not been tried and the v sequence is converging

                } while (vpass && !vtry); // End do-while loop

                // Re-compute qp and scalar values to continue the second stage

                QuadSD_ak1(NN, u, v, p, qp, &a, &b);
                tFlag = calcSC_ak1(N, a, b, &a1, &a3, &a7, &c, &d, &e, &f, &g, &h, K, u, v, qk);

            } // End if ((spass) || (vpass))

        } // End if ((j != 0) && (tFlag != 3))

        ovv = vv;
        oss = ss;
        otv = tv;
        ots = ts;
    } // End for j

    return;
    } // End Fxshfr_ak1

    template <class T>
    void QuadSD_ak1(int NN, T u, T v, T p[MDP1], T q[MDP1], T* a, T* b){

    // Divides p by the quadratic 1, u, v placing the quotient in q and the remainder in a, b

    int i;

    q[0] = *b = p[0];
    q[1] = *a = -((*b)*u) + p[1];

    for (i = 2; i < NN; i++){
        q[i] = -((*a)*u + (*b)*v) + p[i];
        *b = (*a);
        *a = q[i];
    } // End for i

    return;
    } // End QuadSD_ak1

    template <class T>
    int calcSC_ak1(int N, T a, T b, T* a1, T* a3, T* a7, T* c, T* d, T* e, T* f, T* g, T* h, T K[MDP1], T u, T v, T qk[MDP1]){

    // This routine calculates scalar quantities used to compute the next K polynomial and
    // new estimates of the quadratic coefficients.

    // calcSC - integer variable set here indicating how the calculations are normalized
    // to avoid overflow.

    int dumFlag = 3; // TYPE = 3 indicates the quadratic is almost a factor of K

    // Synthetic division of K by the quadratic 1, u, v
    QuadSD_ak1(N, u, v, K, qk, c, d);

    if (fabs((*c)) <= (100.0*DBL_EPSILON*fabs(K[N - 1]))) {
        if (fabs((*d)) <= (100.0*DBL_EPSILON*fabs(K[N - 2])))   return dumFlag;
    } // End if (fabs(c) <= (100.0*DBL_EPSILON*fabs(K[N - 1])))

    *h = v*b;
    if (fabs((*d)) >= fabs((*c))){
        dumFlag = 2; // TYPE = 2 indicates that all formulas are divided by d
        *e = a/(*d);
        *f = (*c)/(*d);
        *g = u*b;
        *a3 = (*e)*((*g) + a) + (*h)*(b/(*d));
        *a1 = -a + (*f)*b;
        *a7 = (*h) + ((*f) + u)*a;
    } // End if(fabs(d) >= fabs(c))
    else {
        dumFlag = 1; // TYPE = 1 indicates that all formulas are divided by c;
        *e = a/(*c);
        *f = (*d)/(*c);
        *g = (*e)*u;
        *a3 = (*e)*a + ((*g) + (*h)/(*c))*b;
        *a1 = -(a*((*d)/(*c))) + b;
        *a7 = (*g)*(*d) + (*h)*(*f) + a;
    } // End else

    return dumFlag;
    } // End calcSC_ak1


    template <class T>
    void nextK_ak1(int N, int tFlag, T a, T b, T a1, T* a3, T* a7, T K[MDP1], T qk[MDP1], T qp[MDP1]){

    // Computes the next K polynomials using the scalars computed in calcSC_ak1

    int i;
    T temp;

    if (tFlag == 3){ // Use unscaled form of the recurrence
        K[1] = K[0] = 0.0;

        for (i = 2; i < N; i++)   K[i] = qk[i - 2];

        return;
    } // End if (tFlag == 3)

    temp = ((tFlag == 1) ? b : a);

    if (fabs(a1) > (10.0*DBL_EPSILON*fabs(temp))){
        // Use scaled form of the recurrence

        (*a7) /= a1;
        (*a3) /= a1;
        K[0] = qp[0];
        K[1] = -((*a7)*qp[0]) + qp[1];

        for (i = 2; i < N; i++)   K[i] = -((*a7)*qp[i - 1]) + (*a3)*qk[i - 2] + qp[i];

    } // End if (fabs(a1) > (10.0*DBL_EPSILON*fabs(temp)))
    else {
        // If a1 is nearly zero, then use a special form of the recurrence

        K[0] = 0.0;
        K[1] = -(*a7)*qp[0];

        for (i = 2; i < N; i++)   K[i] = -((*a7)*qp[i - 1]) + (*a3)*qk[i - 2];
    } // End else

    return;

    } // End nextK_ak1

    template <class T>
    void newest_ak1(int tFlag, T* uu, T* vv, T a, T a1, T a3, T a7, T b, T c, T d, T f, T g, T h, T u, T v, T K[MDP1], int N, T p[MDP1]){
    // Compute new estimates of the quadratic coefficients using the scalars computed in calcSC_ak1

    T a4, a5, b1, b2, c1, c2, c3, c4, temp;

    (*vv) = (*uu) = 0.0; // The quadratic is zeroed

    if (tFlag != 3){

        if (tFlag != 2){
            a4 = a + u*b + h*f;
            a5 = c + (u + v*f)*d;
        } // End if (tFlag != 2)
        else { // else tFlag == 2
            a4 = (a + g)*f + h;
            a5 = (f + u)*c + v*d;
        } // End else tFlag == 2

        // Evaluate new quadratic coefficients

        b1 = -K[N - 1]/p[N];
        b2 = -(K[N - 2] + b1*p[N - 1])/p[N];
        c1 = v*b2*a1;
        c2 = b1*a7;
        c3 = b1*b1*a3;
        c4 = -(c2 + c3) + c1;
        temp = -c4 + a5 + b1*a4;
        if (temp != 0.0) {
            *uu= -((u*(c3 + c2) + v*(b1*a1 + b2*a7))/temp) + u;
            *vv = v*(1.0 + c4/temp);
        } // End if (temp != 0)

    } // End if (tFlag != 3)

    return;
    } // End newest_ak1

    template <class T>
    void QuadIT_ak1(int N, int* NZ, T uu, T vv, T* szr, T* szi, T* lzr, T* lzi, T qp[MDP1], int NN, T* a, T* b, T p[MDP1], T qk[MDP1], T* a1, T* a3, T* a7, T* d, T* e, T* f, T* g, T* h, T K[MDP1]){

    // Variable-shift K-polynomial iteration for a quadratic factor converges only if the
    // zeros are equimodular or nearly so.

    int i, j = 0, tFlag, triedFlag = 0;
    T c, ee, mp, omp, relstp, t, u, ui, v, vi, zm;

    *NZ = 0; // Number of zeros found
    u = uu; // uu and vv are coefficients of the starting quadratic
    v = vv;

    do {
        Quad_ak1(1.0, u, v, szr, szi, lzr, lzi);

        // Return if roots of the quadratic are real and not close to multiple or nearly
        // equal and of opposite sign.

        if (fabs(fabs(*szr) - fabs(*lzr)) > 0.01*fabs(*lzr))   break;

        // Evaluate polynomial by quadratic synthetic division

        QuadSD_ak1(NN, u, v, p, qp, a, b);

        mp = fabs(-((*szr)*(*b)) + (*a)) + fabs((*szi)*(*b));

        // Compute a rigorous bound on the rounding error in evaluating p

        zm = sqrt(fabs(v));
        ee = 2.0*fabs(qp[0]);
        t = -((*szr)*(*b));

        for (i = 1; i < N; i++)   ee = ee*zm + fabs(qp[i]);

        ee = ee*zm + fabs((*a) + t);
        ee = (9.0*ee + 2.0*fabs(t) - 7.0*(fabs((*a) + t) + zm*fabs((*b))))*DBL_EPSILON;

        // Iteration has converged sufficiently if the polynomial value is less than 20 times this bound

        if (mp <= 20.0*ee){
            *NZ = 2;
            break;
        } // End if (mp <= 20.0*ee)

        j++;

        // Stop iteration after 20 steps
        if (j > 20)   break;

        if (j >= 2){
            if ((relstp <= 0.01) && (mp >= omp) && (!triedFlag)){
            // A cluster appears to be stalling the convergence. Five fixed shift
            // steps are taken with a u, v close to the cluster.

            relstp = ((relstp < DBL_EPSILON) ? sqrt(DBL_EPSILON) : sqrt(relstp));

            u -= u*relstp;
            v += v*relstp;

            QuadSD_ak1(NN, u, v, p, qp, a, b);

            for (i = 0; i < 5; i++){
                tFlag = calcSC_ak1(N, *a, *b, a1, a3, a7, &c, d, e, f, g, h, K, u, v, qk);
                nextK_ak1(N, tFlag, *a, *b, *a1, a3, a7, K, qk, qp);
            } // End for i

            triedFlag = 1;
            j = 0;

            } // End if ((relstp <= 0.01) && (mp >= omp) && (!triedFlag))

        } // End if (j >= 2)

        omp = mp;

        // Calculate next K polynomial and new u and v

        tFlag = calcSC_ak1(N, *a, *b, a1, a3, a7, &c, d, e, f, g, h, K, u, v, qk);
        nextK_ak1(N, tFlag, *a, *b, *a1, a3, a7, K, qk, qp);
        tFlag = calcSC_ak1(N, *a, *b, a1, a3, a7, &c, d, e, f, g, h, K, u, v, qk);
        newest_ak1(tFlag, &ui, &vi, *a, *a1, *a3, *a7, *b, c, *d, *f, *g, *h, u, v, K, N, p);

        // If vi is zero, the iteration is not converging
        if (vi != 0){
            relstp = fabs((-v + vi)/vi);
            u = ui;
            v = vi;
        } // End if (vi != 0)
    } while (vi != 0); // End do-while loop

    return;

    } //End QuadIT_ak1

    template <class T>
    void RealIT_ak1(int* iFlag, int* NZ, T* sss, int N, T p[MDP1], int NN, T qp[MDP1], T* szr, T* szi, T K[MDP1], T qk[MDP1]){

    // Variable-shift H-polynomial iteration for a real zero

    // sss - starting iterate
    // NZ - number of zeros found
    // iFlag - flag to indicate a pair of zeros near real axis

    int i, j = 0, nm1 = N - 1;
    T ee, kv, mp, ms, omp, pv, s, t;

    *iFlag = *NZ = 0;
    s = *sss;

    for ( ; ; ) {
        qp[0] = pv = p[0];

        // Evaluate p at s
        for (i = 1; i < NN; i++)   qp[i] = pv = pv*s + p[i];

        mp = fabs(pv);

        // Compute a rigorous bound on the error in evaluating p

        ms = fabs(s);
        ee = 0.5*fabs(qp[0]);
        for (i = 1; i < NN; i++)   ee = ee*ms + fabs(qp[i]);

        // Iteration has converged sufficiently if the polynomial value is less than
        // 20 times this bound

        if (mp <= 20.0*DBL_EPSILON*(2.0*ee - mp)){
            *NZ = 1;
            *szr = s;
            *szi = 0.0;
            break;
        } // End if (mp <= 20.0*DBL_EPSILON*(2.0*ee - mp))

        j++;

        // Stop iteration after 10 steps

        if (j > 10)   break;

        if (j >= 2){
            if ((fabs(t) <= 0.001*fabs(-t + s)) && (mp > omp)){
                // A cluster of zeros near the real axis has been encountered;
                // Return with iFlag set to initiate a quadratic iteration

                *iFlag = 1;
                *sss = s;
                break;
            } // End if ((fabs(t) <= 0.001*fabs(s - t)) && (mp > omp))

        } //End if (j >= 2)

        // Return if the polynomial value has increased significantly

        omp = mp;

        // Compute t, the next polynomial and the new iterate
        qk[0] = kv = K[0];
        for (i = 1; i < N; i++)   qk[i] = kv = kv*s + K[i];

        if (fabs(kv) > fabs(K[nm1])*10.0*DBL_EPSILON){
            // Use the scaled form of the recurrence if the value of K at s is non-zero
            t = -(pv/kv);
            K[0] = qp[0];
            for (i = 1; i < N; i++)   K[i] = t*qk[i - 1] + qp[i];
        } // End if (fabs(kv) > fabs(K[nm1])*10.0*DBL_EPSILON)
        else { // else (fabs(kv) <= fabs(K[nm1])*10.0*DBL_EPSILON)
            // Use unscaled form
            K[0] = 0.0;
            for (i = 1; i < N; i++)   K[i] = qk[i - 1];
        } // End else (fabs(kv) <= fabs(K[nm1])*10.0*DBL_EPSILON)

        kv = K[0];
        for (i = 1; i < N; i++)   kv = kv*s + K[i];

        t = ((fabs(kv) > (fabs(K[nm1])*10.0*DBL_EPSILON)) ? -(pv/kv) : 0.0);

        s += t;

    } // End infinite for loop

    return;

    } // End RealIT_ak1

    template <class T>
    void Quad_ak1(T a, T b1, T c, T* sr, T* si, T* lr, T* li) {
    // Calculates the zeros of the quadratic a*Z^2 + b1*Z + c
    // The quadratic formula, modified to avoid overflow, is used to find the larger zero if the
    // zeros are real and both zeros are complex. The smaller real zero is found directly from
    // the product of the zeros c/a.

        T b, d, e;

        *sr = *si = *lr = *li = 0.0;

        if (a == 0) {
            *sr = ((b1 != 0) ? -(c/b1) : *sr);
            return;
        } // End if (a == 0))

        if (c == 0){
            *lr = -(b1/a);
            return;
        } // End if (c == 0)

        // Compute discriminant avoiding overflow

        b = b1/2.0;
        if (fabs(b) < fabs(c)){
            e = ((c >= 0) ? a : -a);
            e = -e + b*(b/fabs(c));
            d = sqrt(fabs(e))*sqrt(fabs(c));
        } // End if (fabs(b) < fabs(c))
        else { // Else (fabs(b) >= fabs(c))
            e = -((a/b)*(c/b)) + 1.0;
            d = sqrt(fabs(e))*(fabs(b));
        } // End else (fabs(b) >= fabs(c))

        if (e >= 0) {
            // Real zeros

            d = ((b >= 0) ? -d : d);
            *lr = (-b + d)/a;
            *sr = ((*lr != 0) ? (c/(*lr))/a : *sr);
        } // End if (e >= 0)
        else { // Else (e < 0)
            // Complex conjugate zeros

            *lr = *sr = -(b/a);
            *si = fabs(d/a);
            *li = -(*si);
        } // End else (e < 0)

        return;
    } // End Quad_ak1

}//end namespace std

