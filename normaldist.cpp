#if defined( _DEBUG ) || defined( BOUNDS_CHECK_STL_ARRAYS )
#define _SECURE_SCL 1
#define _SECURE_SCL_THROWS 1
#define _SCL_SECURE_NO_WARNINGS
#define _HAS_ITERATOR_DEBUGGING 0
#else
#define _SECURE_SCL 0
#endif
#if defined(_MSC_VER)
# define NOMINMAX // to suppress MSVC's definitions of min() and max()
// These four pragmas are the equivalent to /fp:fast.
# pragma float_control( except, off )
# pragma float_control( precise, off )
# pragma fp_contract( on )
# pragma fenv_access( off )
#endif

#include "normaldist.h"
#include <cfloat>
#include <cmath>

namespace {
   const double norm_cdf_asymptotic_expansion_first_threshold = -10.0;
   const double norm_cdf_asymptotic_expansion_second_threshold = -1/sqrt(DBL_EPSILON);
}

namespace {
   inline double d_int(const double x){ return( (x>0) ? floor(x) : -floor(-x) ); }
}

normaldist::normaldist()
{
}

double normaldist::calerf(double x, const int jint)
{
    static const double a[5] = { 3.1611237438705656,113.864154151050156,377.485237685302021,3209.37758913846947,.185777706184603153 };
    static const double b[4] = { 23.6012909523441209,244.024637934444173,1282.61652607737228,2844.23683343917062 };
    static const double c__[9] = { .564188496988670089,8.88314979438837594,66.1191906371416295,298.635138197400131,881.95222124176909,1712.04761263407058,2051.07837782607147,1230.33935479799725,2.15311535474403846e-8 };
    static const double d__[8] = { 15.7449261107098347,117.693950891312499,537.181101862009858,1621.38957456669019,3290.79923573345963,4362.61909014324716,3439.36767414372164,1230.33935480374942 };
    static const double p[6] = { .305326634961232344,.360344899949804439,.125781726111229246,.0160837851487422766,6.58749161529837803e-4,.0163153871373020978 };
    static const double q[5] = { 2.56852019228982242,1.87295284992346047,.527905102951428412,.0605183413124413191,.00233520497626869185 };

    static const double zero = 0.;
    static const double half = .5;
    static const double one = 1.;
    static const double two = 2.;
    static const double four = 4.;
    static const double sqrpi = 0.56418958354775628695;
    static const double thresh = .46875;
    static const double sixten = 16.;

    double y, del, ysq, xden, xnum, result;

    static const double xinf = 1.79e308;
    static const double xneg = -26.628;
    static const double xsmall = 1.11e-16;
    static const double xbig = 26.543;
    static const double xhuge = 6.71e7;
    static const double xmax = 2.53e307;

    y = fabs(x);

    if (y <= thresh) {
       /* ------------------------------------------------------------------ */
       /*  Evaluate  erf  for  |X| <= 0.46875 */
       /* ------------------------------------------------------------------ */
       /*<             YSQ = ZERO >*/
       ysq = zero;
       /*<             IF (Y .GT. XSMALL) YSQ = Y * Y >*/
       if (y > xsmall) {
          ysq = y * y;
       }
       /*<             XNUM = A(5)*YSQ >*/
       xnum = a[4] * ysq;
       /*<             XDEN = YSQ >*/
       xden = ysq;
       /*<             DO 20 I = 1, 3 >*/
       for (int i__ = 1; i__ <= 3; ++i__) {
          /*<                XNUM = (XNUM + A(I)) * YSQ >*/
          xnum = (xnum + a[i__ - 1]) * ysq;
          /*<                XDEN = (XDEN + B(I)) * YSQ >*/
          xden = (xden + b[i__ - 1]) * ysq;
          /*<    20       CONTINUE >*/
          /* L20: */
       }
       /*<             RESULT = X * (XNUM + A(4)) / (XDEN + B(4)) >*/
       result = x * (xnum + a[3]) / (xden + b[3]);
       /*<             IF (JINT .NE. 0) RESULT = ONE - RESULT >*/
       if (jint != 0) {
          result = one - result;
       }
       /*<             IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT >*/
       if (jint == 2) {
          result = exp(ysq) * result;
       }
       /*<             GO TO 800 >*/
       goto L800;
       /* ------------------------------------------------------------------ */
       /*  Evaluate  erfc  for 0.46875 <= |X| <= 4.0 */
       /* ------------------------------------------------------------------ */
       /*<          ELSE IF (Y .LE. FOUR) THEN >*/
    } else if (y <= four) {
       /*<             XNUM = C(9)*Y >*/
       xnum = c__[8] * y;
       /*<             XDEN = Y >*/
       xden = y;
       /*<             DO 120 I = 1, 7 >*/
       for (int i__ = 1; i__ <= 7; ++i__) {
          /*<                XNUM = (XNUM + C(I)) * Y >*/
          xnum = (xnum + c__[i__ - 1]) * y;
          /*<                XDEN = (XDEN + D(I)) * Y >*/
          xden = (xden + d__[i__ - 1]) * y;
          /*<   120       CONTINUE >*/
          /* L120: */
       }
       /*<             RESULT = (XNUM + C(8)) / (XDEN + D(8)) >*/
       result = (xnum + c__[7]) / (xden + d__[7]);
       /*<             IF (JINT .NE. 2) THEN >*/
       if (jint != 2) {
          /*<                YSQ = AINT(Y*SIXTEN)/SIXTEN >*/
          double d__1 = y * sixten;
          ysq = d_int(d__1) / sixten;
          /*<                DEL = (Y-YSQ)*(Y+YSQ) >*/
          del = (y - ysq) * (y + ysq);
          /*<                RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT >*/
          d__1 = exp(-ysq * ysq) * exp(-del);
          result = d__1 * result;
          /*<             END IF >*/
       }
       /* ------------------------------------------------------------------ */
       /*  Evaluate  erfc  for |X| > 4.0 */
       /* ------------------------------------------------------------------ */
       /*<          ELSE >*/
    } else {
       /*<             RESULT = ZERO >*/
       result = zero;
       /*<             IF (Y .GE. XBIG) THEN >*/
       if (y >= xbig) {
          /*<                IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300 >*/
          if (jint != 2 || y >= xmax) {
             goto L300;
          }
          /*<                IF (Y .GE. XHUGE) THEN >*/
          if (y >= xhuge) {
             /*<                   RESULT = SQRPI / Y >*/
             result = sqrpi / y;
             /*<                   GO TO 300 >*/
             goto L300;
             /*<                END IF >*/
          }
          /*<             END IF >*/
       }
       /*<             YSQ = ONE / (Y * Y) >*/
       ysq = one / (y * y);
       /*<             XNUM = P(6)*YSQ >*/
       xnum = p[5] * ysq;
       /*<             XDEN = YSQ >*/
       xden = ysq;
       /*<             DO 240 I = 1, 4 >*/
       for (int i__ = 1; i__ <= 4; ++i__) {
          /*<                XNUM = (XNUM + P(I)) * YSQ >*/
          xnum = (xnum + p[i__ - 1]) * ysq;
          /*<                XDEN = (XDEN + Q(I)) * YSQ >*/
          xden = (xden + q[i__ - 1]) * ysq;
          /*<   240       CONTINUE >*/
          /* L240: */
       }
       /*<             RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5)) >*/
       result = ysq * (xnum + p[4]) / (xden + q[4]);
       /*<             RESULT = (SQRPI -  RESULT) / Y >*/
       result = (sqrpi - result) / y;
       /*<             IF (JINT .NE. 2) THEN >*/
       if (jint != 2) {
          /*<                YSQ = AINT(Y*SIXTEN)/SIXTEN >*/
          double d__1 = y * sixten;
          ysq = d_int(d__1) / sixten;
          /*<                DEL = (Y-YSQ)*(Y+YSQ) >*/
          del = (y - ysq) * (y + ysq);
          /*<                RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT >*/
          d__1 = exp(-ysq * ysq) * exp(-del);
          result = d__1 * result;
          /*<             END IF >*/
       }
       /*<       END IF >*/
    }
    /* ------------------------------------------------------------------ */
    /*  Fix up for negative argument, erf, etc. */
    /* ------------------------------------------------------------------ */
    /*<   300 IF (JINT .EQ. 0) THEN >*/
 L300:
    if (jint == 0) {
       /*<             RESULT = (HALF - RESULT) + HALF >*/
       result = (half - result) + half;
       /*<             IF (X .LT. ZERO) RESULT = -RESULT >*/
       if (x < zero) {
          result = -(result);
       }
       /*<          ELSE IF (JINT .EQ. 1) THEN >*/
    } else if (jint == 1) {
       /*<             IF (X .LT. ZERO) RESULT = TWO - RESULT >*/
       if (x < zero) {
          result = two - result;
       }
       /*<          ELSE >*/
    } else {
       /*<             IF (X .LT. ZERO) THEN >*/
       if (x < zero) {
          /*<                IF (X .LT. XNEG) THEN >*/
          if (x < xneg) {
             /*<                      RESULT = XINF >*/
             result = xinf;
             /*<                   ELSE >*/
          } else {
             /*<                      YSQ = AINT(X*SIXTEN)/SIXTEN >*/
             double d__1 = x * sixten;
             ysq = d_int(d__1) / sixten;
             /*<                      DEL = (X-YSQ)*(X+YSQ) >*/
             del = (x - ysq) * (x + ysq);
             /*<                      Y = EXP(YSQ*YSQ) * EXP(DEL) >*/
             y = exp(ysq * ysq) * exp(del);
             /*<                      RESULT = (Y+Y) - RESULT >*/
             result = y + y - result;
             /*<                END IF >*/
          }
          /*<             END IF >*/
       }
       /*<       END IF >*/
    }
    /*<   800 RETURN >*/
 L800:
    return result;
}

double normaldist::erf_cody(double x)
{
    return calerf(x, 0);
}

double normaldist::erfc_cody(double x)
{
    return calerf(x, 1);
}

double normaldist::erfcx_cody(double x)
{
    return calerf(x, 2);
}

double normaldist::norm_cdf(double z)
{
    if (z <= norm_cdf_asymptotic_expansion_first_threshold) {
       // Asymptotic expansion for very negative z following (26.2.12) on page 408
       // in M. Abramowitz and A. Stegun, Pocketbook of Mathematical Functions, ISBN 3-87144818-4.
       double sum = 1;
       if (z >= norm_cdf_asymptotic_expansion_second_threshold) {
          double zsqr = z * z, i = 1, g = 1, x, y, a = DBL_MAX, lasta;
          do {
             lasta = a;
             x = (4 * i - 3) / zsqr;
             y = x * ((4 * i - 1) / zsqr);
             a = g * (x - y);
             sum -= a;
             g *= y;
             ++i;
             a = fabs(a);
          } while (lasta > a && a >= fabs(sum * DBL_EPSILON));
       }
       return -norm_pdf(z) * sum / z;
    }
    return 0.5*erfc_cody( -z*ONE_OVER_SQRT_TWO );
}

double normaldist::inverse_norm_cdf(double u)
{
    //
       // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
       //
       // Produces the normal deviate Z corresponding to a given lower
       // tail area of u; Z is accurate to about 1 part in 10**16.
       // see http://lib.stat.cmu.edu/apstat/241
       //
       const double split1 = 0.425;
       const double split2 = 5.0;
       const double const1 = 0.180625;
       const double const2 = 1.6;

       // Coefficients for P close to 0.5
       const double A0 = 3.3871328727963666080E0;
       const double A1 = 1.3314166789178437745E+2;
       const double A2 = 1.9715909503065514427E+3;
       const double A3 = 1.3731693765509461125E+4;
       const double A4 = 4.5921953931549871457E+4;
       const double A5 = 6.7265770927008700853E+4;
       const double A6 = 3.3430575583588128105E+4;
       const double A7 = 2.5090809287301226727E+3;
       const double B1 = 4.2313330701600911252E+1;
       const double B2 = 6.8718700749205790830E+2;
       const double B3 = 5.3941960214247511077E+3;
       const double B4 = 2.1213794301586595867E+4;
       const double B5 = 3.9307895800092710610E+4;
       const double B6 = 2.8729085735721942674E+4;
       const double B7 = 5.2264952788528545610E+3;
       // Coefficients for P not close to 0, 0.5 or 1.
       const double C0 = 1.42343711074968357734E0;
       const double C1 = 4.63033784615654529590E0;
       const double C2 = 5.76949722146069140550E0;
       const double C3 = 3.64784832476320460504E0;
       const double C4 = 1.27045825245236838258E0;
       const double C5 = 2.41780725177450611770E-1;
       const double C6 = 2.27238449892691845833E-2;
       const double C7 = 7.74545014278341407640E-4;
       const double D1 = 2.05319162663775882187E0;
       const double D2 = 1.67638483018380384940E0;
       const double D3 = 6.89767334985100004550E-1;
       const double D4 = 1.48103976427480074590E-1;
       const double D5 = 1.51986665636164571966E-2;
       const double D6 = 5.47593808499534494600E-4;
       const double D7 = 1.05075007164441684324E-9;
       // Coefficients for P very close to 0 or 1
       const double E0 = 6.65790464350110377720E0;
       const double E1 = 5.46378491116411436990E0;
       const double E2 = 1.78482653991729133580E0;
       const double E3 = 2.96560571828504891230E-1;
       const double E4 = 2.65321895265761230930E-2;
       const double E5 = 1.24266094738807843860E-3;
       const double E6 = 2.71155556874348757815E-5;
       const double E7 = 2.01033439929228813265E-7;
       const double F1 = 5.99832206555887937690E-1;
       const double F2 = 1.36929880922735805310E-1;
       const double F3 = 1.48753612908506148525E-2;
       const double F4 = 7.86869131145613259100E-4;
       const double F5 = 1.84631831751005468180E-5;
       const double F6 = 1.42151175831644588870E-7;
       const double F7 = 2.04426310338993978564E-15;

       if (u<=0)
          return log(u);
       if (u>=1)
          return log(1-u);

       const double q = u-0.5;
       if (fabs(q) <= split1)
       {
          const double r = const1 - q*q;
          return q * (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0) /
             (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + 1.0);
       }
       else
       {
          double r = q<0.0 ? u : 1.0-u;
          r = sqrt(-log(r));
          double ret;
          if (r < split2)
          {
             r = r - const2;
             ret = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0) /
                (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + 1.0);
          }
          else
          {
             r = r - split2;
             ret = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0) /
                (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + 1.0);
          }
          return q<0.0 ? -ret : ret;
       }
}
