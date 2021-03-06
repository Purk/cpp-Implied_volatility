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
#include <cmath>
#include <cfloat>

namespace {
   inline double d_int(const double x){ return( (x>0) ? floor(x) : -floor(-x) ); }
}

/*<       SUBROUTINE CALERF(ARG,RESULT,JINT) >*/
double calerf(double x, const int jint) {

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

double erf_cody(double x){

   return calerf(x, 0);
}

double erfc_cody(double x) {

   return calerf(x, 1);

}

double erfcx_cody(double x) {

   return calerf(x, 2);

}
