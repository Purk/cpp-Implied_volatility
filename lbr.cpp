#include "lbr.h"
#include "normaldist.h"
#include "rationalcubic.h"
#include <cfloat>
#include <cmath>
#include <algorithm>

#define TWO_PI                        6.283185307179586476925286766559005768394338798750
#define SQRT_PI_OVER_TWO              1.253314137315500251207882642405522626503493370305  // sqrt(pi/2) to avoid misinterpretation.
#define SQRT_THREE                    1.732050807568877293527446341505872366942805253810
#define SQRT_ONE_OVER_THREE           0.577350269189625764509148780501957455647601751270
#define TWO_PI_OVER_SQRT_TWENTY_SEVEN 1.209199576156145233729385505094770488189377498728 // 2*pi/sqrt(27)
#define PI_OVER_SIX                   0.523598775598298873077107230546583814032861566563

normaldist normalDistObj;

namespace {
   static const double SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);
   static const double FOURTH_ROOT_DBL_EPSILON = sqrt(SQRT_DBL_EPSILON);
   static const double EIGHTH_ROOT_DBL_EPSILON = sqrt(FOURTH_ROOT_DBL_EPSILON);
   static const double SIXTEENTH_ROOT_DBL_EPSILON = sqrt(EIGHTH_ROOT_DBL_EPSILON);
   static const double SQRT_DBL_MIN = sqrt(DBL_MIN);
   static const double SQRT_DBL_MAX = sqrt(DBL_MAX);

   // Set this to 0 if you want positive results for (positive) denormalised inputs, else to DBL_MIN.
   // Note that you cannot achieve full machine accuracy from denormalised inputs!
   static const double DENORMALISATION_CUTOFF = 0;

   static const double VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC = -DBL_MAX;
   static const double VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM = DBL_MAX;

   inline bool is_below_horizon(double x){ return fabs(x) < DENORMALISATION_CUTOFF; } // This weeds out denormalised (a.k.a. 'subnormal') numbers.

   typedef struct {
       #if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
             long data;
       #else
             int data;
       #endif
   }atomic_t;

   static atomic_t implied_volatility_maximum_iterations = { 2 }; // (DBL_DIG*20)/3 ≈ 100 . Only needed when the iteration effectively alternates Householder/Halley/Newton steps and binary nesting due to roundoff truncation.

   #ifdef ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT
   static atomic_t implied_volatility_output_type = { 0 };
   inline double implied_volatility_output(int count, double volatility)
   { return implied_volatility_output_type.data>0 ? count : volatility; }
   #else
   inline double implied_volatility_output(int count, double volatility){ return volatility; }
   #endif

    #ifdef ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER
       static atomic_t implied_volatility_householder_method_order = { 4 };
       inline double householder_factor(double newton, double halley, double hh3){
          return implied_volatility_householder_method_order.data > 3 ?
                      (1+0.5*halley*newton)/(1+newton*(halley+hh3*newton/6)) :
                      ( implied_volatility_householder_method_order.data > 2 ? 1/(1+0.5*halley*newton) : 1 );
       }
    #else
       inline double householder_factor(double newton, double halley, double hh3){ return (1+0.5*halley*newton)/(1+newton*(halley+hh3*newton/6)); }
    #endif

}

lbr::lbr()
{

}

double lbr::set_implied_volatility_maximum_iterations(double t)
{
    int i = (int)t;
    if (i>=0) {
    #if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
        InterlockedExchange(&(implied_volatility_maximum_iterations.data),i);
    #elif defined( __x86__ ) || defined( __x86_64__ )
        implied_volatility_maximum_iterations.data = i;
    #else
        # error Atomic operations not implemented for this platform.
    #endif
    }

    return implied_volatility_maximum_iterations.data;
}

#ifdef ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT
double lbr::set_implied_volatility_output_type(double t)
{
    int i = (int)t;
    #if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
       InterlockedExchange(&(implied_volatility_output_type.data),i);
    #elif defined( __x86__ ) || defined( __x86_64__ )
       implied_volatility_output_type.data = i;
    #else
    # error Atomic operations not implemented for this platform.
    #endif
       return implied_volatility_output_type.data;
}
#endif

#ifdef ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER
double lbr::set_implied_volatility_householder_method_order(double t)
{
    int i = (int)t;
       if (i>=0) {
    #if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
          InterlockedExchange(&(implied_volatility_householder_method_order.data),i);
    #elif defined( __x86__ ) || defined( __x86_64__ )
          implied_volatility_householder_method_order.data = i;
    #else
    # error Atomic operations not implemented for this platform.
    #endif
       }
       return implied_volatility_householder_method_order.data;
}
#endif

double normalised_intrinsic(double x, double q /* q=±1 */)
{
   if (q*x<=0)
      return 0;
   const double x2=x*x;
   if (x2<98*FOURTH_ROOT_DBL_EPSILON ) {
      return fabs( std::max( (q<0?-1:1)*x*(1+x2*((1.0/24.0)+x2*((1.0/1920.0)+x2*
                            ((1.0/322560.0)+(1.0/92897280.0)*x2)))) , 0.0 ) );
   }

   const double b_max = exp(0.5*x), one_over_b_max = 1 / b_max;
   return fabs(std::max((q<0?-1:1)*(b_max-one_over_b_max),0.));
}

double normalised_intrinsic_call(double x)
{
    return normalised_intrinsic(x,1);
}

double asymptotic_expansion_of_normalised_black_call(double h, double t){
   const double e=(t/h)*(t/h), r=((h+t)*(h-t)), q=(h/r)*(h/r);

   const double asymptotic_expansion_sum = (2.0+q*(-6.0E0-2.0*e+3.0*q*(1.0E1+e*(2.0E1+2.0*e)+
                                            5.0*q*(-1.4E1+e*(-7.0E1+e*(-4.2E1-2.0*e))+7.0*q*
                                            (1.8E1+e*(1.68E2+e*(2.52E2+e*(7.2E1+2.0*e)))+9.0*q*
                                            (-2.2E1+e*(-3.3E2+e*(-9.24E2+e*(-6.6E2+e*(-1.1E2-2.0*e))))+
                                            1.1E1*q*(2.6E1+e*(5.72E2+e*(2.574E3+e*(3.432E3+e*(1.43E3+e*
                                            (1.56E2+2.0*e)))))+1.3E1*q*(-3.0E1+e*(-9.1E2+e*(-6.006E3+e*
                                            (-1.287E4+e*(-1.001E4+e*(-2.73E3+e*(-2.1E2-2.0*e))))))+
                                            1.5E1*q*(3.4E1+e*(1.36E3+e*(1.2376E4+e*(3.8896E4+e*(4.862E4+
                                            e*(2.4752E4+e*(4.76E3+e*(2.72E2+2.0*e)))))))+1.7E1*q*
                                            (-3.8E1+e*(-1.938E3+e*(-2.3256E4+e*(-1.00776E5+e*
                                            (-1.84756E5+e*(-1.51164E5+e*(-5.4264E4+e*(-7.752E3+e*
                                            (-3.42E2-2.0*e))))))))+1.9E1*q*(4.2E1+e*(2.66E3+e*
                                            (4.0698E4+e*(2.3256E5+e*(5.8786E5+e*(7.05432E5+e*
                                            (4.0698E5+e*(1.08528E5+e*(1.197E4+e*(4.2E2+2.0*e)))))))))+
                                            2.1E1*q*(-4.6E1+e*(-3.542E3+e*(-6.7298E4+e*(-4.90314E5+e*
                                            (-1.63438E6+e*(-2.704156E6+e*(-2.288132E6+e*(-9.80628E5+e*
                                            (-2.01894E5+e*(-1.771E4+e*(-5.06E2-2.0*e))))))))))+2.3E1*q*
                                            (5.0E1+e*(4.6E3+e*(1.0626E5+e*(9.614E5+e*(4.08595E6+e*
                                            (8.9148E6+e*(1.04006E7+e*(6.53752E6+e*(2.16315E6+e*
                                            (3.542E5+e*(2.53E4+e*(6.0E2+2.0*e)))))))))))+2.5E1*q*
                                            (-5.4E1+e*(-5.85E3+e*(-1.6146E5+e*(-1.77606E6+e*
                                            (-9.37365E6+e*(-2.607579E7+e*(-4.01166E7+e*(-3.476772E7+
                                            e*(-1.687257E7+e*(-4.44015E6+e*(-5.9202E5+e*(-3.51E4+e*
                                            (-7.02E2-2.0*e))))))))))))+2.7E1*q*(5.8E1+e*(7.308E3+e*
                                            (2.3751E5+e*(3.12156E6+e*(2.003001E7+e*(6.919458E7+e*
                                            (1.3572783E8+e*(1.5511752E8+e*(1.0379187E8+e*(4.006002E7+e*
                                            (8.58429E6+e*(9.5004E5+e*(4.7502E4+e*
                                            (8.12E2+2.0*e)))))))))))))+2.9E1*q*(-6.2E1+e*(-8.99E3+e*
                                            (-3.39822E5+e*(-5.25915E6+e*(-4.032015E7+e*(-1.6934463E8+e*
                                            (-4.1250615E8+e*(-6.0108039E8+e*(-5.3036505E8+e*
                                            (-2.8224105E8+e*(-8.870433E7+e*(-1.577745E7+e*
                                            (-1.472562E6+e*(-6.293E4+e*(-9.3E2-2.0*e))))))))))))))+
                                            3.1E1*q*(6.6E1+e*(1.0912E4+e*(4.74672E5+e*(8.544096E6+e*
                                            (7.71342E7+e*(3.8707344E8+e*(1.14633288E9+e*(2.07431664E9+
                                            e*(2.33360622E9+e*(1.6376184E9+e*(7.0963464E8+e*(1.8512208E8+
                                            e*(2.7768312E7+e*(2.215136E6+e*(8.184E4+e*(1.056E3+2.0*e)))))))))))))))+
                                            3.3E1*(-7.0E1+e*(-1.309E4+e*(-6.49264E5+e*(-1.344904E7+e*
                                            (-1.4121492E8+e*(-8.344518E8+e*(-2.9526756E9+e*(-6.49588632E9+e*
                                            (-9.0751353E9+e*(-8.1198579E9+e*(-4.6399188E9+e*(-1.6689036E9+e*
                                            (-3.67158792E8+e*(-4.707164E7+e*(-3.24632E6+e*(-1.0472E5+e*
                                            (-1.19E3-2.0*e)))))))))))))))))*q)))))))))))))))));

   const double b = ONE_OVER_SQRT_TWO_PI*exp((-0.5*(h*h+t*t)))*(t/r)*asymptotic_expansion_sum;
   return fabs(std::max(b , 0.0));
}

namespace { /* η */ static const double asymptotic_expansion_accuracy_threshold = -10; }

double normalised_black_call_using_erfcx(double h, double t)
{
    const double b = 0.5 * exp(-0.5*(h*h+t*t)) * (normalDistObj.erfcx_cody(-ONE_OVER_SQRT_TWO*(h+t)) -
                                                  normalDistObj.erfcx_cody(-ONE_OVER_SQRT_TWO*(h-t)) );
    return fabs(std::max(b,0.0));
}

double small_t_expansion_of_normalised_black_call(double h, double t){

   const double a = 1+h*(0.5*SQRT_TWO_PI)*normalDistObj.erfcx_cody(-ONE_OVER_SQRT_TWO*h), w=t*t, h2=h*h;
   const double expansion = 2*t*(a+w*((-1+3*a+a*h2)/6+w*((-7+15*a+h2*(-1+10*a+a*h2))/120+w*
                            ((-57+105*a+h2*(-18+105*a+h2*(-1+21*a+a*h2)))/5040+w*((-561+945*a+h2*
                            (-285+1260*a+h2*(-33+378*a+h2*(-1+36*a+a*h2))))/362880+w*((-6555+10395*
                            a+h2*(-4680+17325*a+h2*(-840+6930*a+h2*(-52+990*a+h2*(-1+55*a+a*h2)))))/39916800+
                            ((-89055+135135*a+h2*(-82845+270270*a+h2*(-20370+135135*a+h2*(-1926+25740*a+h2*
                            (-75+2145*a+h2*(-1+78*a+a*h2))))))*w)/6227020800.0))))));

   const double b = ONE_OVER_SQRT_TWO_PI*exp((-0.5*(h*h+t*t)))*expansion;
   return fabs(std::max(b,0.0));
}

namespace {/*τ*/static const double small_t_expansion_of_normalised_black_threshold = 2*SIXTEENTH_ROOT_DBL_EPSILON;}

double normalised_black_call_using_norm_cdf(double x, double s){
   const double h = x/s, t = 0.5*s, b_max = exp(0.5*x), b = normalDistObj.norm_cdf(h + t) * b_max - normalDistObj.norm_cdf(h - t) / b_max;
   return fabs(std::max(b,0.0));
}

double normalised_black_call_with_optimal_use_of_codys_functions(double x, double s){
   const double codys_threshold = 0.46875, h = x/s, t = 0.5*s, q1 = -ONE_OVER_SQRT_TWO*(h+t), q2 = -ONE_OVER_SQRT_TWO*(h-t);
   double two_b;
   if ( q1 < codys_threshold ) {
       if ( q2 < codys_threshold ) {
           two_b = exp(0.5*x)*normalDistObj.erfc_cody(q1) - exp(-0.5*x)*normalDistObj.erfc_cody(q2);
       }else {
           two_b = exp(0.5*x)*normalDistObj.erfc_cody(q1) - exp(-0.5*(h*h+t*t))*normalDistObj.erfcx_cody(q2);
       }
   }else {
       if ( q2 < codys_threshold ) {
           two_b =  exp(-0.5*(h*h+t*t))*normalDistObj.erfcx_cody(q1) - exp(-0.5*x)*normalDistObj.erfc_cody(q2);
       }else {
           two_b =  exp(-0.5*(h*h+t*t)) * ( normalDistObj.erfcx_cody(q1) - normalDistObj.erfcx_cody(q2) );
       }
   }
   return fabs(std::max(0.5*two_b,0.0));
}

double lbr::normalised_black_call(double x, double s)
{
    if (x>0)
       return normalised_intrinsic_call(x)+normalised_black_call(-x,s); // In the money.

    if (s<=fabs(x)*DENORMALISATION_CUTOFF)
       return normalised_intrinsic_call(x); // sigma=0 -> intrinsic value.

    if ( x < s*asymptotic_expansion_accuracy_threshold  &&  0.5*s*s+x <
         s*(small_t_expansion_of_normalised_black_threshold+asymptotic_expansion_accuracy_threshold) )
       return asymptotic_expansion_of_normalised_black_call(x/s,0.5*s);

    if ( 0.5*s < small_t_expansion_of_normalised_black_threshold )
       return small_t_expansion_of_normalised_black_call(x/s,0.5*s);

 #ifdef DO_NOT_OPTIMISE_NORMALISED_BLACK_IN_REGIONS_3_AND_4_FOR_CODYS_FUNCTIONS
    // When b is more than, say, about 85% of b_max=exp(x/2), then b is dominated by the first of the two terms in the Black formula, and we retain more accuracy by not attempting to combine the two terms in any way.
    // We evaluate the condition h+t>0.85  avoiding any divisions by s.
    if ( x+0.5*s*s > s*0.85 )
       return normalised_black_call_using_norm_cdf(x,s);
    return normalised_black_call_using_erfcx(x/s,0.5*s);
 #else
    return normalised_black_call_with_optimal_use_of_codys_functions(x,s);
 #endif
}

inline double square(double x){ return x*x; }

double lbr::normalised_vega(double x, double s)
{
    const double ax = fabs(x);
    return (ax<=0) ? ONE_OVER_SQRT_TWO_PI*exp(-0.125*s*s) : ( (s<=0 || s<=ax*SQRT_DBL_MIN) ? 0 :
                                    ONE_OVER_SQRT_TWO_PI*exp(-0.5*(square(x/s)+square(0.5*s))) );
}

double lbr::normalised_black(double x, double s, double q)
{
    return normalised_black_call(q<0?-x:x,s); /* Reciprocal-strike call-put equivalence */
}

double lbr::black(double F, double K, double sigma, double T, double q)
{
    const double intrinsic = fabs(std::max((q<0?K-F:F-K),0.0));
    // Map in-the-money to out-of-the-money
    if (q*(F-K)>0)
       return intrinsic + black(F,K,sigma,T,-q);

    return std::max(intrinsic,(sqrt(F)*sqrt(K))*normalised_black(log(F/K),sigma*sqrt(T),q));
}

#ifdef COMPUTE_LOWER_MAP_DERIVATIVES_INDIVIDUALLY
double f_lower_map(const double x,const double s){
   if (is_below_horizon(x))
      return 0;
   if (is_below_horizon(s))
      return 0;
   const double z=SQRT_ONE_OVER_THREE*fabs(x)/s, Phi=norm_cdf(-z);
   return TWO_PI_OVER_SQRT_TWENTY_SEVEN*fabs(x)*(Phi*Phi*Phi);
}

double d_f_lower_map_d_beta(const double x,const double s){
   if (is_below_horizon(s))
      return 1;
   const double z=SQRT_ONE_OVER_THREE*fabs(x)/s, y = z*z, Phi=norm_cdf(-z);
   return TWO_PI*y*(Phi*Phi) * exp(y+0.125*s*s);
}
double d2_f_lower_map_d_beta2(const double x,const double s){
   const double ax=fabs(x), z=SQRT_ONE_OVER_THREE*ax/s, y = z*z, s2=s*s, Phi=norm_cdf(-z), phi=norm_pdf(z);
   return PI_OVER_SIX * y/(s2*s) * Phi * ( 8*SQRT_THREE*s*ax + (3*s2*(s2-8)-8*x*x)*Phi/phi ) * exp(2*y+0.25*s2);
}
void compute_f_lower_map_and_first_two_derivatives(const double x,const double s,double &f,double &fp,double &fpp){
   f   = f_lower_map(x,s);
   fp  = d_f_lower_map_d_beta(x,s);
   fpp = d2_f_lower_map_d_beta2(x,s);
}
#else
void compute_f_lower_map_and_first_two_derivatives(const double x,const double s,double &f,double &fp,double &fpp){
   const double ax=fabs(x), z=SQRT_ONE_OVER_THREE*ax/s, y = z*z, s2=s*s, Phi=normalDistObj.norm_cdf(-z), phi=normalDistObj.norm_pdf(z);
   fpp = PI_OVER_SIX * y/(s2*s) * Phi * ( 8*SQRT_THREE*s*ax + (3*s2*(s2-8)-8*x*x)*Phi/phi ) * exp(2*y+0.25*s2);
   if (is_below_horizon(s)) {
      fp = 1;
      f = 0;
   } else {
      const double Phi2=Phi*Phi;
      fp = TWO_PI*y*Phi2*exp(y+0.125*s*s);
      if (is_below_horizon(x))
         f = 0;
      else
         f = TWO_PI_OVER_SQRT_TWENTY_SEVEN*ax*(Phi2*Phi);
   }
}
#endif

double inverse_f_lower_map(const double x,const double f){
   return is_below_horizon(f) ? 0 : fabs(x/(SQRT_THREE * normalDistObj.inverse_norm_cdf( std::pow( f/(TWO_PI_OVER_SQRT_TWENTY_SEVEN*fabs(x)) , 1./3.) )));
}

#ifdef COMPUTE_UPPER_MAP_DERIVATIVES_INDIVIDUALLY
double f_upper_map(const double s){
   return norm_cdf(-0.5*s);
}
double d_f_upper_map_d_beta(const double x,const double s){
   return is_below_horizon(x) ? -0.5 : -0.5*exp(0.5*square(x/s));
}
double d2_f_upper_map_d_beta2(const double x,const double s){
   if (is_below_horizon(x))
      return 0;
   const double w = square(x/s);
   return SQRT_PI_OVER_TWO*exp(w+0.125*s*s)*w/s;
}
void compute_f_upper_map_and_first_two_derivatives(const double x,const double s,double &f,double &fp,double &fpp){
   f   = f_upper_map(s);
   fp  = d_f_upper_map_d_beta(x,s);
   fpp = d2_f_upper_map_d_beta2(x,s);
}
#else
void compute_f_upper_map_and_first_two_derivatives(const double x,const double s,double &f,double &fp,double &fpp){
   f = normalDistObj.norm_cdf(-0.5*s);
   if (is_below_horizon(x)) {
      fp = -0.5;
      fpp = 0;
   } else {
      const double w = square(x/s);
      fp = -0.5*exp(0.5*w);
      fpp = SQRT_PI_OVER_TWO*exp(w+0.125*s*s)*w/s;
   }
}
#endif

double inverse_f_upper_map(double f){
   return -2.*normalDistObj.inverse_norm_cdf(f);
}

double unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
    (double beta, double x, double q /* q=±1 */, int N)
{
    lbr lbrObj;
    rationalcubic rCubicObj;
   // Subtract intrinsic.
   if (q*x>0) {
      beta = fabs(std::max(beta-normalised_intrinsic(x, q),0.));
      q = -q;
   }
   // Map puts to calls
   if (q<0){
      x = -x;
      q = -q;
   }
   if (beta<=0) // For negative or zero prices we return 0.
      return implied_volatility_output(0,0);

   if (beta<DENORMALISATION_CUTOFF) // For positive but denormalised (a.k.a. 'subnormal') prices, we return 0 since it would be impossible to converge to full machine accuracy anyway.
      return implied_volatility_output(0,0);

   const double b_max = exp(0.5*x);

   if (beta>=b_max)
      return implied_volatility_output(0,VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM);

   int iterations=0, direction_reversal_count = 0;
   double f=-DBL_MAX, s=-DBL_MAX, ds=s, ds_previous=0, s_left=DBL_MIN, s_right=DBL_MAX;
   // The temptation is great to use the optimised form b_c = exp(x/2)/2-exp(-x/2)·Phi(sqrt(-2·x)) but that would require implementing all of the above types of round-off and over/underflow handling for this expression, too.
   const double s_c=sqrt(fabs(2*x)), b_c = lbrObj.normalised_black_call(x,s_c), v_c = lbrObj.normalised_vega(x, s_c);
   // Four branches.
   if ( beta<b_c ) {
      const double s_l = s_c - b_c/v_c, b_l = lbrObj.normalised_black_call(x,s_l);
      if (beta<b_l){
         double f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2;
         compute_f_lower_map_and_first_two_derivatives(x,s_l,f_lower_map_l,d_f_lower_map_l_d_beta,d2_f_lower_map_l_d_beta2);
         const double r_ll=rCubicObj.convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0.,b_l,0.,f_lower_map_l,1.,d_f_lower_map_l_d_beta,d2_f_lower_map_l_d_beta2,true);
         f = rCubicObj.rational_cubic_interpolation(beta,0.,b_l,0.,f_lower_map_l,1.,d_f_lower_map_l_d_beta,r_ll);
         if (!(f>0)) { // This can happen due to roundoff truncation for extreme values such as |x|>500.
            // We switch to quadratic interpolation using f(0)≡0, f(b_l), and f'(0)≡1 to specify the quadratic.
            const double t = beta/b_l;
            f = (f_lower_map_l*t + b_l*(1-t)) * t;
         }
         s = inverse_f_lower_map(x,f);
         s_right = s_l;

         for (; iterations < N && fabs(ds) > (DBL_EPSILON*s); ++iterations){
            if (ds*ds_previous<0)
               ++direction_reversal_count;
            if ( iterations>0 && ( 3==direction_reversal_count || !(s>s_left && s<s_right) ) ) {
               // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
               // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
               s = 0.5*(s_left+s_right);
               if (s_right-s_left <=(DBL_EPSILON*s)) break;
               direction_reversal_count = 0;
               ds = 0;
            }
            ds_previous=ds;
            const double b = lbrObj.normalised_black_call(x,s), bp = lbrObj.normalised_vega(x, s);
            if ( b>beta && s<s_right ) s_right=s; else if ( b<beta && s>s_left ) s_left=s; // Tighten the bracket if applicable.
            if (b<=0||bp<=0) // Numerical underflow. Switch to binary nesting for this iteration.
               ds = 0.5*(s_left+s_right)-s;
            else {
               const double ln_b=log(b), ln_beta=log(beta), bpob=bp/b, h=x/s, b_halley = h*h/s-s/4, newton = (ln_beta-ln_b)*ln_b/ln_beta/bpob, halley = b_halley-bpob*(1+2/ln_b);
               const double b_hh3 = b_halley*b_halley-3*square(h/s)-0.25, hh3 = b_hh3+2*square(bpob)*(1+3/ln_b*(1+1/ln_b))-3*b_halley*bpob*(1+2/ln_b);
               ds = newton * householder_factor(newton,halley,hh3);
            }
            s += ds = std::max(-0.5*s , ds );
         }
         return implied_volatility_output(iterations,s);
      } else {
         const double v_l = lbrObj.normalised_vega(x, s_l), r_lm = rCubicObj.convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b_l,b_c,s_l,s_c,1/v_l,1/v_c,0.0,false);
         s = rCubicObj.rational_cubic_interpolation(beta,b_l,b_c,s_l,s_c,1/v_l,1/v_c,r_lm);
         s_left = s_l;
         s_right = s_c;
      }
   } else {
      const double s_h = v_c>DBL_MIN ? s_c+(b_max-b_c)/v_c : s_c, b_h = lbrObj.normalised_black_call(x,s_h);
      if(beta<=b_h){
         const double v_h = lbrObj.normalised_vega(x, s_h), r_hm = rCubicObj.convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_c,b_h,s_c,s_h,1/v_c,1/v_h,0.0,false);
         s = rCubicObj.rational_cubic_interpolation(beta,b_c,b_h,s_c,s_h,1/v_c,1/v_h,r_hm);
         s_left = s_c;
         s_right = s_h;
      } else {
         double f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2;
         compute_f_upper_map_and_first_two_derivatives(x,s_h,f_upper_map_h,d_f_upper_map_h_d_beta,d2_f_upper_map_h_d_beta2);
         if ( d2_f_upper_map_h_d_beta2>-SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2<SQRT_DBL_MAX ){
            const double r_hh = rCubicObj.convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_h,b_max,f_upper_map_h,0.,d_f_upper_map_h_d_beta,-0.5,d2_f_upper_map_h_d_beta2,true);
            f = rCubicObj.rational_cubic_interpolation(beta,b_h,b_max,f_upper_map_h,0.,d_f_upper_map_h_d_beta,-0.5,r_hh);
         }
         if (f<=0) {
            const double h=b_max-b_h, t=(beta-b_h)/h;
            f = (f_upper_map_h*(1-t) + 0.5*h*t) * (1-t); // We switch to quadratic interpolation using f(b_h), f(b_max)≡0, and f'(b_max)≡-1/2 to specify the quadratic.
         }
         s = inverse_f_upper_map(f);
         s_left = s_h;
         if (beta > (0.5*b_max)) { // Else we better drop through and let the objective function be g(s) = b(x,s)-beta.

            for (; iterations<N && fabs(ds) > (DBL_EPSILON*s); ++iterations){
               if (ds*ds_previous<0)
                  ++direction_reversal_count;
               if ( iterations>0 && ( 3==direction_reversal_count || !(s>s_left && s<s_right) ) ) {
                  // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
                  // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
                  s = 0.5*(s_left+s_right);
                  if (s_right-s_left <= (DBL_EPSILON*s)) break;
                  direction_reversal_count = 0;
                  ds = 0;
               }
               ds_previous=ds;
               const double b = lbrObj.normalised_black_call(x,s), bp = lbrObj.normalised_vega(x, s);
               if ( b>beta && s<s_right ) s_right=s; else if ( b<beta && s>s_left ) s_left=s; // Tighten the bracket if applicable.
               if (b>=b_max||bp<=DBL_MIN) // Numerical underflow. Switch to binary nesting for this iteration.
                  ds = 0.5*(s_left+s_right)-s;
               else {
                  const double b_max_minus_b = b_max-b, g = log((b_max-beta)/b_max_minus_b), gp = bp/b_max_minus_b;
                  const double b_halley = square(x/s)/s-s/4, b_hh3 = b_halley*b_halley-3*square(x/(s*s))-0.25;
                  const double newton = -g/gp, halley = b_halley+gp, hh3 = b_hh3+gp*(2*gp+3*b_halley);
                  ds = newton * householder_factor(newton,halley,hh3);
               }
               s += ds = std::max(-0.5*s , ds );
            }
            return implied_volatility_output(iterations,s);
         }
      }
   }

   for (; iterations<N && fabs(ds) > (DBL_EPSILON*s); ++iterations){
      if (ds*ds_previous<0)
         ++direction_reversal_count;
      if ( iterations>0 && ( 3==direction_reversal_count || !(s>s_left && s<s_right) ) ) {
         // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
         // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
         s = 0.5*(s_left+s_right);
         if (s_right-s_left <= DBL_EPSILON*s) break;
         direction_reversal_count = 0;
         ds = 0;
      }
      ds_previous=ds;
      const double b = lbrObj.normalised_black_call(x,s), bp = lbrObj.normalised_vega(x, s);
      if ( b>beta && s<s_right ) s_right=s; else if ( b<beta && s>s_left ) s_left=s; // Tighten the bracket if applicable.
      const double newton = (beta-b)/bp, halley = square(x/s)/s-s/4, hh3 = halley*halley-3*square(x/(s*s))-0.25;
      s += ds = std::max(-0.5*s , newton * householder_factor(newton,halley,hh3) );
   }
   return implied_volatility_output(iterations,s);
}

double lbr::implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
    (double price, double F, double K, double T, double q, int N)
{
    const double intrinsic = fabs(std::max((q<0?K-F:F-K),0.0));
    if (price<intrinsic)
       return implied_volatility_output(0,VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC);

    const double max_price = (q<0?K:F);
    if (price>=max_price)
       return implied_volatility_output(0,VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM);

    const double x = log(F/K);
    // Map in-the-money to out-of-the-money
    if (q*x>0) {
       price = fabs(std::max(price-intrinsic,0.0));
       q = -q;
    }

    return unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
            (price/(sqrt(F)*sqrt(K)), x, q, N)/sqrt(T);
}

double lbr::implied_volatility_from_a_transformed_rational_guess(double price, double F, double K, double T, double q)
{
return implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
        (price,F,K,T,q,implied_volatility_maximum_iterations.data);
}

double lbr::normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
    (double beta, double x, double q, int N)
{
    // Map in-the-money to out-of-the-money
    if (q*x>0) {
       beta -= normalised_intrinsic(x, q);
       q = -q;
    }
    if (beta<0)
       return implied_volatility_output(0,VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC);

    return unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
            (beta, x, q, N);
}

double lbr::normalised_implied_volatility_from_a_transformed_rational_guess(double beta, double x, double q)
{
return normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
        (beta,x,q,implied_volatility_maximum_iterations.data);
}
