#ifndef NORMALDIST_H
#define NORMALDIST_H

#include <cmath>

#define ONE_OVER_SQRT_TWO     0.7071067811865475244008443621048490392848359376887
#define ONE_OVER_SQRT_TWO_PI  0.3989422804014326779399460599343818684758586311649
#define SQRT_TWO_PI           2.506628274631000502415765284811045253006986740610

class normaldist
{
public:
    normaldist();

    double calerf(double x, const int jint);
    double erf_cody(double x);
    double erfc_cody(double x);
    double erfcx_cody(double x);
    double norm_cdf(double z);
    inline double norm_pdf(double x){ return ONE_OVER_SQRT_TWO_PI*exp(-.5*x*x); }
    double inverse_norm_cdf(double u);
};

#endif // NORMALDIST_H
