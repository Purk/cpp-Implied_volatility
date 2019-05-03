#ifndef LBR_H
#define LBR_H

#define ENABLE_SWITCHING_THE_OUTPUT_TO_ITERATION_COUNT
#define ENABLE_CHANGING_THE_HOUSEHOLDER_METHOD_ORDER

#include "normaldist.h"

class lbr
{
public:
    lbr();

    double set_implied_volatility_maximum_iterations(double t);
    double set_implied_volatility_output_type(double t);
    double set_implied_volatility_householder_method_order(double t);
    double normalised_black_call(double x, double s);
    double normalised_vega(double x, double s);
    double normalised_black(double x, double s, double q /* q=±1 */);
    double black(double F, double K, double sigma, double T, double q /* q=±1 */);


    double implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
        (double price, double F, double K, double T, double q /* q=±1 */, int N);

    double implied_volatility_from_a_transformed_rational_guess(double price, double F,
                                                                double K, double T,
                                                                double q /* q=±1 */);

    double normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations
        (double beta, double x, double q /* q=±1 */, int N);

    double normalised_implied_volatility_from_a_transformed_rational_guess
        (double beta, double x, double q /* q=±1 */);
};

#endif // LBR_H
