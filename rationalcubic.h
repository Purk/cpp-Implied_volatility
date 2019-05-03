#ifndef RATIONALCUBIC_H
#define RATIONALCUBIC_H


class rationalcubic
{
public:
    rationalcubic();

    double rational_cubic_interpolation(double x, double x_l, double x_r, double y_l,
                                        double y_r, double d_l, double d_r, double r);

    double rational_cubic_control_parameter_to_fit_second_derivative_at_left_side
        (double x_l, double x_r, double y_l, double y_r, double d_l,
         double d_r, double second_derivative_l);

    double rational_cubic_control_parameter_to_fit_second_derivative_at_right_side
        (double x_l, double x_r, double y_l, double y_r,
         double d_l, double d_r, double second_derivative_r);

    double minimum_rational_cubic_control_parameter
        (double d_l, double d_r, double s,
         bool preferShapePreservationOverSmoothness);

    double convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side
        (double x_l, double x_r, double y_l,
         double y_r, double d_l, double d_r,
         double second_derivative_l, bool preferShapePreservationOverSmoothness);

    double convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side
        (double x_l, double x_r, double y_l, double y_r,
         double d_l, double d_r, double second_derivative_r,
         bool preferShapePreservationOverSmoothness);

};

#endif // RATIONALCUBIC_H
