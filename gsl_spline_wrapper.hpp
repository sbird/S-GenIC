/** A short class to wrap the GSL spline routines in a more C++ friendly format
Avoids all the tedious freeing and allocation.
Uses an opaque shared pointer type to implement implicitly shared data and so requires C++11
All functions defer to the private wrapper.
*/

#include <memory>

class gsl_spline_wrapper_private;

class gsl_spline_wrapper
{
    public:
        //Initialize everything: note that we keep a copy of the arrays so you can free them after calling this
        gsl_spline_wrapper(const double *xval, const double *yval, const int nval);

        /**Evaluate the spline at a particular x value:
         * NOT checked that this is within the original range! */
        double eval(const double x);
        //Reset the accelerator for the spline
        void reset(void);

        //Get the bounds
        double xmax() const;

        double xmin() const;
    private:
        //d-pointer: note this is implicitly shared!
        std::shared_ptr<gsl_spline_wrapper_private> d;
};
