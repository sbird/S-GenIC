#include "gsl_spline_wrapper.hpp"
#include <gsl/gsl_interp.h>
#include <cassert>

//A simple class to wrap a gsl spline in a C++ class and thus avoid all the tedious freeing and allocation.
class gsl_spline_wrapper_private
{
    public:
        //Initialize everything: note that we keep a copy of the arrays so you can free them after calling this
        gsl_spline_wrapper_private(const double *xval, const double *yval, const int nval): nval(nval)
        {
            do_alloc(nval);
            for(int i=0; i< nval; i++){
                m_xval[i] = xval[i];
                m_yval[i] = yval[i];
            }
            gsl_interp_init(spline_pointer, m_xval, m_yval, nval);
        }

        //Copy constructor
        gsl_spline_wrapper_private(const gsl_spline_wrapper_private &obj): nval(obj.nval)
        {
            do_alloc(nval);
            for(int i=0; i< nval; i++){
                m_xval[i] = obj.m_xval[i];
                m_yval[i] = obj.m_yval[i];
            }
            gsl_interp_init(spline_pointer, m_xval, m_yval, nval);
        }

        //Copy assignment
        gsl_spline_wrapper_private & operator=(const gsl_spline_wrapper_private &obj)
        {
            if (this == &obj) return *this;
            nval = obj.nval;
            gsl_interp_free(spline_pointer);
            gsl_interp_accel_free(spline_accel);
            free(m_xval);
            free(m_yval);
            do_alloc(nval);
            for(int i=0; i< nval; i++){
                m_xval[i] = obj.m_xval[i];
                m_yval[i] = obj.m_yval[i];
            }
            return *this;
        }

        //Destructor
         ~gsl_spline_wrapper_private()
        {
            gsl_interp_free(spline_pointer);
            gsl_interp_accel_free(spline_accel);
            free(m_xval);
            free(m_yval);
        }

        /**Evaluate the spline at a particular x value:
         * NOT checked that this is within the original range!
        */
        double eval(const double x)
        {
            return gsl_interp_eval(spline_pointer, m_xval, m_yval, x,spline_accel);
        }
        //Reset the accelerator for the spline
        void reset(void)
        {
            gsl_interp_accel_reset(spline_accel);
        }

        //Get the bounds
        double xmax()
        {
            return m_xval[nval-1];
        }

        double xmin()
        {
            return m_xval[0];
        }

    private:
        gsl_interp * spline_pointer;
        gsl_interp_accel* spline_accel;
        double * m_xval, * m_yval;
        //Cannot be const because operator=needs to change it
        int nval;
        //Little private function to do allocations common to all constructors
        void do_alloc(int nval)
        {
            spline_pointer = gsl_interp_alloc(gsl_interp_cspline, nval);
            assert(spline_pointer);
            spline_accel = gsl_interp_accel_alloc();
            assert(spline_accel);
            m_xval = (double *) calloc(nval, sizeof(double));
            assert(m_xval);
            m_yval = (double *) calloc(nval, sizeof(double));
            assert(m_yval);
        }

};

//Here come definitions for the gsl_spline_wrapper classes

gsl_spline_wrapper::gsl_spline_wrapper(const double *xval, const double *yval, const int nval): d(new gsl_spline_wrapper_private(xval, yval, nval) )
{  };

double gsl_spline_wrapper::eval(const double x)
{
    return d->eval(x);
}

//Reset the accelerator for the spline
void gsl_spline_wrapper::reset(void)
{
    d->reset();
};

//Get the bounds
double gsl_spline_wrapper::xmax()
{
    return d->xmax();
}

double gsl_spline_wrapper::xmin()
{
    return d->xmin();
};
