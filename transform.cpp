// file: transform.hpp
// author: Will Huie

#include <stdio.h>
#include <iostream>
#include <tuple>
#include <cmath>

#include "transform.hpp"

using namespace std;

tuple<double, double> transform_vector(double xcomp, double ycomp, double r, double phi){
    /*
     * Transforms the given x-, y-components of a vector in 
     * Cartesian coordinates to the r- and phi-components of a 
     * vector in polar coordinates.
     *
     * Args:
     * xcomp, ycomp - the x-, y-components of the vector in 
     *     Cartesian coordinates
     * r, phi - the coordinates of the tail of the vector
     *
     * Returns:
     * newvals - tuple containing the transformed 2D vector
     *     components
     */

    //initialize our values
    double x = r*cos(phi);
    double y = r*sin(phi);
    double rcomp, phicomp;
    tuple<double, double> newvals;

    //dot the cartesian vectors with the polar unit vectors
    rcomp = xcomp*x/r + ycomp*y/r;
    phicomp = -1.0*xcomp*y/r + ycomp*x/r;

    //pack it into a tuple and return
    newvals = make_tuple(rcomp, phicomp);
    return newvals;
}
