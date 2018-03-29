// file: interpolate.cpp
// author: Will Huie

#include <stdio.h>
#include <iostream>
#include <map>
#include <tuple>
#include <set>
#include <cmath>
#include <limits>

#include "interpolate.hpp"

using namespace std;

tuple<double, double, double, double> find_neighbors(set<double> xvals, set<double> yvals, double r, double phi){
    /*
     * Finds the four nearest cartesian grid points to the
     * specified polar grid point. Assumes that (r,phi) is
     * in bounds.
     *
     * Args:
     * xvals, yvals - double arrays of the x- and y- values of
     *     the cartesian grid points
     * r, phi - double polar coordinates of the target point
     *     on the polar grid
     *
     * Returns:
     * box - tuple containing the x-, and y-values of the
     *     cartesian grid points which lie closest to the
     *     target point (xmin, xmax, ymin, ymax)
     */

    // initialize
    double xmin, xmax, ymin, ymax;
    xmax = INF;
    xmin = -1.0*INF;
    ymax = INF;
    ymin = -1.0*INF;

    // compute the target values
    double polarx, polary;
    polarx = r*cos(phi);
    polary = r*sin(phi);

    // find bounding values
    double d;
    for(set<double>::iterator it = xvals.begin(); it != xvals.end(); it++){
        d = polarx - *(it);
        if(d < 0){
            break;
        }else if(d < (polarx - xmin)){
            xmin = *it;
        }
    }
    for(set<double>::reverse_iterator it = xvals.rbegin(); it != xvals.rend(); it++){
        d = *(it) - polarx;
        if(d < 0){
            break;
        }else if(d < (xmax - polarx)){
            xmax = *it;
        }
    }
    for(set<double>::iterator it = yvals.begin(); it != yvals.end(); it++){
        d = polary - *(it);
        if(d < 0){
            break;
        }else if(d < (polary - ymin)){
            ymin = *it;
        }
    }
    for(set<double>::reverse_iterator it = yvals.rbegin(); it != yvals.rend(); it++){
        d = *(it) - polary;
        if(d < 0){
            break;
        }else if(d < (ymax - polary)){
            ymax = *it;
        }
    }

    // pack values into a tuple and return
    tuple<double, double, double, double> box;
    box = make_tuple(xmin, xmax, ymin, ymax);
    return box;
}

double interpolate_space(tuple<double, double, double, double> box, map<double, map<double, double>> datamap, double r, double phi){
    /*
     * For four given bounding points found in box, interpolate
     * the datamap value at (r,phi).
     *
     * Args:
     * box - gives the coordinates of the four bounding points
     *     in the form (xmin, xmax, ymin, ymax)
     * datamap - gives the value to be interpolated
     * r, phi - the polar coordinates of the target point
     *
     * Returns:
     * result - the interpolated value
     */

    // calculate target point
    double polarx, polary;
    polarx = r*cos(phi);
    polary = r*sin(phi);

    // unpack the bounding box
    double xmin, xmax, ymin, ymax;
    tie(xmin, xmax, ymin, ymax) = box;

    // determine the exact type of linear interpolation to use
    // based on the number of grid lines on which (r,phi) lies
    int type = 0;
    if(xmin == xmax){
        type++;
    }
    if(ymin == ymax){
        type+=2;
    }

    // type 0: (r,phi) lies on no grid lines
    if(type == 0){
        double area, areaxy, areaXy, areaxY, areaXY;
        area = (xmax - xmin)*(ymax - ymin);
        areaxy = (polarx - xmin)*(polary - ymin)/area;
        areaXy = (xmax - polarx)*(polary - ymin)/area;
        areaxY = (polarx - xmin)*(ymax - polary)/area;
        areaXY = (xmax - polarx)*(ymax - polary)/area;

        return (areaxy*datamap[xmax][ymax])+
               (areaXy*datamap[xmin][ymax])+
               (areaxY*datamap[xmax][ymin])+
               (areaXY*datamap[xmin][ymin]);
    // type 1: (r,phi) lies on a vertical grid line
    }else if(type == 1){
        double dist, disty, distY;
        dist = ymax - ymin;
        disty = (polary - ymin)/dist;
        distY = (ymax - polary)/dist;

        return (disty*datamap[xmin][ymax])+
               (distY*datamap[xmin][ymin]);
    // type 2: (r,phi) lies on a horizontal grid line
    }else if(type == 2){
        double dist, distx, distX;
        dist = xmax - xmin;
        distx = (polarx - xmin)/dist;
        distX = (xmax - polarx)/dist;

        return (distx*datamap[xmax][ymin])+
               (distX*datamap[xmin][ymin]);
    // type 3: (r,phi) lies on the intersection of grid lines
    }else if(type == 3){
        return datamap[xmin][ymin];
    }
    return 0;
}

double interpolate_time(double target_time, double time1, double val1, double time2, double val2){
    /*
     * Interpolates between temporal values.
     *
     * Args:
     * target_time - the time gridpoint between time1 and time2
     * time1, val1 - the time and value of the first gridpoint
     * time2, val2 - the time and value of the second gridpoint
     *
     * Returns:
     * result - the interpolated value
     */

    double result;
    if(time1 == time2){
        result = val1;
    }else{
        result = (target_time - time1)*val2 + (time2 - target_time)*val1;
    }
    return result;
}
