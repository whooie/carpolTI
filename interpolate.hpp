// file: interpolate.hpp
// author: Will Huie

#ifndef INTERPOLATE_HPP
#define INTERPOLATE_HPP

#include <map>
#include <tuple>
#include <limits>
#include <set>

#define INF numeric_limits<double>::infinity()

std::tuple<double, double, double, double> find_neighbors(std::set<double> xvals, std::set<double> yvals, double r, double phi);

double interpolate_space(std::tuple<double, double, double, double> box, std::map<double, std::map<double, double>> datamap, double r, double phi);

double interpolate_time(double target_time, double time1, double val1, double time2, double val2);

#endif
