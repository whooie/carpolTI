// file: readfile.hpp
// author: Will Huie

#ifndef READFILE_HPP
#define READFILE_HPP

#include <map>
#include <tuple>
#include <set>
#include <string>

void trim_data(int head, int foot);

std::set<double> unique_column_vals(std::string filename, double minval, double maxval, int colnum);

std::map<double, std::map<double, double>> make_map(std::string filename, int xcolnum, int ycolnum, int depcolnum);

std::map<double, std::set<double>> make_grid_map(std::string filename, int xcolnum, double xminval, double xmaxval, int ycolnum, double yminval, double ymaxval);

double timeof(std::string filename, int hidx, int midx, int sidx, int msidx);

std::set<double> read_original_tvals(int hidx, int midx, int sidx, int msidx);

#endif
