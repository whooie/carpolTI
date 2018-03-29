// file: readfile.hpp
// author: Will Huie

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <string>
#include <dirent.h>

#include "readfile.hpp"

using namespace std;

void trim_data(int head, int foot){
    /*
     * For each file specified per-line in the filelist, trim 
     * the file down to the lines between the head-th and 
     * tail-th lines (inclusive), and copy the remaining data 
     * to the trimmed_data directory. Original files should be
     * placed in the data directory to be copied to the 
     * trimmed_data directory for convenience.
     *
     * Args:
     * head, foot - line numbers of the first and last lines of
     *     the body of data
     */

    // get a set of filenames
    set<string> filelist;
    DIR* data_raw = opendir("./data_raw");
    struct dirent* file;
    while((file = readdir(data_raw)) != NULL){
        string filename = string(file -> d_name);
        if(filename.compare("..") != 0 && filename.compare(".") != 0){
            filelist.insert(filename);
        }
    }
    closedir(data_raw);
    
    // go through the file list
    for(string filename : filelist){
        ifstream infile;
        infile.open("data_raw/"+filename);

        // open the trimmed data file
        string newfilename = "data_trimmed/"+filename;
        ofstream outfile;
        outfile.open(newfilename);

        // go line-by-line through the data file and write the
        // appropriate lines to the trimmed data file
        int linenum = 1;
        string line;
        while(getline(infile, line)){
            if(linenum >= head && linenum <= foot){
                outfile << line << endl;
            }
            linenum++;
        }

        // clean up
        infile.close();
        outfile.close();
    }
}

set<double> unique_column_vals(string filename, double minval, double maxval, int colnum){
    /* 
     * find all unique values from a given column number in a 
     * given range, and sort them into an array.
     *
     * Args:
     * filename - name of the file as a string
     * min, max - the range (inclusive) of values you want
     *     from the data
     * colnum - the nth column of data to take from
     *     (leftmost column is 0)
     *
     * Returns:
     * unique_vals - set containing the pulled values
     */

    // use a set to filter/sort values
    set<double> unique_vals;

    // open the data file
    ifstream infile;
    infile.open(filename);

    // go line-by-line through the data file
    string line;
    while(getline(infile, line)){
        stringstream ssin(line);

        // repeatedly overwrite colval with numbers on the line
        // until you reach the target colnum, then add its
        // value to the set
        double colval;
        int i = 0;
        while(ssin.good() && i < colnum+1){
            ssin >> colval;
            i++;
        }
        if(colval <= maxval && colval >= minval){
            unique_vals.insert(colval);
        }
    }

    // clean up and return
    infile.close();
    return unique_vals;
}

map<double, map<double, double>> make_map(string filename, int xcolnum, int ycolnum, int depcolnum){
    /*
     * create a nested map assigning x,y coords to some given
     * dependent variable.
     *
     * Args:
     * filename - name of the file as a string
     * xcolnum, ycolnum - the column numbers of the x and y
     *     variables, respectively
     * dep_colnum - column number of the dependent variable
     *
     * Returns:
     * datamap - once-nested map of all double keys and values
     *
     * note: column numbering begins at 0 on the left.
     */

    // initialize the map
    map<double, map<double, double>> datamap;

    // open the data file
    ifstream infile;
    infile.open(filename);

    // go line-by-line through the data file
    string line;
    while(getline(infile, line)){
        stringstream ssin(line);

        // iterate through the colnums until we get to the
        // target columns, at which point set the approprate
        // value and add them all to the map
        double x, y, d, dummy;
        int i = 0;
        while(ssin.good() && i < depcolnum+1){
            ssin >> dummy;
            if(i == xcolnum){
                x = dummy;
            }else if(i == ycolnum){
                y = dummy;
            }else if(i == depcolnum){
                d = dummy;
            }
            i++;
        }
        datamap[x][y] = d;
    }

    // clean up and return
    infile.close();
    return datamap;
}

map<double, set<double>> make_grid_map(string filename, int xcolnum, double xminval, double xmaxval, int ycolnum, double yminval, double ymaxval){
    /*
     * Create a map to assign valid y-values to x-values.
     *
     * Args:
     * filename - name of the file as a string
     * xcolnum, ycolnum - the column numbers of the x- and y-
     *     variables, respectively
     * xmin, xmax, ymin, ymax - the lower and upper bounds for
     *     x- and y-values
     *
     * Returns:
     * gridmap - map of doubles assigned to double sets
     */

    // initialize the map
    map<double, set<double>> gridmap;
    
    // open the data file
    ifstream infile;
    infile.open(filename);

    // go line-by-line through the data file
    string line;
    while(getline(infile, line)){
        stringstream ssin(line);

        // iterate throughthe colnums until we get to the
        // target columns, at which point set the appropriate
        // value and add them all to the map
        double x, y, dummy;
        int i = 0;
        while(ssin.good() && i < (ycolnum > xcolnum ? ycolnum : xcolnum)){
            ssin >> dummy;
            if(i == xcolnum){
                x = dummy;
            }else if(i == ycolnum){
                y = dummy;
            }
        }

        // if x and y are within bounds, then add them to the
        // map
        if(xminval <= x && x <= xmaxval && yminval <= y && y <= ymaxval){
            gridmap[x].insert(y);
        }
    }

    // clean up and return
    infile.close();
    return gridmap;
}

double timeof(string filename, int hidx, int midx, int sidx, int msidx){
    /*
     * Reads the timestamp from a file's name and returns the 
     * time as a double.
     *
     * Args:
     * filename - the name of the file
     * *idx - the indices of the hour, minute, second, and
     *     millisecond digits in filename, respectively; hours,
     *     minutes, and seconds are assumed to be two digits 
     *     long, while milliseconds is assumed to be three
     *     digits long
     *
     *
     * Returns:
     * dtime - the time from filename in seconds as a double
     */

    string hh = filename.substr(hidx, 2);
    string mm = filename.substr(midx, 2);
    string ss = filename.substr(sidx, 2);
    string mss = filename.substr(msidx, 3);

    double dtime = 3600*stod(hh) + 60*stod(mm) + stod(ss) + 0.001*stod(mss);

    return dtime;
}

set<double> read_original_tvals(int hidx, int midx, int sidx, int msidx){
    /*
     * Creates an array of time values from the original grid.
     *
     * Args:
     * filelist - the path to the filelist
     * *idx - the indicies of the hour, minute, second, and
     *     millisecond digits in the filenames, respectively;
     *     hours, minutes, seconds, are assumed to be two
     *     digits long, while milliseconds is assumed to be
     *     three digits long
     * tmin, tmax - the bounds (inclusive) of the range of
     *     desired values
     *
     * Returns:
     * unique_vals - tuple containing the values pulled
     */

    //construct a set to filter out repeats and for easy sorting
    set<double> unique_vals;
    
    // get a set of filenames
    set<string> filelist;
    DIR* data_trimmed = opendir("./data_trimmed");
    struct dirent* file;
    while((file = readdir(data_trimmed)) != NULL){
        string filename = string(file -> d_name);
        if(filename.compare("..") != 0 && filename.compare(".") != 0){
            filelist.insert(filename);
        }
    }
    closedir(data_trimmed);

    // go through the file list
    for(string filename : filelist){
        ifstream infile;
        infile.open("data_trimmed"+filename);

        // calculate the file's time in seconds and add it to
        // unique_vals
        double timeval = timeof(filename, hidx, midx, sidx, msidx);
        unique_vals.insert(timeval);
    }

    // clean up and return
    return unique_vals;
}
