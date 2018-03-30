// file: min.cpp
// author: Will Huie

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <map>
#include <cmath>
#include <dirent.h>

#include "readfile.hpp"
#include "interpolate.hpp"
#include "transform.hpp"

using namespace std;

int main(){
    /*
     * Interpolates spacial grid points given in the rlist and 
     * philist files from the data in the files in the directory
     * data_raw (in Raluca Ilie's format), and writes them to 
     * an output file in the format useful to Scot Elkington. 
     * Assumes the data has been trimmed. This version does not
     * perform interpolation in time, and will write r, phi,
     * and time values per-line to file in the output directory.
     */

    // initialize our constants up here for easy adjustment;
    // remember to edit these lines to your specifications!

    // first and last line numbers of the data (inclusive)
    const int head = 24;
    const int foot = 92687;

    // r-value range (inclusive)
    const double rmin = 0.0;
    const double rmax = 8.0;

    // phi-value range (inclusive)
    const double phimin = 0;
    const double phimax = 2*M_PI;

    // x-value range (inclusive)
    const double xmin = -8.0;
    const double xmax = 8.0;

    // y-value range (inclusive)
    const double ymin = -8.0;
    const double ymax = 8.0;

    // indices for hour, minute, second, and millisecond values
    // in the titles of each file (starting at 0)
    const int hidx = 20;
    const int midx = 22;
    const int sidx = 24;
    const int msidx = 27;

    // column numbers of the independent variables (starts at 0
    // on the left)
    const int xcolnum = 0;
    const int ycolnum = 1;

    // column numbers of the dependent variables (starts at 0
    // on the left)
    const int Excolnum = 11; // E_x
    const int Eycolnum = 12; // E_y
    const int Bzcolnum = 9; // B_z

    // paths to necessary input files
    const string rlist = "gridlists/rlist.txt";
    const string philist = "gridlists/philist.txt";

    // paths to necessary output files
    const string outfile = "output/output.dat";
    const string rkey = "output/rvalues.txt";
    const string phikey = "output/phivalues.txt";
    const string tkey = "output/tvalues.txt";

    /***************/

    // trim the data
    cout << "Using files in ./data_raw/" << endl;
    /*
     * trimming only needs to be done once if you're working
     * with the same data, so comment out these two lines for
     * subsequent run over the same set of files
     */
    cout << "Trimming data with head=" << head << " and foot=" << foot <<endl;
    trim_data(head, foot);

    // open the output file
    ofstream output;
    output.open(outfile);
    output << fixed << setprecision(7) << scientific;

    cout << fixed << setprecision(7) << scientific;
    cout << "Getting polar gridpoints in the ranges" << endl;
    cout << "    r: [" << rmin << "," << rmax << "]" << endl;
    cout << "  phi: [" << phimin << "," << phimax << "]" << endl;
    // get out polar gridpoints and target times
    set<double> rvals = unique_column_vals(rlist, rmin, rmax, 0);
    vector<double> rvalsvec;
    for(set<double>::iterator it = rvals.begin(); it != rvals.end(); it++){
        rvalsvec.push_back(*it);
    }

    set<double> phivals = unique_column_vals(philist, phimin, phimax, 0);
    vector<double> phivalsvec;
    for(set<double>::iterator it = phivals.begin(); it != phivals.end(); it++){
        phivalsvec.push_back(*it);
    }

    // get a list of the times we'll be operating over
    set<double> tvals = read_original_tvals(hidx, midx, sidx, msidx);
    
    cout << "Accepted polar gridpoints can be found in ./output/" << endl;
    // write the r, phi, and t values to their respective files
    // and write the accepted r, phi, and t values to the output
    // files with associated headers
    ofstream rfile, phifile, tfile;
    rfile.open(rkey);
    phifile.open(phikey);
    tfile.open(tkey);
    bool line;

    rfile << "INDEX   R-VALUES (Re)" << endl;
    rfile << fixed << setprecision(7) << scientific;
    output << "RADIAL DISTANCES (Re) - " << rvals.size() << endl;
    set<double>::iterator r_it = rvals.begin();
    for(unsigned int i = 0; i < rvals.size(); i++,r_it++){
        rfile << setw(5) << i+1 << setw(17) << *r_it << endl;
        output << setw(17) << *r_it;
        line = true;
        if((i+1)%5 == 0){
            output << endl;
            line = false;
        }
    }
    if(line){
        output << endl;
    }

    phifile << "INDEX   PHI-VALUES (rad)" << endl;
    phifile << fixed << setprecision(7) << scientific;
    output << "ANGULAR POINTS (rad) - " << phivals.size() << endl;
    set<double>::iterator phi_it = phivals.begin();
    for(unsigned int i = 0; i < phivals.size(); i++,phi_it++){
        phifile << setw(5) << i+1 << setw(17) << *phi_it << endl;
        output << setw(17) << *phi_it;
        line = true;
        if((i+1)%5 == 0){
            output << endl;
            line = false;
        }
    }
    if(line){
        output << endl;
    }

    tfile << "INDEX   TIME VALUES (s)" << endl;
    tfile << fixed << setprecision(7) << scientific;
    output << "TIME VALUES (s) - " << tvals.size() << endl;
    set<double>::iterator t_it = tvals.begin();
    for(unsigned int i = 0; i < tvals.size(); i++,t_it++){
        tfile << setw(5) << i+1 << setw(17) << *t_it << endl;
        output << setw(17) << *t_it;
        line = true;
        if((i+1)%5 == 0){
            output << endl;
            line = false;
        }
    }
    if(line){
        output << endl;
    }

    rfile.close();
    phifile.close();
    tfile.close();

    // get a set of all the files in data_trimmed/
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

    cout << "Getting rectangular girdpoints in the ranges" << endl;
    cout << "    x: [" << xmin << "," << xmax << "]" << endl;
    cout << "    y: [" << ymin << "," << ymax << "]" << endl;

    // get (x,y) gridpoints, add them to a vector for indexing
    set<double> xvals = unique_column_vals("data_trimmed/"+(*filelist.begin()), xmin, xmax, xcolnum);
    set<double> yvals = unique_column_vals("data_trimmed/"+(*filelist.begin()), ymin, ymax, ycolnum);

    // write the main header of the output file
    output << "FIRST LINE: i = radial index, j = angular index, k = time index; values are listed in accompanying files." << endl;
    output << "SECOND LINE: Ephi (V/m), Er (V/m), Bz (nT), d(Bz)/dr (nT/Re), d(Bz)/dphi (nT/rad)" << endl;

    cout << "Begin calculating values" << endl;
    // for each time file, do
    int k = 0;
    for(string filename : filelist){
        cout << "  " << filename << " (" << k+1 << "/" << tvals.size() << ")" << endl;

        cout << "    Generating maps for Ex, Ey, and Bz" << endl;
        // make maps for Ex, Ey, and Bz
        map<double, map<double, double>> Ex_vals, Ey_vals, Bz_vals;
        Ex_vals = make_map("data_trimmed/"+filename, xcolnum, ycolnum, Excolnum);
        Ey_vals = make_map("data_trimmed/"+filename, xcolnum, ycolnum, Eycolnum);
        Bz_vals = make_map("data_trimmed/"+filename, xcolnum, ycolnum, Bzcolnum);

        cout << "    Interpolating and transforming" << endl;
        // for each (r,phi) gridpoint, do
        unsigned int i, j;
        //vector<double>::iterator r_it = rvalsvec.begin();
        //vector<double>::iterator phi_it = phivalsvec.begin();
        for(j = 0; j < phivalsvec.size(); j++){
            for(i = 0; i < rvalsvec.size(); i++){
                //printf("      Grid point %i/%i\r", j*rvals.size()+(i+1), phivals.size()*rvals.size());
                cout << "      Grid point " << j*rvals.size()+(i+1) << "/" << phivals.size()*rvals.size() << "\r";
                // get the (r,phi) coords, as well as hi/lo
                // coords for calculating the Bz derivatives
                double r = rvalsvec[i];
                double r_lo = i == 0 ? rvalsvec[i] : rvalsvec[i-1];
                double r_hi = i == rvalsvec.size() - 1 ? rvalsvec[i] : rvalsvec[i+1];

                double phi = phivalsvec[j];
                double phi_lo = j == 0 ? phivalsvec[j] : phivalsvec[j-1];
                double phi_hi = j == phivalsvec.size() - 1 ? phivalsvec[j] : phivalsvec[j+1];

                // find the bounding (x,y) coords
                tuple<double, double, double, double> box, box_r_lo, box_r_hi, box_phi_lo, box_phi_hi;
                box = find_neighbors(xvals, yvals, r, phi);
                box_r_lo = find_neighbors(xvals, yvals, r_lo, phi);
                box_r_hi = find_neighbors(xvals, yvals, r_hi, phi);
                box_phi_lo = find_neighbors(xvals, yvals, r, phi_lo);
                box_phi_hi = find_neighbors(xvals, yvals, r, phi_hi);

                // interpolate our target values
                double Ex, Ey, Bz, Bz_r_lo, Bz_r_hi, Bz_phi_lo, Bz_phi_hi;
                Ex = interpolate_space(box, Ex_vals, r, phi);
                Ey = interpolate_space(box, Ey_vals, r, phi);
                Bz = interpolate_space(box, Bz_vals, r, phi);
                Bz_r_lo = interpolate_space(box_r_lo, Bz_vals, r_lo, phi);
                Bz_r_hi = interpolate_space(box_r_hi, Bz_vals, r_hi, phi);
                Bz_phi_lo = interpolate_space(box_phi_lo, Bz_vals, r, phi_lo);
                Bz_phi_hi = interpolate_space(box_phi_hi, Bz_vals, r, phi_hi);

                // transform <Ex,Ey> to <Er,Ephi>
                double Er, Ephi;
                tie(Er, Ephi) = transform_vector(Ex, Ey, r, phi);
                // convert from mV/m to V/m
                Er = Er/1000;
                Ephi = Ephi/1000;

                // calculate d(Bz)/dr and d(Bz)/dphi
                double dBzdr, dBzdphi;
                dBzdr = (Bz_r_hi - Bz_r_lo) / (r_hi - r_lo);
                dBzdphi = (Bz_phi_hi - Bz_phi_lo) / (phi_hi - phi_lo);

                // write to the output file
                output << setw(8) << i+1 << setw(8) << j+1 << setw(8) << k+1 << endl;
                output << setw(17) << Ephi;
                output << setw(17) << Er;
                output << setw(17) << Bz;
                output << setw(17) << dBzdr;
                output << setw(17) << dBzdphi << endl;
            }
        }
        k++;
        cout << endl;
    }
    output.close();
    return 0;
}
