Summary of carpolTI
Will Huie

=== 1 Introduction ===
    carpolTI, short for Cartesian-to-Polar Transformer and Interpolator (tentatively), is a set of
methods written in C++ with the original intention of translating data from a uniform Cartesian grid
to a polar one. While this document will provide a quick run-down of what carpolTI was originally
written to di, it should be noted that these methods were written with flexibility in mind, and
should be easily adaptable to arbitrary cases.

=== 2 Program Overview and Initial Setup ===
    This section goes over the general steps of the process and particular arrangement of files that
carpolTI requires.

= 2.1 The General Method =
    Broadly speaking, the overall procedure of carpolTI can be given by the following pseudocode.

----------------------------------------------------------------------------------------------------
Algorithm: carpolTI
----------------------------------------------------------------------------------------------------
1 Trim headers and footers off the input data;
2 Generate set PHI of (r,phi,t) polar grid points;
3 foreach (r,phi,t) in PHI do
4   Generate maps E:(x,y,t) -> (Ex,Ey), B:(x,y,t) -> Bz of electric and magnetic field values,
      respectively;
5   Using E and B, interpolate values at (r,phi,t);
6   Transform <Ex,Ey> -> <Er,Ephi>;
7   Using B, calculate d(Bz)/dr, d(Bz)/dphi;
8   Write results to output file;
----------------------------------------------------------------------------------------------------

= 2.2 Preparation =
    Before beginning, however, carpolTI requires a specific file structure. To start, run in a Bash
terminal to generate the necessary directories:

  $ cd /path/to/carpolTI
  $ make install

Then obey the following guidelines:
  - data files should be placed in the data_raw/ directory
  - r- and phi-values should be listed per-line in separate files in the gridlists/ directory

= 2.3 Setting Constants and Compiling =
    After files are in their places, main() (in main.cpp) should be edited to set the following
constants:
  - head, foot: starting and ending line numbers (with 1 at the top) of the desired data in each file.
  - rmin, rmax: the range (inclusive) of desired r-bounds in the polar grid.
  - phimin, phimax: the desired bounds (inclusive) in the polar grid.
  - xmin, xmax, ymin, ymax: the desired bounds (inclusive) in the Cartesian grid.
  - hidx, midx, sidx, msidx: the string indices (with 0 on the far left) of the hour, minute, second,
    and millisecond numbers in the filenames' timestamps.
  - xcolnum, ycolnum, Excolnum, Eycolnum, Bzcolnum: the column indices (with 0 on the far left) of
    x-, y-, Ex-, Ey-, and Bz-values.
  - outfile: path to the output file for the run.
  - rlist, philist: paths to the files containing the r- and phi-values of the polar grid points.
  - rkey, phikey, tkey: paths to the files which list the r-, phi-, and t-values used in the run.
    Note that these values will also be printed in the main output file.
Once these are in order, compile and run with
  
  $ cd /path/to/carpolTI
  $ make
  $ ./main

=== 3 The Complete Procedure ===
    This section will go into greater detail on the general method and the individual tools it uses.
The complete signatures for each method mentioned below can be found below in Section 4.

= 3.1 Trimming the Data Files =
    This action is performed by the single method trim_data(). For each file listed in data_raw/,
trim_data() discards all lines outside of the specified range and copies the remaining part to the
trimmed_data/ directory under the same filename. As most text editors consider the first line of a
file to be line 1, trim_data() will do so as well.

= 3.2 Reading Values =
    After the data has been trimmed, each dile is to be looped over and needs to be scanned for a
few things. That is, for each of the files in data_raw/, do everything that follows in the rest of 
Section 3.
 1. Get sorted arrays of unique coordinate values in the Cartesian grid; accomplished via 
    unique_column_vals(). This returns a sorted array of doubles pointing to all the unique values
    in the n-th column (with 0 being the leftmost column) of a file in a given range, along with its
    size for convenience. This should be called twice for both x- and y-values.
 2. Get maps assigning dependent variable values to (x,y) grid points; accomplished via make_map().
    This returns a once-nested map object assigning values in two columns to those in a third.
    make_map() should be called for as many dependent variables as are desired - in this case, once
    for each of Ex, Ey, and Bz.
Once this is done, unique_column_vals() should be called a few more times on the r- and phi-value
lists in gridlists/.

= 3.3 Doing the Math =
    To aid in the process of eventually writing calculated values to a file along with associated
indices, the arrays of r- and phi-values are to be looped over for each file in data_raw/. The 
values of the dependent variables read above are to be inteprolated linearly at each (r,phi) grid
point. In the case of Ex and Ey, the interpolated values are subsequently transformed to their
corresponding values in polar coordinates. Since the partial derivatives of Bz in r and phi are
additionally required, we also need four other (r,phi) points at which Bz must be calculated. This
process is as follows for each (r_i, phi_i):
 1. Use find_neighbors() to get the values of the four (x,y) grid points closest to (r_i, phi_i).
 2. Use interpolate_space() to calculate Ex(r_i, phi_i), Ey(r_i, phi_i), and Bz(r_i, phi_i).
 3. From the arrays of r- and phi-values, get r_i-1, r_i+1, phi_i-1, and phi_i+1.
 4. Use find_neighbors() and interpolate_space() to calculate Bz(r_i-1, phi_i), Bz(r_i+1, phi_i),
    Bz(r_i, phi_i-1), and Bz(r_i, phi_i+1).
 5. Use these values to calculate d(Bz)/dr and d(Bz)/dphi using a first-order central finite
    difference:

    d(Bz)/dr = (Bz(r_i+1, phi_i) - Bz(r_i-1, phi_i)) / (r_i+1 - r_i-1)
    d(Bz)/dphi = (Bz(r_i, phi_i+1) - Bz(r_i, phi_i-1)) / (phi_i+1 - phi_i-1)

 6. Write the r, phi, and t indices as well as Ephi, Er, Bz, d(Bz)/dr, and d(Bz)/dphi to the output
    file.

=== 4 Source Files ===
    This section will list the complete signatures of all written methods and their locations in
the source header files, should there ever be a need to adapt the main method. For extended
documentation, see their respective implementation files.

= 4.1 readfile.hpp =
    readfile methods are used to parse input files. Here, line numbers start at 1 from the top while
column numbers start at 0 from the left.

void trim_data(int head, int foot);

std::set<double> unique_column_vals(std::string filename, double minval, double maxval, int colnum);

std::map<double, std::map<double, double>> make_map(std::string filename, int xcolnum, int ycolnum, int depcolnum);

std::map<double, std::set<double>> make_grid_map(std::string filename, int xcolnum, double xminval, double xmaxval, int ycolnum, double yminval, double ymaxval);

double timeof(std::string filename, int hidx, int midx, int sidx, int msidx);

std::set<double> read_original_tvals(int hidx, int midx, int sidx, int msidx);

= 4.2 interpolate.hpp =
    interpolate methods deal with everything that has to do with the interpolation process. While a
method for interpolation in time is included, it is not currently used in the main method.

std::tuple<double, double, double, double> find_neighbors(std::set<double> xvals, std::set<double> yvals, double r, double phi);

double interpolate_space(std::tuple<double, double, double, double> box, std::map<double, std::map<double, double>> datamap, double r, double phi);

double interpolate_time(double target_time, double time1, double val1, double time2, double val2);

= 4.3 transform.hpp =
    The only transform method is the one used to transform from Cartesian coordinates to polar 
coordinates.

std::tuple<double, double> transform_vector(double xcomp, double ycomp, double r, double phi);

