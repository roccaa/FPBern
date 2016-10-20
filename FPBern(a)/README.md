# FPBern(a)
## Description
`FPBern(a)` is a tool for computation of roundoff error bounds using Bernstein expansions.

We consider the roundoff error r of a program implementing a polynomial function f(x):

			r(x,e) = l(x, e) + h(x, e)

with l(x, e) the part of r(x, e) linear w.r.t. e and h(x, e) the part  of r(x, e) non-linear w.r.t. e.

`FPBern(a)` gives an upper bound of `|l(x,e)|` with x lying in a box, and e being bounded by a given machine epsilon. 

For instance, one can consider f(x) = x1^3 + (3/4) x1 x2^2, with x1 in [-1, 1] and x2 in [-2, 2].

###Programs representation
FPBern(a) handles input files with .ini extension (mandatory) with the following structure:

 + OPTIONS
- name = `name of the program`
- precision = `allowing to define a machine epsilon equal to 2^(-precision)`
- nbvars = `dimension of x`
- nberrors = `dimension of e`
 + Programs
- function = ` polynomials function. The operations +,*,- and ^ (and / in the coefficients) are supported`
- input_bl = `lower bounds of the input values`
- input_bu = `upper bounds of the input values`

## Installation instructions
### Prerequisites
FPBern(a) is implemented in C++. Thus, a C++ compiler is required.
FPBern(a) was tested on Ubuntu 14.04 LTS.

Moreover, FPBern(a) relies on three external libraries:
- GiNaC (GiNaC is Not a CAS), for the symbolic manipulation of polynomials
- GLPK (GNU Linear Programming Kit), for solving linear programming problems
- boost (boost/property_tree/ptree.hpp and boost/property_tree/ini_parser.hpp and boost/lexical_cast.hpp)

###Download
FPBern(a) is maintained as a GitHub repository at the address https://github.com/roccaa/FPBern.git

It can be obtained either by typing the following command:

$ git clone https://github.com/roccaa/FPBern.git

###Installation
To install from the source type:

	$ make

This creates a binary called FPBern in /bin.

To run FPBern, move to /bin and launch the FPBern binary with the command:

	$ ./FPBern file_name1 file_name2 ...    (without the .ini)

To run aset of benchmarks, launch the command

 	$ ./FPBern rigidbody1 rigidbody2 kepler0 kepler1 kepler2 sineTaylor sineOrder3 sqroot himmilbeau schwefel magnetism caprasse
