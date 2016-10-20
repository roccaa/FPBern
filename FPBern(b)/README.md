# FPBern(b)
## Description
`FPBern(b)` is a tool for computation of roundoff error bounds using Bernstein expansions.

We consider the roundoff error r of a program implementing a polynomial function f(x):

			r(x,e) = l(x, e) + h(x, e)

with `l(x,e)` the part of `r(x,e)` linear w.r.t. `e` and `h(x,e)` the part  of `r(x,e)` non-linear w.r.t. `e`.

`FPBern(b)` gives an upper bound of `|l(x,e)|` with x lying in a box, and e being bounded by a given machine epsilon. 

For instance, one can consider f(x) = x1^3 + (3/4) x1 x2^2, with x1 in [-1, 1] and x2 in [-2, 2].

###Program representation
To use FPBern(b), you need to perform the following:

- Declare all variables `(x,e)` symbols with `syms` from `Matlab Symbolic Toolbox`
- Define `l`, the l(x,e) polynomial
- Define `vars`, the variable `x` as an array of the above declared symbols
- Define `params`, the errors terms `e` (in the same way as above)
- Define `X`, the input box, e.g. by `X = [-5,5,-5,5];`, the 1st (resp. 2nd) variable lies in `[-5,5]`
- Call [bmax,bmin] = bernstein_method(l,params,vars,X,'polynomial') to compute the lower and the upper bounds of l(x,e)

## Installation instructions
### Prerequisites
FPBern(b) is implemented in Matlab2015a, thus it is not guaranteed it will work on a previous (or later) version.

Moreover, FPBern(b) relies on a module of Matlab:
- `Matlab symbolic toolbox`

### Download
FPBern(b) is maintained as a GitHub repository at the address https://github.com/roccaa/FPBern.git

It can be obtained by typing the following command:

$ git clone https://github.com/roccaa/FPBern.git

### Benchmarks
To run a set of benchmarks (and any further tests), launch the following commands (on Matlab):

    `$ addpath(genpath(benchmarks/)); addpath(genpath(include/)); addpath(genpath(src/));`

Eventually you can run the following commands:

 	$ do_all;
	$ result
