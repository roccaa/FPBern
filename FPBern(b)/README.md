# FPBern(b)
## Description
`FPBern(b)` is a tool for computation of roundoff error bounds using Bernstein expansion.

If the roundoff error of a program implementing a polynomail function f(x)  is given by r(x,e) s.t.:

			r(x,e) = l(x,e)+h(x,e)

with `l(x,e)` the part of `r(x,e)` linear w.r.t `e` and `h(x,e)` the part  of `r(x,e)` non-linear in `e`.

Then `FPBern(b)` gives an upper bound of `|l(x,e)|` for `x` laying inside a box, and `e` enclose by a given epsilon. 

Thus, the semantic of the handled programs is:

- polynomial functions taken into isolation

ex: x1^3 + (3/4)*x1*x2^2
- inputs laying inside a box

ex: x1 = [-1 1] and x2 = [-2 2]

###Programs representation
FPBern(b) is use in Matlab script files:

- Declare all variables `(x,e)` symbols with `syms` from `Symbolic Toolbox`
- Define `l`, the l(x,e) polynomial
- Define `vars`, the variable `x` as an array of the previsously declared symbols
- Define `params`, the errors terms `e` in the same way
- Define `X`, the input as a box: ex `X = [-5,5,-5,5];` means 1st variable in in `[-5,5]` as well as the second variables.
- call [bmax,bmin] = bernstein_method(l,params,vars,X,'polynomial') to compute the lower and the upper bounds of l(x,e).

## Installation instructions
### Prerequisites
FPBern(b) is implemented in Matlab2015a, thus it is not guaranteed it will work on a previous (or later) version.

Moreover, FPBern(b) relies on a module of Matlab:
- `Matlab symbolic toolbox`

###Download
FPBern(b) is maintained as a GitHub repository at the address https://github.com/roccaa/FPBern.git

It can be obtained either by typing the shell command:

$ git clone https://github.com/roccaa/FPBern.git

or by downloading the ZIP archive at https://github.com/roccaa/FPBern.git

### Benchmarks
To run a set of benchmarks (and any further tests), launch the command (on Matlab):

    `$ addpath(genpath(benchmarks/));addpath(genpath(include/));addpath(genpath(src/));`

Then to lauch the benchmarks:

 	$ do_all;
	$ result
    




