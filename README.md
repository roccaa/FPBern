# FPBern
## Description
`FPBern` is a tool for computation of roundoff error bounds using Bernstein expansion.

We consider the roundoff error r of a program implementing a polynomial function f(x):

			r(x,e) = l(x, e) + h(x, e)

with `l(x, e)` the part of `r(x, e)` linear w.r.t. `e` and `h(x, e)` the part  of `r(x, e)` non-linear w.r.t. `e`.

`FPBern` computes an upper bound of `|l(x, e)|` with `x` lying in a box, and `e` being bounded by a given machine `epsilon`.  

For instance, one can consider f(x) = x1^3 + (3/4) x1 x2^2, with x1 in [-1, 1] and x2 in [-2, 2].

## Sub modules

FPBern is divided in two independant modules:

- FPBern(a): computes the roundoff error with a C++ implementation.

- FPBern(b): computes the roundoff error with a Matlab2015a implementation.

Please refer to the the specific README of each module for details on the installation and usage.
