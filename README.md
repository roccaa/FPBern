# FPBern
## Description
`FPBern` is a tool for computation of roundoff error bounds using Bernstein expansion.

If the roundoff error of a program implementing a polynomail function f(x)  is given by r(x,e) s.t.:

			r(x,e) = l(x,e)+h(x,e)

with `l(x,e)` the part of `r(x,e)` linear w.r.t `e` and `h(x,e)` the part  of `r(x,e)` non-linear in `e`.

Then `FPBern` gives an upper bound of `|l(x,e)|` for `x` laying inside a box, and `e` enclose by a given epsilon. 

Thus, the semantic of the handled programs is:

- polynomial functions taken into isolation

ex: x1^3 + (3/4)*x1*x2^2
- inputs laying inside a box

ex: x1 = [-1 1] and x2 = [-2 2]

## Sub modules

FPBern is divide in two modules:

- FPBern(a): computes the roundoff error with a C++ implementation.

- FPBern(b): computes the roundoff error with a Matlab2015a implementation.

Please refer to the the specific README of each module for details on the use and installation. 





