# CRM-CFP 
## CRM to solve the Convex Feasibility Problem

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> CRM-CFP

CRM-CFP concerns numerical results presented in [[Behling2020]](#1) and [[Arefidamghani2020]](#2) where the Circumcentered-Reflection Method (CRM) was used to solve the Convex Feasibility Problem (CFP) of finding a common point to the nonempty intersection of closed and convex sets.


It is authored by Luiz-Rafael Santos in co-authorship with Reza Arefidamghani, Roger Behling, Yunier Bello-Cruz and Alfredo N. Iusem. 

## How to use CRM-CFP

To (locally) reproduce this project, do the following:

0. Clone the respository locally

1. Open a Julia console and do:
   ```julia
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

## References


<a id="1">[Behling2020]</a>  R. Behling, J.-Y. Bello-Cruz, e L.-R. Santos, “On the Circumcentered-Reflection Method for the Convex Feasibility Problem”, Numer. Algorithms, jul. 2020, doi: [10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6).

<a id="2">[Arefidamghani2020]</a>  R. Arefidamghani, R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, e L.-R. Santos, 
“The circumcentered-reflection method achieves better rates than alternating projections”, 
arXiv:2007.14466, ago. 2020, (online).  [https://arxiv.org/abs/2007.14466].
