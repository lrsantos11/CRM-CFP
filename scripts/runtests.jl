using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))
using CSV, DataFrames
using BenchmarkProfiles, LaTeXStrings, Plots
pgfplotsx()
include(scriptsdir("plots_util.jl"))
##
"""
* Use the following to run the tests and generate the pictures that appear in [^Behling2020].
```
julia> include(scriptsdir("runtests_BBS20.jl"))
```
* Use the following to run the tests and generate the pictures that appear in  [^Arefidamghani21].
```
julia> include(scriptsdir("runtests_ABBIS20.jl"))
```
* Use the following to run the tests and generate the pictures that appear in [^Araujo21].
```
julia> include(scriptsdir("runtests_AABBIS20.jl"))
```

* References
[^Araujo21] G. Araújo, R. Arefidamghani, R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, e L.-R. Santos, 

[^Arefidamghani21] R. Arefidamghani, R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, e L.-R. Santos, 
“The circumcentered-reflection method achieves better rates than alternating projections”, 
Comp Optim Appl, [Online], apr. 2021, . doi: [10.1007/s10589-021-00275-6](https://doi.org/10.1007/s10589-021-00275-6)

[^Behling2020] R. Behling, J.-Y. Bello-Cruz, e L.-R. Santos, “On the Circumcentered-Reflection Method 
for the Convex Feasibility Problem”, Numer. Algorithms, jul. 2020, 
doi: [10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6). 
"""

# include(scriptsdir("runtests_BBS20.jl"))
# include(scriptsdir("runtests_ABBIS20.jl"))
# include(scriptsdir("runtests_AABBIS21.jl"))