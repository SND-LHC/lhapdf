unset key

cutoff = 0.83

theta(x) = x<0 ? 0 : 1
r(x) = x * theta(x-cutoff)
g(x) = (1-x/cutoff)*(1-theta(x-cutoff))
b(x) = x/cutoff*(1-theta(x-cutoff))

set xr [1e-8:1.0]
set yr [10:1e8]

set logscale xy
set palette model RGB functions r(gray), g(gray), b(gray)

p "accuracy.dat" using 1:2:4 with image
