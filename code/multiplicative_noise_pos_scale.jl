
# Import 3rd party packages

using DifferentialEquations, StochasticDiffEq, LinearAlgebra, Calculus, Discretizers, KernelDensity, Distributions
using JLD2, FileIO, DataFrames, LaTeXStrings;
using DelimitedFiles,CSV, Query;
using Plots;

# v(x) = a1*(x-1.5)^4/4 - b1*(x-1.5)^3/3 - c1*(x-1.5)^2/2 - d1
# v'(x) = a1*(x-1.5)^3 - b1*(x-1.5)^2 - c1*(x-1.5) - d1
# f(x) = -v'(x)

#######################################################
# f(x) = -a1*(x-1.5)^3 + b1*(x-1.5)^2 + c1*(x-1.5) + d1
# g(x) = α(x-1.5)
#######################################################

JSON.parse(document.getElementById('jupyter-config-data').textContent).token

a1=0.7
b1=0.0
c1=1.0
d1=0.1

v(x) = (a1*(x-1.5)^4)/4 - (b1*(x-1.5)^3)/3 - (c1*(x-1.5)^2)/2 - d1*(x-1.5)

plot(v,-1,4)

# Compute the stable stationary points and then unstable stationary points

# stable=0.35832
# stable=2.7424
# unstable=1.39928

# g(x) = sqrt(f(x))

test(x) = exp(-log(abs(-a1*x^3+b1*x^2+c1*x+d1)) + 2*x - 1)
plot(test)

test(x) = exp(-log(abs(-a1*(x-1.5)^3+b1*(x-1.5)^2+c1*(x-1.5)+d1)) + 2*(x-1.5))
plot(test)

land_test(x) = log(abs(-a1*x^3+b1*x^2+c1*x+d1)) - 2*x
plot(land_test)

land_test(x) = log(abs(-a1*(x-1.5)^3+b1*(x-1.5)^2+c1*(x-1.5)+d1)) - 2*(x-1.5)
plot(land_test,0,5)

##################################################
# ANALYTICAL STEADY STATE PROBABILITY DISTRIBUTION
##################################################

α=0.04
psx(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )

plot(psx)

α=0.001
psx(x) = exp( (-2log(abs(x-1.5)))
    - big((a1*x^2)/(α^2)) + big((2*x*(1.5*a1+b1))/(α^2)) + big(2*c1*log(abs(x-1.5))/α^2) - big((2*d1)/((α^2)*(x-1.5))) )

plot(psx)

d1=0.0
################################################################
# ADJUSTING NOISE LEVELS AND EFFECT ON NORMALISED SSD WITH D=0.0
################################################################

# Solution is undefined for x=[1.47:1.50]

α=0.04
psx_1(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )
p1=plot(psx_1,0,4,legend=false)

α=0.1
psx_2(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )
p2=plot(psx_2,0,4,legend=false)

α=0.5
psx_3(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )
p3=plot(psx_3,0,4,legend=false)

α=0.99
psx_4(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )
p4=plot(psx_4,0,4,legend=false)

plot(p1,p2,p3,p4,layout=(4,1))

# THIS CREATES NICE CONTINUOUS PLOTS SO THE ISSUE MUST BE WITH THE LAST TERM
# IF X EVALUATES TO ZERO THEN WE WILL BE DIVIDING BY ZERO AND THEREFORE SOLUTION IS UNDEFINED AT THIS POINT

# WHEN NOISE IS HIGH THEN SOLUTION DIPS AT UNSTABLE STATIONARY POINT

#####################################################
# ADJUSTING NOISE LEVELS AND EFFECT ON NORMALISED SSD
#####################################################

# Solution is undefined for x=[1.47:1.50]

x_left = collect(0.0:0.01:1.46)
x_left1 = collect(-0.5:0.01:1.4)
x_left2 = collect(-1.0:0.01:1.2)
x_right = collect(1.51:0.01:4.0)
x_right1 = collect(1.51:0.01:4.0)

α=0.04
psx_1(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )
psx_norm = psx_1.(x_left)
psx_norm = psx_norm./(0.01*sum(psx_norm))
plot(x_left,psx_norm)
psx_norm = psx_1.(x_right)
psx_norm = psx_norm./(0.01*sum(psx_norm))
p1=plot!(x_right,psx_norm)

α=0.1
psx_2(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )
psx_norm = psx_2.(x_left)
psx_norm = psx_norm./(0.01*sum(psx_norm))
plot(x_left,psx_norm)
psx_norm = psx_2.(x_right)
psx_norm = psx_norm./(0.01*sum(psx_norm))
p2=plot!(x_right,psx_norm)

α=0.5
psx_3(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )
psx_norm = psx_3.(x_left1)
psx_norm = psx_norm./(0.01*sum(psx_norm))
plot(x_left1,psx_norm)
psx_norm = psx_3.(x_right)
psx_norm = psx_norm./(0.01*sum(psx_norm))
p3=plot!(x_right,psx_norm)

α=1.0
psx_4(x) = exp( (-2log(abs(x-1.5)))
    - ((a1*x^2)/(α^2)) + ((2*x*(1.5*a1+b1))/(α^2)) + (2*c1*log(abs(x-1.5))/α^2) - ((2*d1)/((α^2)*(x-1.5))) )
psx_norm = psx_4.(x_left2)
psx_norm = psx_norm./(0.01*sum(psx_norm))
plot(x_left2,psx_norm)
psx_norm = psx_4.(x_right1)
psx_norm = psx_norm./(0.01*sum(psx_norm))
p4=plot!(x_right1,psx_norm)

all_psx = plot(p1,p2,p3,p4,layout=(4,1),size=(600,800))

##################################################
# LANDSCAPE RECOVER FROM U = -LOG(Ps(x))
##################################################

# Computed analytically

α=0.1
land(x) = 2log(abs(x-1.5)) + ((a1*x^2)/(α^2)) - ((2*x*(1.5*a1+b1))/(α^2)) - (2*c1*log(abs(x-1.5))/α^2) + ((2*d1)/((α^2)*(x-1.5)))

# plot(land)
plot(land, -1,1.46)
plot!(land, 1.9,5)









########## TESTING ##########
#############################

α=0.2

# excl integral evaluated at x=1
psx(x) = exp(-2log(abs(x)) - ((a1*x^2)/(α^2)) + ((2*x*(4.5*a1+b1))/(α^2)) - ((2*(6.75*a1-3*b1+c1))*log(abs(x))/(α^2)) - ((2*(3.375*a1+2.25*b1-1.5*c1+d1))/(α^2*x)))
plot(psx)

x = collect(0.01:0.01:5.0)
psx_norm = psx.(x)
psx_norm = psx_norm./(0.01*sum(psx_norm))
plot(x,psx_norm)



# eliminating b1 at beginning as it evaluates to 0

α=1.5

# excl integral evaluated at x=1
psx(x) = exp(-2log(x) + 2*( (-a1*x^2)/(2*α^2) + (4.5*a1*x)/α^2 - (6.75*a1*log(abs(x)))/α^2 - (3.375*a1-1.5+d1)/(α^2*x)))

plot(psx)

# incl integral evaluated at x=1
psx_1(x) = exp(-2log(x)
    + 2*( ((-a1*x^2)/(2*α^2) + (4.5*a1*x)/α^2 - (6.75*a1*log(abs(x)))/α^2 - (3.375*a1-1.5+d1)/(α^2*x))) -
(-a1/(2*α^2) + (4.5*a1)/α^2 - (3.375*a1-1.5+d1)/α^2 ))

plot(psx_1)

# Keeping the b1 term

α=0.2

t1=-a1
t2=(4.5*a1+b1)
t3=(6.75*a1-3*b1+c1)
t4=(3.375*a1+2.25*b1-1.5*c1+d1)

# excl integral evaluated at x=1
test(x) = exp(-2log(x) + 2*( (t1*x^2)/(2*α^2) + (t2*x)/α^2  - t3*log(abs(x))/α^2 - t4/(α^2*x)))
plot(test)

# incl integral evaluated at x=1
test_1(x) = exp(-2log(x) + 2*( (t1*x^2)/(2*α^2) + (t2*x)/α^2  - t3*log(abs(x))/α^2 - t4/(α^2*x))
- (-a1/(2*α^2) + t2/α^2 - t4/(α^2)))
plot(test_1)

#############################################################
# NORMALISED ANALYTICAL STEADY STATE PROBABILITY DISTRIBUTION
#############################################################

x = collect(0.01:0.01:5.0)
psx_norm = psx_1.(x)
psx_norm = psx_norm./(0.01*sum(psx_norm))
plot(x,psx_norm)

x = collect(0.01:0.01:5.0)
test_1_norm = test_1.(x)
test_1_norm = test_1_norm./(0.01*sum(test_1_norm))
plot(x,test_1_norm)



α=1.0

land_psx_1(x) = 2log(x) -2*((-a1*x^2)/(2*α^2) + (4.5*a1*x)/α^2 - (6.75*a1*log(abs(x)))/α^2 - (3.375*a1-1.5+d1)/(α^2*x))
+ (-a1/(2*α^2) + (4.5*a1)/α^2 - (3.375*a1-1.5+d1)/α^2)

plot(land_psx_1)

land_test_1(x) = 2log(x) - 2*((t1*x^2)/(2*α^2) + (t2*x)/α^2  - t3*log(abs(x))/α^2 - t4/(α^2*x)) + (-a1/(2*α^2) + t2/α^2 - t4/(α^2))

plot(land_test_1)
