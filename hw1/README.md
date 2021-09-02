## Homework 1, Questions 1,2,7,8 of fundamentals of Engineering Numerical Analysis 
### Question 1 
(a) The interpolatd value at x=0.9 is 0.0235
(b) As you can see below results for the lagrange curve are too far from the
smooth real curve since it is a 20th order polynomial which meets all the 21
equaly distanced points on the smooth Runge's curve"

<img src="./1_b.png?raw=true" width="500">

### Question 2 
I thought about the chain law and came up with the solution which you can see
below.
[fig 1_b](https://latex.codecogs.com/gif.latex?%5Csum_%7Bk%3D0%7D%5E%7Bn%7D%5Cleft%20%28%20y_%7Bk%7D%5Csum_%7Bj%3D0%2Cj%5Cneq%20k%7D%5E%7Bn%7D%5Cleft%20%28%20%5Cdisplaystyle%20%5Cfrac%7B1%7D%7Bx_%7Bk%7D-x%7Bj%7D%7D%5Cprod_%7Bi%3D0%2Ci%5Cneq%20j%2C%20i%5Cneq%20k%7D%5E%7Bn%7D%5Cfrac%7Bx-x_%7Bi%7D%7D%7Bx_%7Bk%7D-x_%7Bi%7D%7D%20%5Cright%20%29%20%5Cright%20%29)

### Question 7 
(a,b,c) Becouse of large error in results I conclude that the global \
interpolation of Lagrange is not suitable for exterpolation \
and is not capable of capturing smooth curves beyond the data limits.

<img src="./7.png?raw=true" width="500">

### Question 8 
(a) The prediction by Lagrange polynomial is way off and does not fallow the
trend of the data points at all.

<img src="./8_a.png?raw=true" width="500">
(b,c) cubic spline interpolation is much closer to the real data than the
lagrange polynomial interpolation as expected. cubic spline is a local method
for interpolation afterall.

<img src="./8_bc.png?raw=true" width="500">
