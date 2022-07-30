# Bi-Objective Simplex Algorithm

#### Aim

This code is an academic and minimalist implementation of the simplex algorithm for two (linear) objectives. It has been done for teaching purposes.

This is an extended version of the single objective simplex algorithm's implementation described in chapter 5 of [Julia programming for Operations Research (Changhyun Kwon)](http://www.chkwon.net/julia/). 

The version implemented corresponds to the parametric simplex described in chapter 6 of [Multicriteria Optimization (Matthias Ehrgott)](http://https://link.springer.com/book/10.1007/3-540-27659-9).

The code is compliant with julia 1.x, stable but must be polished (comments are partially written in french).

#### Example of a run

```
julia> c1 = [ 3;  1; 0; 0]
julia> c2 = [-1; -2; 0; 0]
julia> T  = [ 0  1  1  0 ;  3 -1  0  1 ]
julia> d = [3; 6]


julia> c1 = Array{Float64}(c1)
julia> c2 = Array{Float64}(c2)
julia> T  = Array{Float64}(T)
julia> d  = Array{Float64}(d)

julia> @time  parametric_simplex_method(c1, c2, T, d)
------+----------------------------+-------
      | -3.00  -1.00   0.00   0.00 |   0.00
------+----------------------------+-------
x[ 3] |  0.00   1.00   1.00   0.00 |   3.00
x[ 4] |  3.00  -1.00   0.00   1.00 |   6.00
------+----------------------------+-------
 
------+----------------------------+-------
      | -3.00  -1.00   0.00   0.00 |   0.00
      |  1.00   2.00   0.00   0.00 |   0.00
------+----------------------------+-------
x[ 3] |  0.00   1.00   1.00   0.00 |   3.00
x[ 4] |  3.00  -1.00   0.00   1.00 |   6.00
------+----------------------------+-------
Jc = Any[1, 2]
j=1lambda[j]=0.25
j=2lambda[j]=0.6666666666666666
[0.25 0.6666666666666666]
2
1 3
Pivoting: entering = x_2, exiting = x_3
------+----------------------------+-------
      | -3.00   0.00   1.00   0.00 |   3.00
      |  1.00   0.00  -2.00   0.00 |  -6.00
------+----------------------------+-------
x[ 2] |  0.00   1.00   1.00   0.00 |   3.00
x[ 4] |  3.00   0.00   1.00   1.00 |   9.00
------+----------------------------+-------
Jc = Any[1]
j=1lambda[j]=0.25
[0.25]
1
2 4
Pivoting: entering = x_1, exiting = x_4
------+----------------------------+-------
      |  0.00   0.00   2.00   1.00 |  12.00
      |  0.00   0.00  -2.33  -0.33 |  -9.00
------+----------------------------+-------
x[ 2] |  0.00   1.00   1.00   0.00 |   3.00
x[ 1] |  1.00   0.00   0.33   0.33 |   3.00
------+----------------------------+-------
  4.977332 seconds (13.13 M allocations: 635.104 MiB, 6.83% gc time)

julia>
```
