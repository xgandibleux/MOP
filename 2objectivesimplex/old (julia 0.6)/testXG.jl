include("simplex_XG.jl")
using SimplexMethod

include("simplexparam_XG.jl")
#using ParametricSimplexMethod

#c1 = [-3; -2; -1; -5; 0; 0; 0]
#c2 = [-3; -2; -1; -5; 0; 0; 0]
#T = [7 3 4 1 1 0 0 ;
#     2 1 1 5 0 1 0 ;
#     1 4 5 2 0 0 1 ]
#d = [7; 3; 8]

c1 = [ 3;  1; 0; 0]
c2 = [-1; -2; 0; 0]
T  = [ 0  1  1  0 ;
       3 -1  0  1 ]
d = [3; 6]

c1 = Array{Float64}(c1)
c2 = Array{Float64}(c2)
T  = Array{Float64}(T)
d  = Array{Float64}(d)

#using MathProgBase, GLPK
#@time sol = linprog(c, A, '=', b, GLPKSolverLP) #GurobiSolver(Method=0))

#include("search_bfs.jl")
#@time opt_x1, obj1 = searchBFS(c, A, b)

println("Objective 1: ")
@time opt_x1, obj1 = simplex_method(c1, T, d)
println(" ")
println("Objective 2: ")
@time opt_x2, obj2 = simplex_method(c2, T, d)
println(" ")
println("Objective 1&2: ")
@time  parametric_simplex_method(c1, c2, T, d)


#println()
#println("obj by MathProgBase.linprog  = ", sol.objval)
#println("obj by search_extreme_points = ", obj1)
#println("obj by simplex_method        = ", obj2)
#println()
#println("x*  by MathProgBase.linprog  = ", sol.sol)
#println("x*  by search_extreme_points = ", opt_x1)
#println("x*  by simplex_method        = ", opt_x2)
#println()
