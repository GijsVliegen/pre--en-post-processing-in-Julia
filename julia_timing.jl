include("rotations.jl")
using BenchmarkTools #installatie vanuit https://github.com/JuliaCI/BenchmarkTools.jl

function operationOm(a)
    b = as_matrix(a)
    c = from_matrix(b)
end

function operationEu(a)
    b = as_euler_angle(a)
    c = from_Euler_angles(b)
end

n = 10000
a = from_random(n, n, Float64)
@btime operationOm($a)

