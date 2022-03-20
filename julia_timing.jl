include("rotations.jl")
using BenchmarkTools #installatie vanuit https://github.com/JuliaCI/BenchmarkTools.jl

function operationOm(a)
    b = as_matrix(a)
    c = from_matrix(b)
end

function operationEu(a)
    b = as_euler_angle(a)
    c = from_euler_angle(b)
end
function operationAx(a)
    b = as_axis_angle(a)
    c = from_axis_angle(b)
end
function operationHo(a)
    b = as_homochoric(a)
    c = from_homochoric(b)
end
function operationRo(a)
    b = as_rodriguesfrank(a)
    c = from_rodriguesfrank(b)
end

n = 200000
a = from_random(n, n, Float64)
println("Om:")
@btime operationOm($a)
println("Eu:")
@btime operationEu($a)
println("Ax:")
@btime operationAx($a)
println("Ho:")
@btime operationHo($a)
println("Ro:")
@btime operationRo($a)

