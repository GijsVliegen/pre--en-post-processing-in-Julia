include("rotations.jl")
using BenchmarkTools #installatie vanuit https://github.com/JuliaCI/BenchmarkTools.jl

function operation(a)
    b = as_matrix(a)
    c = from_matrix(b)
end

function time(a)
    @time operation(a)
    #@benchmark operation(a) 
end
function timing(n)
    a = from_random(n, n, Float64)
    #println(a)
    time(a)
end

function timingB(a)
    b = as_matrix(a)
    c = from_matrix(b)
    return c
end
