include("rotations.jl")
using BenchmarkTools #installatie vanuit https://github.com/JuliaCI/BenchmarkTools.jl

function operationOm(a)
    b = as_matrix(a)
    c = from_matrix(b)
end

function operationEu(a)
    b = as_Euler_angles(a)
    c = from_Euler_angles(b)
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
    b = as_Rodrigues_vector(a)
    c = from_Rodrigues_vector(b)
end

function oneTest(a)
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
end

function testAll()
    for i in [(17: 22)...]
        n = 2^i
        a = from_random(n)
        println("voor n = ", n)
        oneTest(a)
    end
end

function testInPlaceVSNormalOnce(a)
    println("Eu:")
    @btime operationEu($a)
    println("Eu!:")
    @btime operationEu!($a)
end

function testInPlaceVSNormal()
    for i in [(1:16)...]
        n = 2^i
        a = from_random(n, n)
        println("voor n = ", n)
        testInPlaceVSNormalOnce(a)
    end
end

function testMultDim()
    a = from_random(10000)
    oneTest(a)
    a = from_random(10000, (10, 10, 10, 10))
    oneTest(a)
end

testAll();