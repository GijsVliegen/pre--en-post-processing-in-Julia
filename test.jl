using Test
using Random
include("rotations.jl")

# generate a array of N uniformly distributed random quaternions: http://planning.cs.uiuc.edu/node198.html
N=1000
quaternions = from_random(N)

@testset "rotations.jl" begin
    # conversion between eulerAngle and quaternion
    @test begin
                result = true
                for q in quaternions
                    e = as_euler_angle(q)
                    q1 = from_Euler_angles(e)
                    if !isClose(q,q1,10^(-13),10^(-13))
                        result = false
                        break
                    end
                end
            result
        end

    # conversion between rotationMatrix and quaternion
    @test begin
                result = true
                for q in quaternions
                    m = as_matrix(q)
                    q1 = from_matrix(m)
                    if !isClose(q,q1,10^(-13),10^(-13))
                        result = false
                        break
                    end
                end
            result
        end

    # conversion between axisAngle and quaternion
    @test begin
                result = true
                for q in quaternions
                    a = as_axis_angle(q)
                    q1 = from_axis_angle(a)
                    if !isClose(q,q1,10^(-13),10^(-13))
                        result = false
                        break
                    end
                end
            result
        end

    # conversion between rodriguesFrank and quaternion (not yet implemented)
    @test begin
                result = true
                for q in quaternions
                    r = as_rodriguesfrank(q)
                    q1 = from_rodriguesfrank(r)
                    if !isClose(q,q1,10^(-13),10^(-13))
                        result = false
                        break
                    end
                end
            result
        end

end
