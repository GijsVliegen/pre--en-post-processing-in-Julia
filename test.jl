using Test
using Random
include("rotations.jl")

# quOne = quaternion(1, 1, 1, 0)
# quTwo = quaternion(1+1e-1, 1, 1, 0)
# if isClose(quOne, quTwo, 0.01, 0.01)
#     println("yeet")
# else
#     println("niet yeet")
# end
#
# eA = eulerAngle(1, 1, 1)
# rotMatrix = eulerAngleToRotationMatrix(eA)
# print(rotMatrix)

# generate a array of N uniformly distributed random quaternions: http://planning.cs.uiuc.edu/node198.html
N=1000
quaternions = quaternion[]
for i = 1:N
    u1 = rand()
    u2 = rand()
    u3 = rand()
    h = quaternion(sqrt(1-u1)*sin(2*pi*u2), sqrt(1-u1)*cos(2*pi*u2), sqrt(u1)*sin(2*pi*u3), sqrt(u1)*cos(2*pi*u3))
    push!(quaternions, copy(h))
end

eulerAngles = eulerAngle[]
for i = 1:N

end


@testset "rotations.jl" begin
    # conversion between eulerAngle and quaternion
    @test begin
                result = true
                for q in quaternions
                    e = toEulerAngle(q)
                    q1 = toQuaternion(e)
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
                    m = toRotationMatrix(q)
                    q1 = toQuaternion(m)
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
                    a = toAxisAngle(q)
                    q1 = toQuaternion(a)
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
                    r = toRodriguesFrank(q)
                    q1 = toQuaternion(r)
                    if !isClose(q,q1,10^(-13),10^(-13))
                        result = false
                        break
                    end
                end
            result
        end

end
