using Test
using Random
include("rotations.jl")

N=1000
quaternions = from_random(N)

@testset "rotations.jl" begin

    # conversion between quaternion and euler angle
    @test begin
        euler = as_euler_angle(quaternions)
        quat = from_euler_angle(euler)
        isClose(quaternions, quat, 10^(-13), 10^(-13))
    end

    # conversion between quaternion and axis angle
    @test begin
        axis = as_axis_angle(quaternions)
        quat = from_axis_angle(axis)
        isClose(quaternions, quat, 10^(-13), 10^(-13))
    end

    # conversion between quaternion and rodrigues frank
    @test begin
        rodr = as_rodriguesfrank(quaternions)
        quat = from_rodriguesfrank(rodr)
        isClose(quaternions, quat, 10^(-13), 10^(-13))
    end

    # conversion between quaternion and homochoric
    @test begin
        homo = as_homochoric(quaternions)
        quat = from_homochoric(homo)
        isClose(quaternions, quat, 10^(-13), 10^(-13))
    end

    # conversion between quaternion and rotation matrix
    @test begin
        matr = as_matrix(quaternions)
        quat = from_matrix(matr)
        isClose(quaternions, quat, 10^(-13), 10^(-13))
    end

end



    # # conversion between eulerAngle and quaternion
    # @test begin
    #             result = true
    #             for q in quaternions
    #                 e = as_euler_angle(q)
    #                 q1 = from_euler_angle(e)
    #                 if !isClose(q,q1,10^(-13),10^(-13))
    #                     result = false
    #                     break
    #                 end
    #             end
    #         result
    #     end
    #
    # # conversion between rotationMatrix and quaternion
    # @test begin
    #             result = true
    #             for q in quaternions
    #                 m = as_matrix(q)
    #                 q1 = from_matrix(m)
    #                 if !isClose(q,q1,10^(-13),10^(-13))
    #                     result = false
    #                     break
    #                 end
    #             end
    #         result
    #     end
    #
    # # conversion between axisAngle and quaternion
    # @test begin
    #             result = true
    #             for q in quaternions
    #                 a = as_axis_angle(q)
    #                 q1 = from_axis_angle(a)
    #                 if !isClose(q,q1,10^(-13),10^(-13))
    #                     result = false
    #                     break
    #                 end
    #             end
    #         result
    #     end
    #
    # # conversion between rodriguesFrank and quaternion (not yet implemented)
    # @test begin
    #             result = true
    #             for q in quaternions
    #                 r = as_rodriguesfrank(q)
    #                 q1 = from_rodriguesfrank(r)
    #                 if !isClose(q,q1,10^(-13),10^(-13))
    #                     result = false
    #                     break
    #                 end
    #             end
    #         result
    #     end
