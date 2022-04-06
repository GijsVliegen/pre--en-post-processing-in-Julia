using Test
using Random
include("rotations.jl")
std_err_rate = 10^(-8)
N=1000
quaternions = from_random(N)

@testset "rotations.jl, n = 1000" begin

    #------------euler angle---------------

    # conversion between quaternion and euler angle
    @test begin
        euler = as_Euler_angles(quaternions)
        quat = from_Euler_angles(euler)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #euler_angle_degrees
    @test begin
        euler = as_Euler_angles(quaternions)
        euler = euler*180/pi
        quat = from_Euler_angles(euler, true)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #------------axis angle---------------

    # conversion between quaternion and axis angle
    @test begin
        axis = as_axis_angle(quaternions)
        quat = from_axis_angle(axis)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #normalize
    @test begin
        axis = as_axis_angle(quaternions)
        for i in 1:4:length(axis)
            axis[i:i+2] = axis[i:i+2] * 2.71
        end
        quat = from_axis_angle(axis, false, true)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #degrees
    @test begin
        axis = as_axis_angle(quaternions, true)
        quat = from_axis_angle(axis, true)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #pair
    @test begin
        axis = as_axis_angle(quaternions, false, true)
        together = zeros(4, length(quaternions))
        for i in 1:length(axis[2])
            one_ax_angle = vcat(axis[1][i*3-2:i*3], axis[2][i])
            together[i*4-3:i*4] = one_ax_angle
        end
        quat = from_axis_angle(together, false)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end


    #------------basis---------------
    """
    @test begin
        matr = as_basis(quaternions)
        quat = from_basis(matr)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end"""

    #------------rodriguez frank---------------

    # conversion between quaternion and rodrigues frank
    @test begin
        rodr = as_Rodrigues_vector(quaternions)
        quat = from_Rodrigues_vector(rodr)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #normalize
    @test begin
        rodr = as_Rodrigues_vector(quaternions)
        for i in 1:4:length(rodr)
            rodr[i:i+2] = rodr[i:i+2] * 2.71
        end
        quat = from_Rodrigues_vector(rodr)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #compact
    @test begin
        rodr = as_Rodrigues_vector(quaternions, true)
        together = zeros(4, length(quaternions))
        for i in 1:length(quaternions)
            norme = norm(rodr[i*3-2:i*3])
            together[i*4-3:i*4] = vcat(rodr[i*3-2:i*3]/norme, norme)
        end
        quat = from_Rodrigues_vector(together)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end


    #------------homochoric---------------

    # conversion between quaternion and homochoric
    @test begin
        homo = as_homochoric(quaternions)
        quat = from_homochoric(homo)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #------------matrix---------------

    # conversion between quaternion and rotation matrix
    @test begin
        matr = as_matrix(quaternions)
        quat = from_matrix(matr)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

end

N = 1
quaternions = from_random(N)

@testset "rotations.jl, n = 1" begin

    #------------euler angle---------------

    # conversion between quaternion and euler angle
    @test begin
        euler = as_Euler_angles(quaternions)
        quat = from_Euler_angles(euler)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #euler_angle_degrees
    @test begin
        euler = as_Euler_angles(quaternions)
        euler = euler*180/pi
        quat = from_Euler_angles(euler, true)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #------------axis angle---------------

    # conversion between quaternion and axis angle
    @test begin
        axis = as_axis_angle(quaternions)
        quat = from_axis_angle(axis)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #normalize
    @test begin
        axis = as_axis_angle(quaternions)
        for i in 1:4:length(axis)
            axis[i:i+2] = axis[i:i+2] * 2.71
        end
        quat = from_axis_angle(axis, false, true)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #degrees
    @test begin
        axis = as_axis_angle(quaternions, true)
        quat = from_axis_angle(axis, true)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #pair
    @test begin
        axis = as_axis_angle(quaternions, false, true)
        together = zeros(4, length(quaternions))
        for i in 1:length(axis[2])
            one_ax_angle = vcat(axis[1][i*3-2:i*3], axis[2][i])
            together[i*4-3:i*4] = one_ax_angle
        end
        quat = from_axis_angle(together, false)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end


    #------------basis---------------
    """
    @test begin
        matr = as_basis(quaternions)
        quat = from_basis(matr)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end"""

    #------------rodriguez frank---------------

    # conversion between quaternion and rodrigues frank
    @test begin
        rodr = as_Rodrigues_vector(quaternions)
        quat = from_Rodrigues_vector(rodr)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #normalize
    @test begin
        rodr = as_Rodrigues_vector(quaternions)
        for i in 1:4:length(rodr)
            rodr[i:i+2] = rodr[i:i+2] * 2.71
        end
        quat = from_Rodrigues_vector(rodr)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #compact
    @test begin
        rodr = as_Rodrigues_vector(quaternions, true)
        together = zeros(4, length(quaternions))
        for i in 1:length(quaternions)
            norme = norm(rodr[i*3-2:i*3])
            together[i*4-3:i*4] = vcat(rodr[i*3-2:i*3]/norme, norme)
        end
        quat = from_Rodrigues_vector(together)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end


    #------------homochoric---------------

    # conversion between quaternion and homochoric
    @test begin
        homo = as_homochoric(quaternions)
        quat = from_homochoric(homo)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
    end

    #------------matrix---------------

    # conversion between quaternion and rotation matrix
    @test begin
        matr = as_matrix(quaternions)
        quat = from_matrix(matr)
        isClose(quaternions, quat, std_err_rate, std_err_rate)
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
