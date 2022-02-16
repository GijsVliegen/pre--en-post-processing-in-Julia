using LinearAlgebra

P = 1

struct quaternion
    i
    j
    k
    angle
end

struct eulerAngle
    phi1
    PHI
    phi2
end

struct rotationMatrix
    matrix# ::Array{Int64, 2}(undef, 2, 2)
end

struct axisAnglePair
    n
    omega
end

struct rodriguesFrank
end

struct homochoric
end

function normalize(rot ::quaternion)
    abs_val = sqrt(rot.angle^2 + rot.i^2 + rot.j^2 + rot.k^2)
    return quaternion(rot.angle/abs_val, rot.i/abs_val, rot.j/abs_val, rot.k/abs_val)
end

function getComponents(rot ::quaternion)
    return [rot.i, rot.j, rot.k, rot.angle]
end

function toQuaternion(rot ::quaternion)
    return rot
end

function copy(rot ::quaternion)
    return quaternion(rot.i, rot.j, rot.k, rot.angle)

end
#in tegenstelling tot in damask wordt er niet rekening gehouden met NaN values
#dit doen we omdat het niet is ingebouwd in isapprox()

function isClose(first ::quaternion, other ::quaternion, rtol, atol, nanEquals = true)
    fcom = getComponents(toQuaternion(first))
    scom = getComponents(toQuaternion(other))
    for i in 1:4
        if !isapprox(fcom[i], scom[i], rtol = rtol , atol = atol, nans = nanEquals)
            return false
        end
    end
    return true
end

#wat volgt zijn algoritmes voor het omzetten van rotaties naar een bepaalde representatie
#enkel de directe omzettingen zijn gegeven, andere omzettingen kunnen opgebouwd worden
#volgens de tabel pagina 11 van de paper "Consistent representations of and conversions
#between 3D rotations"


#voor de gebruiksvriendelijkheid kan het idee van bovenstaande code ook gebruikt worden.

function eulerAngleToRotationMatrix(rotation ::eulerAngle)
    c1 = cos(rotation.phi1)
    s1 = sin(rotation.phi1)
    c2 = cos(rotation.phi2)
    s2 = sin(rotation.phi2)
    C = cos(rotation.PHI)
    S = sin(rotation.PHI)
    matrix = [c1*c2-s1*C*s2 s1*c2+c1*C*s2 S*s2; -c1*s2-s1*C*c2 -s1*s2+c1*C*c2 S*c2; s1*S -c1*S C]
    return rotationMatrix(matrix)
end

function eulerAngleToAxisAngle(rotation ::eulerAngle)::axisAnglePair
    t = tan(rotation.PHI/2)
    sigma = (rotation.phi1 + rotation.phi2)/2.0
    delta = (rotation.phi1 - rotation.phi2)/2.0
    tau = sqrt(t^2 + sin(sigma)^2)
    alpha = 2*atan(tau/(cos(sigma)))
    if (alpha > pi)
        alpha = 2*pi - alpha
    end
    return axisAnglePair([P/t*cos(delta) P/t*sin(delta) P/t*sin(sigma)], a)
end

function eulerAngleToRodriguesFrank(rotation ::eulerAngle)
    #EA naar AA en dan AA naar RF
    axisAngle = eulerAngleToAxisAngle(rotation)
    return rodriguesFrankVector(axisAngle.n, tan(axisAngle.w/2))
end

function eulerAngleToQuaternion(rotation ::eulerAngle)
    sigma = (rotation.phi1 + rotation.phi2)/2.0
    delta = (rotation.phi1 - rotation.phi2)/2.0
    c = cos(rotation.PHI/2)
    s = sin(rotation.PHI/2)
    return quaternion(c*cos(sigma), -P*s*cos(delta), -P*s*sin(delta), -P*c*sin(sigma))
end

function rotationMatrixToEulerAngle(rotation ::rotationMatrix)
    a_33 = rotation.matrix[3, 3]
    if abs(a_33) != 1
        zeta = 1/sqrt(1-a_33^2)
        phi1 = atan(rotation.matrix[3, 1]*zeta, -rotation.matrix[3, 2]*zeta) 
        PHI = acos(a_33)
        phi2 = atan(rotation.matrix[1, 3]*zeta, rotation.matrix[2, 3]*zeta)
        return eulerAngle(phi1, PHI, phi2)
    else
        phi1 = atan(rotation.matrix[1, 2], rotation.matrix[1, 1])
        PHI = pi/2*(1 - a_33)
        return eulerAngle(phi1, PHI, 0)
    end
end

function rotationMatrixToAxisAngle(rotation ::rotationMatrix)
    omega = acos((tr(rotation.matrix) - 1)/2)
    #is voor een volgende keer
end

function rotationMatrixToQuaternion(rotation ::rotationMatrix)
    a = rotation.matrix
    q0 = 1/2*sqrt(1 + a[1, 1] + a[2, 2] + a[3,3])
    q1 = P/2*sqrt(1 + a[1, 1] - a[2, 2] - a[3,3])
    q2 = P/2*sqrt(1 - a[1, 1] + a[2, 2] - a[3,3])
    q3 = P/2*sqrt(1 - a[1, 1] - a[2, 2] + a[3,3])
    if (a[3,2] < a[2, 3])
        q1 = -q1
    end
    if (a[1, 3] < a[3, 1])
        q2 = -q2
    end
    if (a[2, 1] < a[1, 2])
        q3 = -q3
    end
    return normalize(quaternion(q0, q1, q2, q3))
end

function axisAngleToRotationMatrix(rotation ::axisAnglePair)
    n = rotation.n
    c = cos(rotation.omega)
    s = sin(rotation.omega)
    a = [c+(1-c)*n[1]^2 (1-c)*n[1]*n[2]+s*n[3] (1-c)*n[1]*n[3]-s*n[2]; 
    (1-c)*n[1]*n[2]-s*n[3] c+(1-c)*n[2]^2 (1-c)*n[2]*n[3]+s*n[1];
    (1-c)*n[1]*n[3]+s*n[2] (1-c)*n[2]*n[3]-sn[1] c+(1-c)*n[3]^2]
    if (P == 1)
        a = transpose(a)
    end
    return a
end

function axisAngleToRodriguesFrank(rotation ::axisAnglePair)
    #wat als omega = pi???
    return rodriguesFrank(rotation.n*tan(rotation.omega/2))
end

function axisAngleToQuaternion(rotation ::axisAnglePair)
    n = rotation.n*sin(rotation.omega/2)
    return quaternion(n[1], n[2], n[3], cos(rotation.omega/2))
end

function axisAngleToHomochoric(rotation ::axisAnglePair)
    f = (3/4(rotation.omega - sin(rotation.omega)))^(1/3)
    return homochoric(rotation.n*f)
end

function rodriguesFrankToAxisAngle(rotation ::rodriguesFrank)
    rho = abs(rotation.rho)
    return axisAnglePair(rotation.rho/rho, 2*arctan(rho))
end

function rodriguesFrankToHomochoric(rotation ::rodriguesFrank)
    rho = abs(rotation.rho)
    if (rho == 0) 
        return homochoric([0, 0, 0])
    end
    #is het nodig om rho == oneindig te behandelen
    w = 2*atan()
    f = 3(rotation.omega - sin(rotation.omega))/4 
end

function quaternionToEulerAngle()
end

function quaternionToRotationMatrix()
end

function quaternionToAxisAngle()
end

function quaternionToRodriguesFrank()
end

function quaternionToHomochoric()
end

function homochoricToAxisAngle()
end

function homochoricToCubochoric()
end

function qubochoricToHomochoric()
end