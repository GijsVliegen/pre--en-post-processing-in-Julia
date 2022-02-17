using LinearAlgebra
import Base: copy

P = 1

struct quaternion
    i
    j
    k
    angle
end

struct eulerAngle
    #ZXZ structuur
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
    n
    f
end

struct homochoric
    h
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

    #nog eens nakijken
    axisAngle = eulerAngleToAxisAngle(rotation)
    return rodriguesFrank(axisAngle.n, tan(axisAngle.w/2))
end

function eulerAngleToQuaternion(rotation ::eulerAngle)
    sigma = (rotation.phi1 + rotation.phi2)/2.0
    delta = (rotation.phi1 - rotation.phi2)/2.0
    c = cos(rotation.PHI/2)
    s = sin(rotation.PHI/2)
    q0 = c*cos(sigma)
    if q0 < 0
        return quaternion(P*s*cos(delta), P*s*sin(delta), P*c*sin(sigma), q0)
    else
        return quaternion(-P*s*cos(delta), -P*s*sin(delta), -P*c*sin(sigma), q0)
    end
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
    #inspiratie voor de code: https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
    #dit algo komt niet uit de paper, die onduidelijk was, opletten voor omega = 0 of omega = pi
    #omega = 0 -> x y en z zijn arbitrair
    #omega = pi -> x y en z maken uit dus moeten berekend worden
    #handelen van singularities: https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/david.htm
    
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
    return rotationMatrix(a)
end

function axisAngleToRodriguesFrank(rotation ::axisAnglePair)
    #wat als omega = pi???
    f = tan(rotation.omega/2)
    return rodriguesFrank(rotation.n*f, f)
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
    rho = abs(rotation.n)
    return axisAnglePair(rotation.n/rho, 2*arctan(rho))
end

function rodriguesFrankToHomochoric(rotation ::rodriguesFrank)
    rho = abs(rotation.rho)
    if (rho == 0)
        return homochoric([0, 0, 0])
    end
    #is het nodig om rho == oneindig te behandelen
    w = 2*atan(rho)
    f = 3(w - sin(w))/4
    return homochoric(rotation.n*f^(1/3))
end

function quaternionToEulerAngle(rotation ::quaternion)
    #moet unit quaternion zijn
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    q03 = q0^2 + q3^2
    q12 = q1^2 + q2^2
    x = sqrt(q03 * q12)
    if (x == 0 && q12 == 0)
        return eulerAngle(atan(-2*P*q0*q3, q0^2 - q3^2), 0, 0)
    elseif (x == 0 && q03 == 0)
        return eulerAngle(atan(2*q1*q2, q1^2 - q2^2), pi, 0)
    else #als x != 0
        return eulerAngle(atan((q1*q3 - P*q0*q2)/x,(-P*q0*q1 - q2*q3)/x),
            atan(2*x, q03-q12), atan((P*q0*q2 + q1*q3)/x, (q2*q3 - P*q0*q1)/x))
    end
end

function quaternionToRotationMatrix(rotation ::quaternion)
    #moet unit quaternion zijn
    #geeft passieve interpretatie
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    q = q0^2 - (q1^2 + q2^2 + q3^2)
    a = [q+2*q1^2 2*(q1*q2-P*q0*q3) 2*(q1*q3+P*q0*q2);
        2(q1*q2+P*q0*q3) q+2*q2^2 2*(q2*q3-P*q0*q1);
        2*(q1*q3-P*q0*q2) 2*(q2*q3+P*q0*q1) q+2*q3^2]
    return rotationMatrix(a)
end

function quaternionToAxisAngle(rotation ::quaternion)
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    w = 2*acos(q0)
    if (w == 0)
        return axisAnglePair([0, 0, 1], 0)
    end
    if (q0 == 0)
        return axisAnglePair([q1, q2, q3], pi)
    end
    s = sign(q0)/sqrt(q1^2 + q2^2 + q3^2)
    return axisAnglePair([s*q1, s*q2, s*q3], w)
end

function quaternionToRodriguesFrank(rotation ::quaternion)
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    s = sqrt(q1^2 + q2^2 + q3^2)
    t = tan(acos(q0))
    #oppassen als s klein wordt!!!
    #rodriguesFrank opslaan als vector van 4 elementen
    return rodriguesFrank([q1/s, q2/s, q3/s], t)
end

function quaternionToHomochoric(rotation ::quaternion)
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    w = 2*acos(q0)
    if w == 0
        return homochoric([0, 0, 0])
    end
    s = 1/sqrt(q1^2 + q2^2 + q3^2)
    n = [s*q1, s*q2, s*q3]
    f = 3(w-sin(w))/4
    return homochoric(n*f^(1/3))
end

function homochoricToAxisAngle()
end

function homochoricToCubochoric()
end

function qubochoricToHomochoric()
end
