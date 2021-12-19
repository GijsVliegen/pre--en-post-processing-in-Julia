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

function getComponents(rot ::quaternion)
    return [rot.i, rot.j, rot.k, rot.angle]
end

function toQuaternion(rot ::quaternion)
    return rot
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

function rotationMatrixToEulerAngle()
end

function rotationMatrixToAxisAngle()
end

function rotationMatrixToQuaternion()
end

function axisAngleToRotationMatrix()
end

function axisAngleToRodriguesFrank()
end

function axisAngleToHomochoric()
end

function rodriguesFrankToAxisAngle()
end

function rodriguesFrankToHomochoric()
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