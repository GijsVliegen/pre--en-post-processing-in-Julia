struct quaternion
    i
    j
    k
    angle
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