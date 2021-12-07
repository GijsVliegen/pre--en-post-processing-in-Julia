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

function isClose(first ::quaternion, other ::quaternion, rtol, atol)
    fcom = getComponents(toQuaternion(first))
    scom = getComponents(toQuaternion(other))
    for i in 1:4
        if !isapprox(fcom[i], scom[i], rtol = rtol , atol = atol)
            return false
        end
    end
    return true
end    
quOne = quaternion(1, 1, 1, 0)
quTwo = quaternion(1+1e-1, 1, 1, 0)
if isClose(quOne, quTwo, 0.01, 0.01)
    println("yeet")
else
    println("niet yeet")
end    