struct quaternion
    i
    j
    k
    angle
end
function isClose(first, other, rtol, atol, equal_nan = True)
    if isnan(first)
        if isnan(other)
            return equal_nan
        return false
        end
    end    
    fQu = toQuaternion(first)
    sQu = toQuaternion(other)
    return absolute(fQu - sQu) <= (atol + rtol * absolute(b))
end    