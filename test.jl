
quOne = quaternion(1, 1, 1, 0)
quTwo = quaternion(1+1e-1, 1, 1, 0)
if isClose(quOne, quTwo, 0.01, 0.01)
    println("yeet")
else
    println("niet yeet")
end    