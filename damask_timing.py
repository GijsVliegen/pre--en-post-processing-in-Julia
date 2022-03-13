import damask, time
startTime = time.time()
processST = time.process_time()

#timeit

a = damask.Rotation.from_random(50000)
b = a.as_Euler_angles()
c = damask.Rotation.from_Euler_angles(b)

# get the end time
endTime = time.time()
processET = time.process_time()

# get the execution time
elapsed_time = endTime - startTime
elapsed_process_time = processET - processST
print("elapsed time: ", elapsed_time * 10**3, "ms")
print('Execution time:', elapsed_process_time * 10**3, 'ms')