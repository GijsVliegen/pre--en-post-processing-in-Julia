import damask, time
startTime = time.time()
processST = time.process_time()

a = damask.Rotation.from_random(10000)
b = a.as_matrix()
c = damask.Rotation.from_matrix(b)

# get the end time
endTime = time.time()
processET = time.process_time()

# get the execution time
elapsed_time = endTime - startTime
elapsed_process_time = processET - processST
print("elapsed time: ", elapsed_time * 10**3, "ms")
print('Execution time:', elapsed_process_time * 10**3, 'ms')