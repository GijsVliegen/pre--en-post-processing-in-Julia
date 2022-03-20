from functools import total_ordering
import damask, time

a = damask.Rotation.from_random(200000)

#time
def operationEu(a):
    totaltime = 0
    for i in range(100):
        processST = time.process_time()
        b = a.as_Euler_angles()
        c = damask.Rotation.from_Euler_angles(b)
        processET = time.process_time()
        totaltime += processET - processST
    print('Execution time Eu:', totaltime/100 * 10**3, 'ms')

def operationOm(a):
    totaltime = 0
    for i in range(100):
        processST = time.process_time()
        b = a.as_matrix()
        c = damask.Rotation.from_matrix(b)
        processET = time.process_time()
        totaltime += processET - processST
    print('Execution time Om:', totaltime/100 * 10**3, 'ms')

def operationAx(a):
    totaltime = 0
    for i in range(100):
        processST = time.process_time()
        b = a.as_axis_angle()
        c = damask.Rotation.from_axis_angle(b)
        processET = time.process_time()
        totaltime += processET - processST
    print('Execution time Ax:', totaltime/100 * 10**3, 'ms')

def operationHo(a):
    totaltime = 0
    for i in range(100):
        processST = time.process_time()
        b = a.as_homochoric()
        c = damask.Rotation.from_homochoric(b)
        processET = time.process_time()
        totaltime += processET - processST
    print('Execution time Ho:', totaltime/100 * 10**3, 'ms')

def operationRo(a):
    totaltime = 0
    for i in range(100):
        processST = time.process_time()
        b = a.as_Rodrigues_vector()
        c = damask.Rotation.from_Rodrigues_vector(b)
        processET = time.process_time()
        totaltime += processET - processST
    print('Execution time Ro:', totaltime/100 * 10**3, 'ms')

operationEu(a)
operationOm(a)
operationAx(a)
operationHo(a)
operationRo(a)