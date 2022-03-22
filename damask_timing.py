from functools import total_ordering
import damask, time


#time
def operationEu(a, l):
    totaltime = 0
    for i in range(l):
        processST = time.process_time()
        b = a.as_Euler_angles()
        c = damask.Rotation.from_Euler_angles(b)
        processET = time.process_time()
        totaltime += processET - processST
    print(totaltime/l * 10**3)

def operationOm(a, l):
    totaltime = 0
    for i in range(l):
        processST = time.process_time()
        b = a.as_matrix()
        c = damask.Rotation.from_matrix(b)
        processET = time.process_time()
        totaltime += processET - processST
    print(totaltime/l * 10**3)

def operationAx(a, l):
    totaltime = 0
    for i in range(l):
        processST = time.process_time()
        b = a.as_axis_angle()
        c = damask.Rotation.from_axis_angle(b)
        processET = time.process_time()
        totaltime += processET - processST
    print(totaltime/l * 10**3)

def operationHo(a, l):
    totaltime = 0
    for i in range(l):
        processST = time.process_time()
        b = a.as_homochoric()
        c = damask.Rotation.from_homochoric(b)
        processET = time.process_time()
        totaltime += processET - processST
    print(totaltime/l * 10**3)

def operationRo(a, l):
    totaltime = 0
    for i in range(l):
        processST = time.process_time()
        b = a.as_Rodrigues_vector()
        c = damask.Rotation.from_Rodrigues_vector(b)
        processET = time.process_time()
        totaltime += processET - processST
    print(totaltime/l * 10**3)

def oneTest(a, l):
    operationEu(a, l)
    operationOm(a, l)
    operationAx(a, l)
    operationHo(a, l)
    operationRo(a, l)

def longTestRun():
    print("EU - OM - AX - HO - RO")
    for i in range (8,9):
        l = 100
        n = 2**i
        print(n)
        a = damask.Rotation.from_random(n)

oneTest(damask.Rotation.from_random(10000), 100)
oneTest(damask.Rotation.from_random((10, 10, 10, 10)), 100)
