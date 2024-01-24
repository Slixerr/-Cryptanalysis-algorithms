import math
import os
import random
import pandas as pd
import time
import psutil

memory = 0

def fermat(n):
    start_timer = time.time()
    a = math.ceil(math.sqrt(n))
    b = a**2 - n
    while not b == math.isqrt(b) ** 2:
        a += 1
        b = a**2 - n
        if(time.time() - start_timer > 1200): return None
    
    global memory
    memory = round(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 2)
    
    return a - math.sqrt(b), a + math.sqrt(b)
# def fermat(n, verbose=False):
#     start_timer = time.time()
#     a = math.isqrt(n) # int(ceil(n**0.5))
#     b2 = a*a - n
#     b = math.isqrt(n) # int(b2**0.5)
#     count = 0
#     while b*b != b2:
#         if(time.time() - start_timer > 1200): return None
#         if verbose:
#             print('Trying: a=%s b2=%s b=%s' % (a, b2, b))
#         a = a + 1
#         b2 = a*a - n
#         b = math.isqrt(b2) # int(b2**0.5)
#         count += 1
#     p=a+b
#     q=a-b
#     global memory
#     memory = round(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 2)
#     assert n == p * q
#     return p, q
#print(fermat(44461))

######################################################################

def pollard_rho_factorizacion(n):
    start_timer = time.time()

    randomNum = random.randint(2, n-1)
    a = randomNum
    b = randomNum
    while True:
        a = (a**2 + 1) % n
        b = (b**2 + 1) % n
        b = (b**2 + 1) % n
        p = math.gcd(a-b, n)

        if(time.time() - start_timer > 1200): return None
        global memory
        memory = round(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 2)
        if(p > 1 and p < n): 
            return p
        if(p == n): 
            return n

#print(pollard_rho(39617))

######################################################################

def pollard_P_1(n):
    start_timer = time.time()
    a = random.randint(2, n-1)
    if(math.gcd(a, n) > 1 and math.gcd(a, n) < n): return math.gcd(a, n)
    k = 2
    while True:
        a = (a**k) % n
        d = math.gcd(a-1, n)
        if(d > 1 and d < n): return d
        if(d == n): return False
        k += 1

        global memory
        memory = round(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 2)
        if(time.time() - start_timer > 1200): return None





######################################################################Logaritmo Discreto############################################################

def babyStepGiantStep(p, alpha, beta):
    start_timer = time.time()
    n = math.ceil(math.sqrt(p))
    t = {}
    for r in range(n):
        t[pow(alpha, r, p)] = r
    
    alpha_invN = pow(alpha, -n, p)
    gamma = beta
    for q in range(n):
        global memory
        memory = round(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 2)
        if gamma in t:
            j = t[gamma]
            k = q * n + j
            return k
        gamma = (gamma * alpha_invN) % p
        if(time.time() - start_timer > 1200): return None
    return False

######################################################################
def f(x, a, b, alpha, beta, p):
    x = ((alpha**a) * (beta**b)) % p
    modulo = x % 3
    if modulo == 1:
        x_ii = (beta * x) % p
        a_ii = a
        b_ii = (b+1) % (p-1)
    elif modulo == 0:
        x_ii = (x**2) % p
        a_ii = (2*a) % (p-1)
        b_ii = (2*b) % (p-1)
    else: #modulo == 2
        x_ii = (alpha * x) % p
        a_ii = (a+1) % (p-1)
        b_ii = b
    
    return x_ii, a_ii, b_ii


def pollard_rho_logDisc(p, alpha, beta, o):
    start_timer = time.time()
    a = b = aa = bb = 0
    i = x = xx = 1
    
    while i < p:
        if(time.time() - start_timer > 2400): return None
        x,a,b = f(x, a, b, alpha, beta, p)
        xx,aa,bb = f(xx, aa, bb, alpha, beta, p)
        xx,aa,bb = f(xx, aa, bb, alpha, beta, p)
        
        if x == xx:
            if math.gcd(b - bb, o) != 1:
                return False
            return ((aa - a) * pow(b - bb, -1, o)) % o
        i += 1
        global memory
        memory = round(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 2)
        
    return False


# p = 853  
# alpha = 9
# beta = 804 
# o = 71

# print(pollard_rho_logDisc(p, alpha, beta, o))

######################################################################

def inverso(a, b):
    if b == 0:
        return 1, 0, a
    q, r = divmod(a, b)
    x, y, g = inverso(b, r)
    return y, x - q * y, g

def sumaEnCurvasElipticas(p, q, a, b, m):
    if p[2] == 0: return q
    if q[2] == 0: return p
    if p[0] == q[0]:
        if (p[1] + q[1]) % m == 0:
            return 0, 1, 0 
        numerador = (3 * pow(p[0],2) + a) % m
        denominador = (2 * p[1]) % m
    else: 
        numerador = (q[1] - p[1]) % m
        denominador = (q[0] - p[0]) % m
    inv, _, gcd = inverso(denominador, m)

    if gcd > 1:
        return 0, 0, denominador  
    
    pendiente = numerador * inv
    xr = (pendiente**2 - p[0] - q[0]) % m 
    yr = (pendiente * (p[0] - xr) - p[1]) % m 
    return xr, yr, 1 



def exponCuadradosSuc(k, p, a, b, m):
    result = (0, 1, 0)  
    while k > 0:
        if p[2] > 1:
            return p
        if k % 2 == 1: 
            result = sumaEnCurvasElipticas(p, result, a, b, m)
        k = k // 2
        p = sumaEnCurvasElipticas(p, p, a, b, m)

    return result 


def lenstra(n, cota):
    timer = 0
    start_timer = time.time()
    
    while timer < 1200:
        gcd = n
        while gcd == n:
            p = random.randint(0, n - 1), random.randint(0, n - 1),1
            a = random.randint(0, n - 1)
            b = (pow(p[1],2) - pow(p[0],3) - a * p[0]) % n
            gcd = math.gcd(4 * pow(a,3) + 27 * pow(b,2), n)

        global memory
        memory = round(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 2)
        if gcd > 1:
            return gcd

        k = 2
        while k < cota:
            p = exponCuadradosSuc(k, p, a, b, n)
            if p[2] > 1:
                return math.gcd(p[2], n)
            k += 1
        timer = time.time() - start_timer
        


if __name__ == '__main__':
    dataFact = pd.read_csv("datosFactorizacionExt.csv")
    dataFact.set_index("Tamaño", inplace=True)

    resultAnalysisCSV = pd.DataFrame()
    for index, row in dataFact.iterrows():
        start_time = time.time()
        result = fermat(int(row['n']))
        end_time = time.time() - start_time

        data = {
            'Tamaño': index,
            'Número': row['n'],
            'Solución': [str(result)],
            'Tiempo': round(end_time, 3),
            'Memoria': memory
        }
        new_data = pd.DataFrame(data)
        resultAnalysisCSV = pd.concat([resultAnalysisCSV, pd.DataFrame(new_data)], ignore_index=True)
        print(resultAnalysisCSV)
        resultAnalysisCSV.to_csv("Resultados/Fermat.csv", index=False)

    
    resultAnalysisCSV = pd.DataFrame()
    for index, row in dataFact.iterrows():
        start_time = time.time()
        result = pollard_rho_factorizacion(int(row['n']))
        end_time = time.time() - start_time

        data = {
            'Tamaño': index,
            'Número': row['n'],
            'Solución': [str(result)],
            'Tiempo': round(end_time, 3),
            'Memoria': memory
        }
        new_data = pd.DataFrame(data)
        resultAnalysisCSV = pd.concat([resultAnalysisCSV, pd.DataFrame(new_data)], ignore_index=True)
        
        print(resultAnalysisCSV)
        resultAnalysisCSV.to_csv("Resultados/Pollard_rho_factorizacion.csv", index=False)
    
    resultAnalysisCSV = pd.DataFrame()
    for index, row in dataFact.iterrows():
        start_time = time.time()
        result = lenstra(int(row['n']), 100)
        end_time = time.time() - start_time

        data = {
            'Tamaño': index,
            'Número': row['n'],
            'Solución': [str(result)],
            'Tiempo': round(end_time, 3),
            'Memoria': memory
        }
        new_data = pd.DataFrame(data)
        resultAnalysisCSV = pd.concat([resultAnalysisCSV, pd.DataFrame(new_data)], ignore_index=True)
        print(resultAnalysisCSV)
        resultAnalysisCSV.to_csv("Resultados/Lenstra.csv", index=False)

    resultAnalysisCSV = pd.DataFrame()
    for index, row in dataFact.iterrows():
        start_time = time.time()
        result = pollard_P_1(int(row['n']))
        end_time = time.time() - start_time

        data = {
            'Tamaño': index,
            'Número': row['n'],
            'Solución': [str(result)],
            'Tiempo': round(end_time, 3),
            'Memoria': memory
        }
        new_data = pd.DataFrame(data)
        resultAnalysisCSV = pd.concat([resultAnalysisCSV, pd.DataFrame(new_data)], ignore_index=True)
        print(resultAnalysisCSV)
        resultAnalysisCSV.to_csv("Resultados/Pollard_P_1.csv", index=False)

    
    

######################################################Logaritmo Discreto###################################################
    
    dataFact = pd.read_csv("datosLogaritmoExt.csv")
    dataFact.set_index("Tamaño", inplace=True)

    resultAnalysisCSV = pd.DataFrame()
    for index, row in dataFact.iterrows():
        start_time = time.time()
        result = babyStepGiantStep(int(row['m']), int(row['alfa']), int(row['beta']))
        end_time = time.time() - start_time

        data = {
            'Tamaño': index,
            'Número': row['m'],
            'Solución': [str(result)],
            'Tiempo': round(end_time, 3),
            'Memoria': memory
        }
        new_data = pd.DataFrame(data)
        resultAnalysisCSV = pd.concat([resultAnalysisCSV, pd.DataFrame(new_data)], ignore_index=True)

        print(resultAnalysisCSV)
        resultAnalysisCSV.to_csv("Resultados/babyStepGiantStepExt.csv", index=False)

    
    resultAnalysisCSV = pd.DataFrame()

    for index, row in dataFact.iterrows():
        start_time = time.time()
        result = pollard_rho_logDisc(int(row['m']), int(row['alfa']), int(row['beta']), int(row['orden']))
        end_time = time.time() - start_time

        data = {
            'Tamaño': index,
            'Número': row['m'],
            'Solución': [str(result)],
            'Tiempo': round(end_time, 3),
            'Memoria': memory
        }
        new_data = pd.DataFrame(data)
        resultAnalysisCSV = pd.concat([resultAnalysisCSV, pd.DataFrame(new_data)], ignore_index=True)
        
        print(resultAnalysisCSV)
        resultAnalysisCSV.to_csv("Resultados/pollard_rho_logDiscExt.csv", index=False)

    
    resultAnalysisCSV = pd.DataFrame()
    

