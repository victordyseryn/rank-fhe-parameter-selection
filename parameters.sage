q = 2

def key_size(m,n,w):
    size = (m*m+n*w)//8
    return size

def ct_size(m,n,w):
    size = (2*m*n)//8
    return size

def add_time(m,n):
    time = (2*m*n)/3000000
    if time < 1:
        return N(time)
    else:
        return floor(time)

def mul_time(m,n):
    time = (3*(m*n)**1.6)/3000000
    if time < 1:
        return N(time)
    else:
        return floor(time)

def boot_time(m,n):
    time = (2*(m**2)*(n**1.6))/3000000
    if time < 1:
        return N(time)
    else:
        return floor(time)

def rGV_bound(q, n, k, m):
    w = 1
    while (q_binomial(m,w,q)*q**(w*n) <= q**(m*(n-k))):
        w += 1
    return w-1

def cpx_comb_RSD(q,n,k,m,w):
    comb_RSD = (((n-k)*m)**3)*(q**(w*ceil((k+1)*m/n)-m))
    return int(log(max(1,comb_RSD),2))

def cpx_alg_RSD(q,n,k,m,w):
    a = 0
    while m*binomial(n-k-1,w) < binomial(n-a,w)-1:
        a += 1
    alg_RSD = q**(a*w)*m*binomial(n-k-1,w)*(binomial(n-a,w)**2)
    return int(log(max(1,alg_RSD),2))

def parameter_selection(d, start_m):
    for m in range(start_m, start_m+1000, 1):
        #print("m = ", m)
        n = 10
        while n < m/2:
            if d == 1:
                w = rGV_bound(q, 2*n, n, m)
            else:
                w = int(m**(1/d))
            #print("w = ", w)
            for l in range(2,int(3*w/4)+2):
                cpx = cpx_comb_RSD(q,n*l,n,m,w)
                if cpx < 143:
                    #print("n = ", n, ",l = ",l)
                    break
                cpx2 = cpx_alg_RSD(q,n*l,n,m,w)
                if cpx2 < 143:
                    #print("n = ", n, ",l = ",l)
                    break

            if l == int(3*w/4)+1:
                found = True
                print("-------------")
                print("m = ", m)
                print("n = ", n)
                print("w = ", w)
                print("l = ", l-1)
                print("-------------")
                print("key = ", key_size(m,n,w), "B")
                print("ct  = ", ct_size(m,n,w), "B")
                print("add = ", add_time(m,n), "ms")
                print("mul = ", mul_time(m,n), "ms")
                print("bts = ", boot_time(m,n), "ms")
                return
            else:
                n += 1

print("========== d = 0 ==========")
parameter_selection(1, 100)
print("")

print("========== d = 1 ==========")
parameter_selection(2, 100)
print("")

print("========== d = 2 ==========")
parameter_selection(3, 250)
print("")

print("========== d = 3 ==========")
parameter_selection(4, 1000)
print("")

print("========== d = 4 ==========")
parameter_selection(5, 2500)
print("")






