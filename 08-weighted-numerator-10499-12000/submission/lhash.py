import struct
M=(1<<64)-1
def murmur(data, seed=11):
    m=0xc6a4a7935bd1e995; r=47
    n=len(data)
    h=(seed ^ ((n*m)&M))&M
    nb=n//8
    for i in range(nb):
        k=struct.unpack_from('<Q',data,i*8)[0]
        k=(k*m)&M; k^=k>>r; k=(k*m)&M
        h^=k; h=(h*m)&M
    tail=data[nb*8:]
    rem=n&7
    if rem:
        for j in range(rem-1,-1,-1):
            h^=tail[j]<<(8*j)
        h=(h*m)&M
    h^=h>>r; h=(h*m)&M; h^=h>>r
    return h
def mix(h,k):
    m=0xc6a4a7935bd1e995
    k=(k*m)&M; k^=k>>47; k^=m
    h^=k; h=(h*m)&M
    return h
def ofString(b): return mix(1723, murmur(b))
def ofBytes(b): return mix(1723, murmur(b))
