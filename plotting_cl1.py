import matplotlib.pyplot as plt
import math
import scipy as sp
import numpy as np

data = np.loadtxt('efgh_k1C0p1.txt');
#dataC = np.loadtxt('rsltsZF.txt');

N = 30;
n1 = 3600;
Nav = 400; 
k = data[0:N,0];
E = data[n1*N:(n1+1)*N,1]; 
F = data[n1*N:(n1+1)*N,2]; 
G = data[n1*N:(n1+1)*N,3]; 
H = data[n1*N:(n1+1)*N,4]; 

for n in range(Nav): 
    E = E + data[(n1+n)*N:(n1+n+1)*N,1];
    F = F + data[(n1+n)*N:(n1+n+1)*N,2]; 
    G = G + abs(data[(n1+n)*N:(n1+n+1)*N,3]);
    H = H + abs(data[(n1+n)*N:(n1+n+1)*N,4]);
   




plt.loglog(k,F,'b')
plt.loglog(k,k**(-5./3.)*math.exp(7),'g')
plt.loglog(k,k**(-1.)*math.exp(5.5),'r')
plt.title('F(k) for $C\ll1$ with lines $k^{-5/3}$ in green and $k^{-1}$ in red')
plt.show() 

 
plt.loglog(k,E,'b') 
plt.loglog(k,k**(-5./3.)*math.exp(7),'g')
plt.loglog(k,k**(-3.)*math.exp(10),'r')
plt.title('E  for $C\ll1$ with lines $k^{-5/3}$ in green and $k^{-3}$ in red')
plt.show()
 
 
plt.loglog(k,G/k,'b')
#plt.loglog(k,k**(-5./3.)*math.exp(10),'r')
plt.loglog(k,k**(-7./3.)*math.exp(5.5),'g')
#plt.loglog(k,k**(-2)*math.exp(15),'r')
plt.loglog(k,k**(-3.)*math.exp(7),'r')
plt.title('G/k  for $C\ll1$ with lines $k^{-7/3}$ in green and $k^{-3}$ in red')

plt.show()



plt.loglog(k,G,'b')
plt.loglog(k,np.sqrt(E*F),'g')
plt.title('G  (in blue) and $\sqrt{EF}$ (in green) for $C\ll 1$') 
plt.show()
