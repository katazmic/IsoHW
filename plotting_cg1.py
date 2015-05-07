import matplotlib.pyplot as plt
import math
import scipy as sp
import numpy as np

data = np.loadtxt('efgh_k1C40.txt');
#dataC = np.loadtxt('rsltsZF.txt');

N = 30;
n1 = 900;
Nav = 150; 
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
plt.loglog(k[5:11],k[5:11]**(-11./3.)*math.exp(4.5),'g')
plt.loglog(k[11:17],k[11:17]**(-5)*math.exp(6.5),'r')
plt.loglog(k[17:25],k[17:25]**(-3.5)*math.exp(1.2),'y')
plt.title('F(k) for $C\gg1$ with lines $k^{-11/3}$ in green and $k^{-5}$ in red $k^{-3.5}$ in yellow' )
plt.show() 
  
plt.loglog(k,E,'b')
plt.loglog(k,k**(-5./3.)*math.exp(4),'g')
plt.loglog(k,k**(-3.5)*math.exp(8),'r')
plt.title('E  for $C\gg1$ with lines $k^{-5/3}$ in green and $k^{-3.5}$ in red')
plt.show()



plt.loglog(k,G,'b')
plt.loglog(k,k**(-5./3.)*math.exp(0.52),'g')
plt.loglog(k,k**(-3.5)*math.exp(4),'r')
plt.title('G  for $C\gg 1$ with lines $k^{-5/3}$ in green and $k^{-3.5}$ in red')

 
plt.show()


plt.loglog(k,G,'b')
plt.loglog(k,np.sqrt(E*F),'g')
plt.title('G  (in blue) and $\sqrt{EF}$ (in green) for $C\gg 1$') 
plt.show()
