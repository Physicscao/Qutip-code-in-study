from qutip import *
import numpy as np

#Grover 
n=int(input('数据个数：'))
m=int(input('目标态的位置：'))
#制作目标态
psi=fock(n,m-1); 
#制作完备基态
h=1
Psi=basis(n)
for h in range(1,n):
    Psi=Psi+fock(n,h)
#完备基态归一化
Psi_in=Psi.unit()
#设定Grover算符
U_theta=qeye(n)-2*(psi*psi.dag()).unit();
U_psi=qeye(n)-2*(Psi*Psi.dag()).unit();
U=-U_psi*U_theta;
#计算Grover迭代次数
N=np.pi*np.sqrt(n)/4;
if N-int(N)<0.5:
    flag=int(N)
else:
    flag=int(N)+1   
#进行Grover迭代
i=1;
while i<= flag:
    Psi_in=U*Psi_in
    i=i+1
#输出搜索概率
a=(Psi_in*Psi_in.dag()).unit().diag()
print(a[m-1])