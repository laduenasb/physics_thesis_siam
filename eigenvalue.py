#########################################################################
############### Autovalores del hamiltoniano de Anderson ################
#########################################################################
from __future__ import division
from sympy import *
pi=pi.evalf()
#########################################################################
########################## Parametros ###################################
#########################################################################
delta_0=0.01
E_f=-10.0*delta_0 # Energia del estado localizado
E_q=0.0 # Energia de la banda de conduccion
U=20.0*delta_0 # Repulsion columbiana
D=1.0
V=sqrt(delta_0*2.0*D/pi) # Hibridacion

x1=E_f # Energia del estado localizado
x2=E_q # Energia de la banda de conduccion
#########################################################################
#########################################################################
U_0=U
V_0=V
a_2=-(3*x1+3*x2+U)
delta=((x1-x2)**2+4*V**2)**0.5
delta_p=((x1+U-x2)**2+4*V**2)**0.5
Q=-(12*V**2+(x1+x2)**2+(2*x1+U)**2+(2*x2)**2-(x2+x1)*(2*x1+U)-(x2+x1)*(2*x2)-(2*x1+U)*(2*x2))/9.0
R=(-3*((x2+x1)**2*(2*x1+U)+(x2+x1)**2*(2*x2)+(2*x1+U)**2*(x2+x1)+(2*x1+U)**2*(2*x2)+(2*x2)**2*(x2+x1)+(2*x2)**2*(2*x1+U))+12*(x2+x1)*(2*x1+U)*(2*x2)+18*V**2*(2*(x2+x1)-(2*x1+U)-2*x2)+2*((x2+x1)**3+(2*x1+U)**3+(2*x2)**3))/54
theta_1=acos(R/sqrt(-Q**3))
E_1=0.0
E_2=0.5*(x1+x2-delta)
E_3=E_2
E_4=0.5*(x1+x2+delta)
E_5=E_4
E_6=x1+x2
E_7=E_6
E_8=E_7
E_9=2*sqrt(-Q)*cos(theta_1/3.0)-a_2/3.0
E_10=2*sqrt(-Q)*cos((theta_1+2*pi)/3.0)-a_2/3.0
E_11=2*sqrt(-Q)*cos((theta_1+4*pi)/3.0)-a_2/3.0
E_12=0.5*(3*x1+3*x2+U+delta_p)
E_13=E_12
E_14=0.5*(3*x1+3*x2+U-delta_p)
E_15=E_14
E_16=2*x1+2*x2+U
E=[E_1,E_2,E_3,E_4,E_5,E_6,E_7,E_8,E_9,E_10,E_11,E_12,E_13,E_14,E_15,E_16]
##################################    
########## Matrix H ##############
##################################
#E_f,E_q,V,U=symbols("epsilon_f epsilon_q V U")
E_f=x1
E_q=x2
U=U_0
V=V_0
def f(i,j):
    return 0
H=Matrix(16,16,f)
H[0,0]=0.0
H[1,1]=E_f
H[2,2]=E_q
H[3,3]=E_f
H[4,4]=E_q
H[5,5]=E_f+E_q
H[6,6]=E_f+E_q
H[7,7]=E_f+E_q
H[8,8]=E_f+E_q
H[9,9]=2*E_q
H[10,10]=2*E_f+U
H[11,11]=2*E_f+E_q+U
H[12,12]=E_f+2*E_q
H[13,13]=2*E_f+E_q+U
H[14,14]=E_f+2*E_q
H[15,15]=2*E_f+2*E_q+U
H[1,2]=V
H[2,1]=H[1,2]
H[3,4]=V
H[4,3]=H[3,4]
H[7,9]=V
H[9,7]=V
H[7,10]=V
H[10,7]=H[7,10]
H[8,9]=-V
H[9,8]=H[8,9]
H[8,10]=-V
H[10,8]=-V
H[11,12]=V
H[12,11]=H[11,12]
H[13,14]=V
H[14,13]=H[13,14]
print "Matriz H\n"
pprint(H)
print ""
print "Eigenvalues\n"
for i in range(0,len(E)):
    print "E"+str(i+1)+":",E[i]
E_S=H.eigenvals()
E_Sn=E_S.items()
E_n=[]
for i in range(0,len(E_Sn)):
    E_n.append((E_Sn[i][0]).evalf())
pprint(E_S)
print "\n \n"
pprint(E_n)
pprint(H.eigenvects())
