###################################
########## Matriz H ###############
###################################
from sympy import *
E_f,E_q,V,U=symbols("epsilon_f epsilon_q V U")
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
