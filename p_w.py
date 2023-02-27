#Calculo de la densidad de estados
from __future__ import division
from math import *
from numpy import *
import numpy as np
#from sympy import atan,acos,sin,cos,sqrt,exp,log,lambdify,Symbol
import time
#x=Symbol("x") # w alternativo
delta_0=0.01
E_f=-10.0*delta_0 # Energia del estado localizado
E_q=-1.168e-6*delta_0
U=(20.0)*delta_0 # Repulsion columbiana
mu=0.0
D=1.0
A=-D # intervalo [A,B]
B=D # intervalo [A,B]
A_n=-1.0
B_n=1.0
V=sqrt(delta_0*2.0*D/pi) # Hibridacion
T=0.001*delta_0
beta=1.0/T
n_puntos=1001
eta_=1e-12
##############################
E_f=np.float64(E_f)
E_q=np.float64(E_q)
#E_q_p=np.float64(E_q_p)
U=np.float64(U)
mu=np.float64(mu)
D=np.float64(D)
A_n=np.float64(A_n)
B_n=np.float64(B_n)
eta_=np.float64(eta_)
#########################################################################
##############################Autovalores################################
#########################################################################
fmin=np.float64(1.0)
Eqmin=np.float64(0.0)
pho_fsmumin=np.float64(0.0)
n_f_pp_1min=np.float64(1.0)
def a(E_i,E_f,E_q,U,V):
    f=1.0/sqrt(2.0+4*V**2*(1.0/(E_i-2*E_f-U)**2+1.0/(E_i-2*E_q)**2))
    return f
def b(E_i,E_f,E_q,U,V):
    f=2.0*V*a(E_i,E_f,E_q,U,V)/(E_i-2*E_f-U)
    return f
def c(E_i,E_f,E_q,U,V,E_9):
    f=2.0*V*a(E_i,E_f,E_q,U,V)/(E_9-2*E_q)
    return f
def gkondo_s(E_q,V):
    x1=E_f
    x2=E_q
    a_2_p=-(3*x1+3*x2+U)
    delta=sqrt((x1-x2)**2+4*V**2)
    delta_p=sqrt((x1+U-x2)**2+4*V**2)
    Q=-(12*V**2+(x1+x2)**2+(2*x1+U)**2+(2*x2)**2-(x2+x1)*(2*x1+U)-(x2+x1)*(2*x2)-(2*x1+U)*(2*x2))/9.0
    R=(-3*((x2+x1)**2*(2*x1+U)+(x2+x1)**2*(2*x2)+(2*x1+U)**2*(x2+x1)+(2*x1+U)**2*(2*x2)+(2*x2)**2*(x2+x1)+(2*x2)**2*(2*x1+U))+12*(x2+x1)*(2*x1+U)*(2*x2)+18*V**2*(2*(x2+x1)-(2*x1+U)-2*x2)+2*((x2+x1)**3+(2*x1+U)**3+(2*x2)**3))/54.0
    theta_1=acos(R/sqrt(-Q**3))
    E_1=np.float64(0.0)
    E_2=0.5*(x1+x2-delta)
    E_3=E_2
    E_4=0.5*(x1+x2+delta)
    E_5=E_4
    E_6=x1+x2
    E_7=E_6
    E_8=E_7
    E_9=2*sqrt(-Q)*cos(theta_1/3.0)-a_2_p/3.0
    E_10=2*sqrt(-Q)*cos((theta_1+2*pi)/3.0)-a_2_p/3.0
    E_11=2*sqrt(-Q)*cos((theta_1+4*pi)/3.0)-a_2_p/3.0
    E_12=0.5*(3*x1+3*x2+U+delta_p)
    E_13=E_12
    E_14=0.5*(3*x1+3*x2+U-delta_p)
    E_15=E_14
    E_16=2*x1+2*x2+U
    E=[E_1,E_2,E_3,E_4,E_5,E_6,E_7,E_8,E_9,E_10,E_11,E_12,E_13,E_14,E_15,E_16]
    E_min=min(E)
    E_1=E_1-E_min
    E_2=E_2-E_min
    E_3=E_3-E_min
    E_4=E_4-E_min
    E_5=E_5-E_min
    E_6=E_6-E_min
    E_7=E_7-E_min
    E_8=E_8-E_min
    E_9=E_9-E_min
    E_10=E_10-E_min
    E_11=E_11-E_min
    E_12=E_12-E_min
    E_13=E_13-E_min
    E_14=E_14-E_min
    E_15=E_15-E_min
    E_16=E_16-E_min
    exp_E=[exp(-beta*E_1),exp(-beta*E_2),exp(-beta*E_3),exp(-beta*E_4),exp(-beta*E_5),exp(-beta*E_6),exp(-beta*E_7),exp(-beta*E_8),exp(-beta*E_9),exp(-beta*E_10),exp(-beta*E_11),exp(-beta*E_12),exp(-beta*E_13),exp(-beta*E_14),exp(-beta*E_15),exp(-beta*E_16)]
    Z_=0.0 #funcion de particion canonica
    for i in range(0,len(E)):
        Z_+=exp_E[i]
    SEXP_1=[exp_E[0]+exp_E[3],exp_E[1]+exp_E[5],exp_E[0]+exp_E[1],exp_E[3]+exp_E[5],exp_E[2]+exp_E[8],exp_E[2]+exp_E[9],exp_E[2]+exp_E[10],exp_E[4]+exp_E[8],exp_E[4]+exp_E[9],exp_E[4]+exp_E[10],exp_E[8]+exp_E[13],exp_E[9]+exp_E[13],exp_E[10]+exp_E[13],exp_E[8]+exp_E[11],exp_E[9]+exp_E[11],exp_E[10]+exp_E[11]]
    SEXP_2=[exp_E[14]+exp_E[15],exp_E[7]+exp_E[11],exp_E[12]+exp_E[15],exp_E[7]+exp_E[13],exp_E[3]+exp_E[8],exp_E[3]+exp_E[9],exp_E[3]+exp_E[10],exp_E[2]+exp_E[8],exp_E[2]+exp_E[9],exp_E[2]+exp_E[10],exp_E[8]+exp_E[11],exp_E[9]+exp_E[11],exp_E[10]+exp_E[11],exp_E[8]+exp_E[13],exp_E[9]+exp_E[13],exp_E[10]+exp_E[13]]
    SEXP_3=[exp_E[4]+exp_E[8],exp_E[4]+exp_E[9],exp_E[4]+exp_E[10],exp_E[2]+exp_E[8],exp_E[2]+exp_E[9],exp_E[2]+exp_E[10],exp_E[8]+exp_E[13],exp_E[9]+exp_E[13],exp_E[10]+exp_E[13],exp_E[8]+exp_E[11],exp_E[9]+exp_E[11],exp_E[10]+exp_E[11]]
    SEXP_4=[exp_E[0]+exp_E[4],exp_E[2]+exp_E[6],exp_E[0]+exp_E[2],exp_E[4]+exp_E[6],exp_E[1]+exp_E[8],exp_E[1]+exp_E[9],exp_E[1]+exp_E[10],exp_E[3]+exp_E[8],exp_E[3]+exp_E[9],exp_E[3]+exp_E[10],exp_E[8]+exp_E[14],exp_E[9]+exp_E[14],exp_E[10]+exp_E[14],exp_E[8]+exp_E[12],exp_E[9]+exp_E[12],exp_E[10]+exp_E[12]]
    SEXP_5=[exp_E[13]+exp_E[15],exp_E[5]+exp_E[11],exp_E[11]+exp_E[15],exp_E[5]+exp_E[13],exp_E[3]+exp_E[8],exp_E[3]+exp_E[9],exp_E[3]+exp_E[10],exp_E[1]+exp_E[8],exp_E[1]+exp_E[9],exp_E[1]+exp_E[10],exp_E[8]+exp_E[12],exp_E[9]+exp_E[12],exp_E[10]+exp_E[12],exp_E[8]+exp_E[14],exp_E[9]+exp_E[14],exp_E[10]+exp_E[14]]
    SEXP_6=[exp_E[1]+exp_E[8],exp_E[1]+exp_E[9],exp_E[1]+exp_E[10],exp_E[3]+exp_E[8],exp_E[3]+exp_E[9],exp_E[3]+exp_E[10],exp_E[8]+exp_E[12],exp_E[9]+exp_E[12],exp_E[10]+exp_E[12],exp_E[8]+exp_E[14],exp_E[9]+exp_E[14],exp_E[10]+exp_E[14]]
    for i in range(16):
        SEXP_1[i]=(SEXP_1[i]/Z_)
        SEXP_2[i]=(SEXP_2[i]/Z_)
        SEXP_4[i]=(SEXP_4[i]/Z_)
        SEXP_5[i]=(SEXP_5[i]/Z_)
        if i<12:
             SEXP_3[i]=(SEXP_3[i]/Z_)
             SEXP_6[i]=(SEXP_6[i]/Z_)
#for i in range(0,16):
#    print i+1,":",E[i]
#def c(T,E):
#    cv=0.0
#    ay=0.0
#    z=0.0
#    for Ei in E:
#        cv+=Ei**2*exp(-Ei/T)
#        ay+=Ei*exp(-Ei/T)
#        z+=exp(-Ei/T)
#    return (cv-ay**2/z)/(z*T**2)
#for i in range(1000,10000):
#    print i*0.01, c(i*0.01,E)
#########################################################################
##############################Residuos###################################
#########################################################################
    phi=atan(2.0*V/(E_q-E_f+delta))
    theta=atan(2.0*V/(E_f+U-E_q-delta_p))
    a_9=a(E_9+E_min,E_f,E_q,U,V)
    a_10=a(E_10+E_min,E_f,E_q,U,V)
    a_11=a(E_11+E_min,E_f,E_q,U,V)
    a_=[a_9,a_10,a_11]
    b_9=b(E_9+E_min,E_f,E_q,U,V)
    b_10=b(E_10+E_min,E_f,E_q,U,V)
    b_11=b(E_11+E_min,E_f,E_q,U,V)
    b_=[b_9,b_10,b_11]
    c_9=c(E_9+E_min,E_f,E_q,U,V,E_9+E_min)
    c_10=c(E_10+E_min,E_f,E_q,U,V,E_10+E_min)
    c_11=c(E_11+E_min,E_f,E_q,U,V,E_11+E_min)
    c_=[c_9,c_10,c_11]
    return E,exp_E,phi,theta,a_,b_,c_,SEXP_1,SEXP_2,SEXP_3,SEXP_4,SEXP_5,SEXP_6
def gkondo_ss(E_q,V):
    S_s=gkondo_s(E_q,V)
    E=S_s[0]
    exp_E=S_s[1]
    phi=S_s[2]
    theta=S_s[3]
    a_=S_s[4]
    b_=S_s[5]
    c_=S_s[6]
    SEXP_1=S_s[7]
    SEXP_2=S_s[8]
    SEXP_3=S_s[9]
    SEXP_4=S_s[10]
    SEXP_5=S_s[11]
    SEXP_6=S_s[12]
    return E,exp_E,phi,theta,a_,b_,c_,SEXP_1,SEXP_2,SEXP_3,SEXP_4,SEXP_5,SEXP_6
def gkondo(w,V,E,exp_E,phi,theta,a_,b_,c_,SEXP_1,SEXP_2,SEXP_3,SEXP_4,SEXP_5,SEXP_6):
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
#Omega=-log(Z_)/beta # energia libre de Helmholtz F=-kBTln(Z)
##############################################################################
###############################g11############################################
##############################################################################
    g11=sin(phi)**2*(SEXP_1[0]/(w+E[0]-E[3])+1.5*SEXP_1[1]/(w+E[1]-E[5]))+cos(phi)**2*(SEXP_1[2]/(w+E[0]-E[1])+1.5*SEXP_1[3]/(w+E[3]-E[5]))
    g33=sin(theta)**2*(SEXP_2[0]/(w+E[14]-E[15])+1.5*SEXP_2[1]/(w+E[7]-E[11]))+cos(theta)**2*(SEXP_2[2]/(w+E[12]-E[15])+1.5*SEXP_2[3]/(w+E[7]-E[13]))
    g13=np.float64(0.0)
#    g22=sin(phi)**2*(SEXP_4[0]/(w+E[0]-E[4])+1.5*SEXP_4[1]/(w+E[2]-E[6]))+cos(phi)**2*(SEXP_4[2]/(w+E[0]-E[2])+1.5*SEXP_4[3]/(w+E[4]-E[6]))
#    g44=sin(theta)**2*(SEXP_5[0]/(w+E[13]-E[15])+1.5*SEXP_5[1]/(w+E[5]-E[11]))+cos(theta)**2*(SEXP_5[2]/(w+E[11]-E[15])+1.5*SEXP_5[3]/(w+E[5]-E[13]))
#    g24=np.float64(0.0)
    for i in range(8,11):
        g11+=(SEXP_1[i-4]/(w+E[2]-E[i]))*(a_[i-8]*sin(phi))**2+(SEXP_1[i-1]/(w+E[4]-E[i]))*(a_[i-8]*cos(phi))**2+(SEXP_1[i+2]/(w+E[i]-E[13]))*(c_[i-8]*sin(theta))**2+(SEXP_1[i+5]/(w+E[i]-E[11]))*(c_[i-8]*cos(theta))**2
        g33+=(SEXP_2[i-4]/(w+E[3]-E[i]))*(b_[i-8]*sin(phi))**2+(SEXP_2[i-1]/(w+E[2]-E[i]))*(b_[i-8]*cos(phi))**2+(SEXP_2[i+2]/(w+E[i]-E[11]))*(a_[i-8]*sin(theta))**2+(SEXP_2[i+5]/(w+E[i]-E[13]))*(a_[i-8]*cos(theta))**2
        g13+=(SEXP_3[i-8]/(w+E[4]-E[i])-SEXP_3[i-5]/(w+E[2]-E[i]))*(a_[i-8]*b_[i-8]*sin(phi)*cos(phi))+(SEXP_3[i-2]/(w+E[i]-E[13])-SEXP_3[i+1]/(w+E[i]-E[11]))*(a_[i-8]*c_[i-8]*sin(theta)*cos(theta))
#        g22+=(SEXP_4[i-4]/(w+E[1]-E[i]))*(a_[i-8]*sin(phi))**2+(SEXP_4[i-1]/(w+E[3]-E[i]))*(a_[i-8]*cos(phi))**2+(SEXP_4[i+2]/(w+E[i]-E[14]))*(c_[i-8]*sin(theta))**2+(SEXP_4[i+5]/(w+E[i]-E[12]))*(c_[i-8]*cos(theta))**2
#        g44+=(SEXP_5[i-4]/(w+E[3]-E[i]))*(b_[i-8]*sin(phi))**2+(SEXP_5[i-1]/(w+E[1]-E[i]))*(b_[i-8]*cos(phi))**2+(SEXP_5[i+2]/(w+E[i]-E[12]))*(a_[i-8]*sin(theta))**2+(SEXP_5[i+5]/(w+E[i]-E[14]))*(a_[i-8]*cos(theta))**2
#        g24+=(SEXP_6[i-8]/(w+E[1]-E[i])-SEXP_6[i-5]/(w+E[3]-E[i]))*(a_[i-8]*b_[i-8]*sin(phi)*cos(phi))+(SEXP_6[i-2]/(w+E[i]-E[12])-SEXP_6[i+1]/(w+E[i]-E[14]))*(a_[i-8]*c_[i-8]*sin(theta)*cos(theta))
    g31=g13
    g11=-g11
    g33=-g33
    g13=-g13
    g31=-g31
#    g42=g24
    #g11=g11.evalf()
#print g11.subs(w,0)
##############################################################################
###############################g33############################################
##############################################################################
        #g33=g33.evalf()
#print g33.subs(w,3)
##############################################################################
###############################g13############################################
##############################################################################        
        #g13=g13.evalf()
#print g13.subs(w,3)
##############################################################################
###############################g31############################################
##############################################################################
##############################################################################
###############################g22############################################
##############################################################################
        #g22=g22.evalf()
#print g22.subs(w,3)
##############################################################################
###############################g44############################################
##############################################################################            
        #g44=g44.evalf()
#print g44.subs(w,3)
##############################################################################
###############################g24############################################
##############################################################################        
        #g24=g24.evalf()
#print g24.subs(w,3)
##############################################################################
###############################g42############################################
##############################################################################
##############################################################################
##########################M_up################################################
##############################################################################
    phi_o_sigma=(-1/(w-E_q+mu))
    #den_=1+delta_0**2*phi_o_sigma*(g11+g33+g13+g31)
    #m_11_up=(g11+delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
    #m_13_up=(g13-delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
    #m_31_up=(g31-delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
    #m_33_up=(g33+delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
##############################################################################
##########################M_down##############################################
##############################################################################
    #den_0=1+delta_0**2*phi_o_sigma*(g22+g44-g24-g42)
    #m_22_down=(g22+delta_0**2*phi_o_sigma*(g22*g44-g24*g42))/den_0
    #m_24_down=(g24+delta_0**2*phi_o_sigma*(g22*g44-g24*g42))/den_0
    #m_42_down=(g42+delta_0**2*phi_o_sigma*(g22*g44-g24*g42))/den_0
    #m_44_down=(g44+delta_0**2*phi_o_sigma*(g22*g44-g24*g42))/den_0
##############################################################################
##########################G_up################################################
##############################################################################
    phi_sigma=log((w-B+mu)/(w-A+mu))/(2*D)
    #den_2=1-V**2*phi_sigma*(m_11_up+m_13_up+m_31_up+m_33_up)
    #G_11_up=(m_11_up-V**2*phi_sigma*(m_11_up*m_33_up-m_13_up*m_31_up))/den_2
    #G_13_up=(m_13_up+V**2*phi_sigma*(m_11_up*m_33_up-m_13_up*m_31_up))/den_2
    #G_31_up=(m_31_up+V**2*phi_sigma*(m_11_up*m_33_up-m_13_up*m_31_up))/den_2
    #G_33_up=(m_33_up-V**2*phi_sigma*(m_11_up*m_33_up-m_13_up*m_31_up))/den_2
##############################################################################
##########################G_down##############################################
##############################################################################
    #den_3=1-V**2*phi_sigma*(m_22_down+m_44_down-m_24_down-m_42_down)
    #G_22_down=(m_22_down-V**2*phi_sigma*(m_22_down*m_44_down-m_24_down*m_42_down))/den_3
    #G_24_down=(m_24_down-V**2*phi_sigma*(m_22_down*m_44_down-m_24_down*m_42_down))/den_3
    #G_42_down=(m_42_down-V**2*phi_sigma*(m_22_down*m_44_down-m_24_down*m_42_down))/den_3
    #G_44_down=(m_44_down-V**2*phi_sigma*(m_22_down*m_44_down-m_24_down*m_42_down))/den_3
############ Gcc #########
    #gamma_a=m_11_up+m_13_up+m_31_up+m_33_up
    #gamma_b=m_22_down+m_44_down-m_24_down-m_42_down
    #A_up=-w-V**2*gamma_a
    #A_down=-w-V**2*gamma_b
    #G_cc_up=(1/(2*D))*log((A_up+D-mu)/(A_up-D-mu))
    #G_cc_down=(1/(2*D))*log((A_down+D-mu)/(A_down-D-mu))
##########################
    #G_ff_at=0.0
    #for k in range(0,16):
    #    G_ff_at=G_ff_at+M[k]/(w-u_p[k])
    #G_ff_at=G_ff_at/Z_
#!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!! M_at !!!!!!!!!!
    #M_at=G_ff_at/(1.0+G_ff_at*delta_0**2*(-1.0/(w-E_q)))
#!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!! G_ff_imp !!!!!!
    #G_ff_imp=-M_at/(1.0+M_at*(V**2/(2.0*D))*log((w+D+mu)/(w-D+mu)))
#!!!!!!!!!!!!!!!!!!!!!!!!!
    #n_F=1/(1+exp(beta*real(w)))
# funcion de green up
    #G_F=-(G_11_up+G_33_up+G_31_up+G_13_up).evalf()
    sg=g11+g33+g31+g13
    det_g=g11*g33-g13*g31
    K_=1+delta_0**2*phi_o_sigma*sg
    V=sqrt(delta_0*2.0*D/pi)
    S_=1-V**2*phi_sigma*sg/K_
    G_F=(-1/(1/sg+delta_0**2*phi_o_sigma-V**2*phi_sigma))
    G_13_up=(g13+det_g*(V**2*phi_sigma/K_-delta_0**2*phi_o_sigma*S_))/(1+delta_0**2*phi_o_sigma*sg-V**2*phi_sigma*sg)
    #print re(w),(-1/pi)*im(G_F)#,(-1/pi)*im(G_F_n)
#   GF2=w/(w*(w-E_f)-V**2)
#   GF3=GF2*V**2/w**2
#   GF2=GF2+GF3
#GF=GF-(G_22_down+G_44_down)#+G_24_down+G_42_down)
#I_=lambda sd:(-(-1/pi)*im((G_11_up.subs(w,sd+I*eta_)).evalf())*(n_F.subs(w,sd)).evalf()).evalf()
#print I_(0.0)
#print integrate.quad(I_,-D,D)
#w_=1e-6
#print (-1/pi)*im(GF.subs(w,w_+I*eta_).evalf())
    #!!!!!!!! numeros de ocupacion !!!!!!!!!
    return G_F,G_13_up
n_F_=[]
for j in range(n_puntos):
    w_=complex(np.float64(A_n+j*(B_n-A_n)/(n_puntos-1)),np.float64(eta_))
    n_F_.append(1/(1+exp(beta*real(w_))))
Eq_puntos=1
Eq_A=-1.0*delta_0
Eq_B=1.0*delta_0
for l in range(0,Eq_puntos):
    #E_q=np.float64((Eq_A+l*(Eq_B-Eq_A)/(Eq_puntos-1)))
    #m_1=cos(phi)**2*(exp(-beta*E_1)+exp(-beta*E_2)+1.5*exp(-beta*E_4)+1.5*exp(-beta*E_6))
    #m_2=sin(phi)**2*(exp(-beta*E_1)+exp(-beta*E_4)+1.5*exp(-beta*E_2)+1.5*exp(-beta*E_6))
    #m_3=(exp(-beta*E_3)+exp(-beta*E_10))*((a_10*sin(phi))**2+(b_10*cos(phi))**2)
    #m_4=(exp(-beta*E_3)+exp(-beta*E_11))*((a_11*sin(phi))**2+(b_11*cos(phi))**2)
    #m_5=(exp(-beta*E_3)+exp(-beta*E_9))*((a_9*sin(phi))**2+(b_9*cos(phi))**2)
    #m_6=(exp(-beta*E_4)+exp(-beta*E_10))*((a_10*cos(phi))**2+(b_10*sin(phi))**2)
    #m_7=sin(theta)**2*(1.5*(exp(-beta*E_8)+exp(-beta*E_12))+exp(-beta*E_15)+exp(-beta*E_16))
    #m_8=(exp(-beta*E_9)+exp(-beta*E_12))*((c_9*cos(theta))**2+(a_9*sin(theta))**2)
    #m_9=(exp(-beta*E_10)+exp(-beta*E_12))*((c_10*cos(theta))**2+(a_10*sin(theta))**2)
    #m_10=(exp(-beta*E_11)+exp(-beta*E_12))*((c_11*cos(theta))**2+(a_11*sin(theta))**2)
    #m_11=(exp(-beta*E_10)+exp(-beta*E_14))*((c_10*sin(theta))**2+(a_10*cos(theta))**2)
    #m_12=(exp(-beta*E_11)+exp(-beta*E_15))*((c_11*sin(theta))**2+(a_11*cos(theta))**2)
    #m_13=(exp(-beta*E_5)+exp(-beta*E_9))*((a_9*cos(phi))**2+(b_9*sin(phi))**2)
    #m_14=(exp(-beta*E_5)+exp(-beta*E_11))*((a_11*cos(phi))**2+(b_11*sin(phi))**2)
    #m_15=cos(theta)**2*(1.5*(exp(-beta*E_8)+exp(-beta*E_14))+exp(-beta*E_13)+exp(-beta*E_16))
    #m_16=(exp(-beta*E_9)+exp(-beta*E_15))*((a_9*cos(theta))**2+(c_9*sin(theta))**2)
    #M=[m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,m_13,m_14,m_15,m_16]
#########################################################################
###############################polos#####################################
#########################################################################
    #u_1=E_2-E_1
    #u_2=E_5-E_1
    #u_3=E_10-E_2
    #u_4=E_11-E_2
    #u_5=E_9-E_2
    #u_6=E_10-E_4
    #u_7=E_12-E_6
    #u_8=E_12-E_9
    #u_9=E_12-E_10
    #u_10=E_12-E_11
    #u_11=E_14-E_10
    #u_12=E_14-E_11
    #u_13=E_9-E_4
    #u_14=E_11-E_4
    #u_15=E_14-E_6
    #u_16=E_14-E_9
    #u_p=[u_1,u_2,u_3,u_4,u_5,u_6,u_7,u_8,u_9,u_10,u_11,u_12,u_13,u_14,u_15,u_16]
    #!!!!!!!! numeros de ocupacion !!!!!!!!                
    #n_f_00_1=0.0
    n_f_s=np.float64(0.0)
    n_f_pp_1=np.float64(0.0)
    #n_f_mm_1=0.0
    #n_f_dd_1=0.0
    #n_f_00_2=0.0
    #n_f_mm_2=0.0
    #n_f_pp_2=0.0
    #n_f_dd_2=0.0
    S_s_00=gkondo_ss(E_q,delta_0)
    E=S_s_00[0]
    exp_E=S_s_00[1]
    phi=S_s_00[2]
    theta=S_s_00[3]
    a_=S_s_00[4]
    b_=S_s_00[5]
    c_=S_s_00[6]
    SEXP_1=S_s_00[7]
    SEXP_2=S_s_00[8]
    SEXP_3=S_s_00[9]
    SEXP_4=S_s_00[10]
    SEXP_5=S_s_00[11]
    SEXP_6=S_s_00[12]
    #S_s_0=gkondo(x,delta_0,E,exp_E,phi,theta,a_,b_,c_,SEXP_1,SEXP_2,SEXP_3,SEXP_4,SEXP_5,SEXP_6)
    #SS=lambdify(x,S_s_0,"numpy")
    for j in range(0,n_puntos):
        w=complex(A_n+j*(B_n-A_n)/(n_puntos-1),eta_)
        S=gkondo(w,delta_0,E,exp_E,phi,theta,a_,b_,c_,SEXP_1,SEXP_2,SEXP_3,SEXP_4,SEXP_5,SEXP_6)
        G_F=S[0]
        G_13_up=S[1]
        #S_=SS(w)
        #G_F=S_[0]
        #G_13_up=S_[1]
        G_F_n_=G_F+2.0*G_13_up
        #G_11_up=lambdify((x,E_q_),S_s_0[1],"numpy")
        #G_33_up=lambdify((x,E_q_),S_s_0[2],"numpy")
        print real(w), (-1/pi)*imag(G_F)
        n_F=n_F_[j]
        if j==0 or j==n_puntos-1:
            #n_f_00_1=n_f_00_1+(-1.0/pi)*im(G_11_up.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            n_f_s=n_f_s+(-1.0/pi)*imag(G_F)*n_F
            n_f_pp_1=n_f_pp_1+(-1.0/pi)*imag(G_F_n_)
            #n_f_mm_1=n_f_mm_1+(-1.0/pi)*im(G_33_up.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_dd_1=n_f_dd_1+(-1.0/pi)*im(G_33_up.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #n_f_00_2=n_f_00_2+(-1.0/pi)*im(G_22_down.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_mm_2=n_f_mm_2+(-1.0/pi)*im(G_22_down.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #n_f_pp_2=n_f_pp_2+(-1.0/pi)*im(G_44_down.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_dd_2=n_f_dd_2+(-1.0/pi)*im(G_44_down.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #print j
        elif j%2==0 and j!=0 and j!=n_puntos-1:
            #n_f_00_1=n_f_00_1+2.0*(-1.0/pi)*im(G_11_up.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            n_f_s=n_f_s+2.0*(-1.0/pi)*imag(G_F)*n_F
            n_f_pp_1=n_f_pp_1+2.0*(-1.0/pi)*imag(G_F_n_)
            #n_f_mm_1=n_f_mm_1+2.0*(-1.0/pi)*im(G_33_up.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_dd_1=n_f_dd_1+2.0*(-1.0/pi)*im(G_33_up.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #n_f_00_2=n_f_00_2+2.0*(-1.0/pi)*im(G_22_down.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_mm_2=n_f_mm_2+2.0*(-1.0/pi)*im(G_22_down.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #n_f_pp_2=n_f_pp_2+2.0*(-1.0/pi)*im(G_44_down.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_dd_2=n_f_dd_2+2.0*(-1.0/pi)*im(G_44_down.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #print j
        elif j%2==1 and j!=0 and j!=n_puntos-1:
            #n_f_00_1=n_f_00_1+4.0*(-1.0/pi)*im(G_11_up.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            n_f_s=n_f_s+4.0*(-1.0/pi)*imag(G_F)*n_F
            n_f_pp_1=n_f_pp_1+4.0*(-1.0/pi)*imag(G_F_n_)
            #n_f_mm_1=n_f_mm_1+4.0*(-1.0/pi)*im(G_33_up.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_dd_1=n_f_dd_1+4.0*(-1.0/pi)*im(G_33_up.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #n_f_00_2=n_f_00_2+4.0*(-1.0/pi)*im(G_22_down.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_mm_2=n_f_mm_2+4.0*(-1.0/pi)*im(G_22_down.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #n_f_pp_2=n_f_pp_2+4.0*(-1.0/pi)*im(G_44_down.subs(w,w_+I*eta_).evalf())*(1.0-n_F).subs(w,w_).evalf()
            #n_f_dd_2=n_f_dd_2+4.0*(-1.0/pi)*im(G_44_down.subs(w,w_+I*eta_).evalf())*n_F.subs(w,w_).evalf()
            #print j
        if abs(real(w))<=0.5*(B_n-A_n)/(n_puntos-1):
            pho_fsmu=(-1.0/pi)*imag(G_F)
    #n_f_00_1=n_f_00_1*(B_n-A_n)/(3.0*n_puntos)
    n_f_s=n_f_s*(B_n-A_n)/(3.0*n_puntos)
    n_f_pp_1=n_f_pp_1*(B_n-A_n)/(3.0*n_puntos)
    n_f_pp_1=np.float64(n_f_pp_1)
    #n_f_mm_1=n_f_mm_1*(B_n-A_n)/(3.0*n_puntos)
    #n_f_dd_1=n_f_dd_1*(B_n-A_n)/(3.0*n_puntos)
    #n_f_00_2=n_f_00_2*(B_n-A_n)/(3.0*n_puntos)
    #n_f_mm_2=n_f_mm_2*(B_n-A_n)/(3.0*n_puntos)
    #n_f_pp_2=n_f_pp_2*(B_n-A_n)/(3.0*n_puntos)
    #n_f_dd_2=n_f_dd_2*(B_n-A_n)/(3.0*n_puntos)
    #n_f_tot_1=n_f_00_1+n_f_pp_1+n_f_mm_1+n_f_dd_1
    #n_f_tot_2=n_f_00_2+n_f_pp_2+n_f_mm_2+n_f_dd_2
    F_Eq=(sin(pi*n_f_s)**2/(delta_0*pi))-pho_fsmu
    print E_q/delta_0,n_f_s,n_f_pp_1,pho_fsmu
    #if abs(F_Eq)<fmin:
    #    fmin=abs(F_Eq)
    #    Eqmin=E_q
    #    pho_fsmumin=pho_fsmu
    if abs(n_f_pp_1-1.0)<1e-6:
        fmin=abs(F_Eq)
        Eqmin=E_q
        pho_fsmumin=pho_fsmu
        n_f_pp_1min=n_f_pp_1
#print Eqmin/delta_0,fmin,n_f_pp_1min, pho_fsmumin