          program densidad
          implicit none    
          real(8) pi,eta,delta_0,E_f,E_q0,E_q1,E_q,U,A,B,D,V,V_0,T,beta,mu
          real(8) x1,x2,a_2_p,delta,delta_p,Q,R,theta_1,Eq_AI,Eq_BI 
          real(8) E_1,E_2,E_3,E_4,E_5,E_6,E_7,E_8,E_9,E_10,E_11,E_12
          real(8) E_13,E_14,E_15,E_16,E_min
          real(8) phi,theta,a_9,a_10,a_11,b_9,b_10,b_11,c_9,c_10,c_11
          complex(8) g11,g33,g13,g31,g22,g44,g24,g42,g11_n,g33_n,g13_n,&
          g31_n,m_11_up_n,m_13_up_n,m_31_up_n,m_33_up_n,phi_sigma_n,phi_o_sigma_n
          complex(8) g11_mu,g33_mu,g13_mu,g31_mu
          complex(8) g44_n
          real(8), dimension(16) :: E
          real(8), dimension(3) :: a_,b_,c_
          real(8) a_0,b_0,c_0,Z_,T_A,T_B,T_n,T_0,T_0_n
          real(8) A_n,B_n,Eq_A,Eq_B,A_I,B_I
          real(8) n_F,diff_n_F,E_q0_,E_q1_,E_q_
          complex(8) phi_sigma,phi_o_sigma,den_,den_mu,den_n,den_0,den_2,den_3
          complex(8) den_2_n,den_2_mu
          complex(8) phi_sigma_mu,phi_o_sigma_mu
          complex(8) m_11_up,m_13_up,m_31_up,m_33_up,m_22_down,m_24_down,m_42_down,m_44_down&
          ,G_11_up,G_13_up,G_31_up,G_33_up,G_22_down,G_24_down,G_42_down,G_44_down
          complex(8) G_11_up_n,G_13_up_n,G_31_up_n,G_33_up_n
          complex(8) m_11_up_mu,m_13_up_mu,m_31_up_mu,m_33_up_mu
          integer i,j,k,l,n_puntos,Eq_puntos,E_f_puntos,U_puntos,co_
          integer p_,EU_,hel_,T_puntos
          real(8) G_,S_,gamma_,c2e2dhb,C_v_lim,V_E,V_E_0,V_E_n,V_E_lim
          real(8) f_E,f_U
          real(8) n_cc_up
          complex(8) G_00_s,G_F_mu,G_cc_up_mu,w_mu
          complex(8) w,G_F,w_n,sg_mu
          complex(8) G_00_11,G_00_13,G_00_31,G_00_33,sg,G_F_n,G_F_n_2,G_F_2
          complex(8) G_F_G,sg_n,G_cc_up,G_cf_0,G_cf_s
          real(8) m_1,m_2,m_3,m_4,m_5,m_6,m_7,m_8,m_9,m_10,m_11,m_12,&
          m_13,m_14,m_15,m_16
          real(8), dimension(16) :: M
          real(8) u_1,u_2,u_3,u_4,u_5,u_6,u_7,u_8,u_9,u_10,u_11,u_12,&
          u_13,u_14,u_15,u_16
          real(8), dimension(16) :: u_p
          complex(8) G_ff_at, M_at, G_ff_imp
          real(8) n_f_00_1,n_f_pp_1,n_f_mm_1,n_f_dd_1
          real(8) n_f_00_2,n_f_pp_2,n_f_mm_2,n_f_dd_2
          real(8) n_f_tot_1,n_f_tot_2,n_f_s
          real(8) F_Eq0,F_Eq1,F_Eq,F_Eq_n,F_Eq_n_0,F_Eq_n_1 ! funcion para calcular E_q tal que la regla de suma de friedel se cumpla
          real(8) pho_fsmu,pho_csmu,tol,tol_n,E_fA,E_fB,U_A,U_B,Fmin_,N_s
          character(20) str
          external Fmin_,F_Eq_n
          common/VARS1/pi,V,E_f,U,x1,x2,delta_0,A,B,mu
          common/VARS2/w_mu,beta,D,A_n,B_n,eta,n_puntos
          common/VARS3/E_q,T_0
          common/VARS4/N_s
          common/VARS5/E_fA,E_fB,E_f_puntos
          common/VARS6/U_A,U_B,U_puntos
          common/T_/T
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!! PARAMETROS DEL PROBLEMA !!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          pi=dacos(-1.0d0)
          EU_=0 ! si EU_=0, entonces U es constante y E varia, si EU_=1,
          !!! entonces E es constante y U varia. Cuando EU_=0 hace la
          !!! grafica de conductancia G vs Ef
          N_s=2.5d0 !! numero particulas banda conduccion?
          delta_0=1.0d-2
          E_f=-10.0d0*delta_0 ! Enegia del estado localizado
          !E_q=0.0d0*delta_0
          U=20.0d0*delta_0 ! Repulsion columbiana
          D=100.0d0*delta_0
          T=1.0d-3*delta_0 ! Temperatura
          beta=1.0d0/T
          eta=1.0d-18
          mu=0.0d0 !!!!! potencial quimico
          n_puntos=10000 !!!!! cantidad de puntos para la densidad
          !!!!! de estados, debe ser un numero par debido al criterio de
          !!!!! integracion por regla de simpson
          tol=1.0d-5 ! tolerancia al calcular E_q y tambien la incertidumbre de E_q, no de F_Eq
          !!!!!!!!!intervalo para buscar el Eq talque FSR se cumpla,
          !donde Eq_puntos es la cantidad de puntos a buscar entre Eq_A
          !y Eq_B!!!!!!!!!!!
          Eq_A=-1.0d0*delta_0 
          Eq_B=1.0d0*delta_0
          Eq_puntos=10 !!!!! E_q_puntos-1 intervalos para encontrar
          !!!!!!!!!!!!!!!!!! E_qmin
          E_fA=-30.0d0*delta_0 !!!!!!! intervalo en Eq cuando este varia 
          !!!!! y U es constante, con E_f_puntos la cantidad de puntos
          E_fB=5.0d0*delta_0
          E_f_puntos=36 !!!!! cantidad de archivos para calcular DOS
          !!! Eq_min
          !!!!!!!!!!!!!!!!!!!!
          U_A=2.0d0*delta_0 !!!!!!! intervalo en U cuando este varia 
          !!!!! y Ef es constante, con U_puntos la cantidad de puntos
          U_B=30.0d0*delta_0
          U_puntos=33 !!!!! cantidad de archivos para calcular DOS
          T_A=T !!! temperatura inicial
          T_B=100.0d0*delta_0 !!! temperatura final
          T_puntos=1000 !!! cantidad de puntos para temperatura
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!! NO EDITAR A PARTIR DE ESTE PUNTO !!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          V=dsqrt(delta_0*2.0d0*D/pi)  ! Hibridacion  
          gamma_=V**2/delta_0
          !##########Autovalores##########
          T_0=T
          A_I=-300.0d0*D*T !!! se elige asi el intervalo para que al
          B_I=300.0d0*D*T !!! hacer el producto beta*w, el intervalo nuevo este entre 
          !!! [-100,100] y asi el resultado no es nulo y la integral no es
          !!! cero y la particion es buena y fina
          A_n=-D !!!! intervalo [A_n,B_n] de w, puede variar, pero es
          ! el intervalo de integracion en los numeros de ocupacion
          ! asi que debe ser igual a -D
          !!!! alrededor del potencial quimico, aqui mu=0
          B_n=D ! intervalo [A,B] de w, puede variar, pero es
          ! el intervalo de integracion en los numeros de ocupacion
          ! asi que debe ser igual a D
          tol_n=tol !!!! variable de la tolerancia auxiliar
          w_mu=DCMPLX(mu,eta) !!!! w=mu para el calculo de pho(w=mu)
          A=-D ! intervalo [A,B] de la banda de conduccion
          B=D ! intervalo [A,B] de la banda de conduccion
          if (EU_==0) then
              open(7,file="U="//trim(str(int(U/delta_0)))//"."//&
              trim(str(int(10000*abs(U/delta_0-int(U/delta_0)))))//&
              "_G_vs_Ef.txt") !!!! si se llama
              !!!! open(6,file="U=20_G_vs_Ef.txt") por alguna
              !!!! razon falla y no avanza el programa
              open(5,file="U="//trim(str(int(U/delta_0)))//"."//&
              trim(str(int(10000*abs(U/delta_0-int(U/delta_0)))))//&
              "_pho_vs_Ef.txt")
              open(4,file="U="//trim(str(int(U/delta_0)))//"."//&
              trim(str(int(10000*abs(U/delta_0-int(U/delta_0)))))//&
              "_nf_vs_Ef.txt")
          end if
          if (EU_==1) then
              open(10,file="Ef="//trim(str(int(E_f/delta_0)))//"."//&
              trim(str(int(10000*abs(E_f/delta_0-int(E_f/delta_0)))))//&
              "_G_vs_U.txt") 
              open(3,file="Ef="//trim(str(int(E_f/delta_0)))//"."//&
              trim(str(int(10000*abs(E_f/delta_0-int(E_f/delta_0)))))//&
              "_pho_vs_U.txt") 
              open(2,file="Ef="//trim(str(int(E_f/delta_0)))//"."//&
              trim(str(int(10000*abs(E_f/delta_0-int(E_f/delta_0)))))//&
              "_nf_vs_U.txt") 
          end if
          !!!!!!! cuando se quiera hacer una grafica con Ef o U constantes
          !!!!!!! simplemente comentar la linea E_f o U
          if (EU_==0) then
              U_puntos=1
          end if
          if (EU_==1) then
              E_f_puntos=1
          end if
          do p_=0,E_f_puntos-1
              E_f=f_E(p_,EU_)
          do k=0,U_puntos-1
              U=f_U(k,EU_)
          !!!! calculo de Eq usando la funcion Fmin_
              E_q0=Fmin_(Eq_A,Eq_A+(Eq_B-Eq_A)/(Eq_puntos-1),F_Eq_n,tol)
              E_q0_=E_q0
              F_Eq_n_0=F_Eq_n(E_q0)
          do l=2,Eq_puntos-1 !!!! l empieza desde 2 porque en la linea
          !!!! anterior ya se calcula E_q en el 1er intervalo
              Eq_AI=Eq_A+(l-1)*(Eq_B-Eq_A)/(Eq_puntos-1)
              Eq_BI=Eq_A+l*(Eq_B-Eq_A)/(Eq_puntos-1)
              E_q1=Fmin_(Eq_AI,Eq_BI,F_Eq_n,tol)
              F_Eq_n_1=F_Eq_n(E_q1)
              if (F_Eq_n_1<F_Eq_n_0) then
                 E_q=E_q1
                 E_q0=E_q1
                 F_Eq_n_0=F_Eq_n_1
              else
                 E_q=E_q0
                 F_Eq_n_0=F_Eq_n(E_q0)
              end if
          end do
          open(13,file="C_v_vs_T_Ef="//trim(str(int(E_f/delta_0)))//"."//&
          trim(str(int(10000*abs(E_f/delta_0-int(E_f/delta_0)))))//& 
          ":U="//trim(str(int(U/delta_0)))//"."//&
          trim(str(int(10000*abs(U/delta_0-int(U/delta_0)))))//& 
          ".txt")
          open(11,file="C_v_vs_T_LIMITE_Ef="//trim(str(int(E_f/delta_0)))//"."//&
          trim(str(int(10000*abs(E_f/delta_0-int(E_f/delta_0)))))//& 
          ":U="//trim(str(int(U/delta_0)))//"."//&
          trim(str(int(10000*abs(U/delta_0-int(U/delta_0)))))//& 
          ".txt")
          open(9,file="V_E_vs_LIMITE_T_Ef="//trim(str(int(E_f/delta_0)))//"."//&
          trim(str(int(10000*abs(E_f/delta_0-int(E_f/delta_0)))))//& 
          ":U="//trim(str(int(U/delta_0)))//"."//&
          trim(str(int(10000*abs(U/delta_0-int(U/delta_0)))))//& 
          ".txt")
          open(8,file="V_E_vs_T_Ef="//trim(str(int(E_f/delta_0)))//"."//&
          trim(str(int(10000*abs(E_f/delta_0-int(E_f/delta_0)))))//& 
          ":U="//trim(str(int(U/delta_0)))//"."//&
          trim(str(int(10000*abs(U/delta_0-int(U/delta_0)))))//& 
          ".txt")
          T_0_n=T_0
          V_E_0=V_E(T_0_n)
          write(8,*) T_0_n,V_E_0
          write(9,*) T_0_n,V_E_lim(T_0_n)
          write(13,*) T_0,C_v_lim(T_0_n)
          do l=1,T_puntos-1
              T_n=T_A+l*(T_B-T_A)/(T_puntos-1)
              V_E_n=V_E(T_n)
              T_n=T_A+l*(T_B-T_A)/(T_puntos-1)
              write(8,*) T_n,V_E_n
              write(13,*) T_n,(V_E_n-V_E_0)/(T_n-T_0_n)
              T_0_n=T_n
              V_E_0=V_E_n
              write(9,*) T_n,V_E_lim(T_n)
              T_n=T_A+l*(T_B-T_A)/(T_puntos-1)
              write(11,*) T_n,C_v_lim(T_n)
          end do
          close(8)
          close(9)
          close(11)
          close(13)
          T_0_n=T_0
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          open(1,file="Ef="//trim(str(int(E_f/delta_0)))//"."//&
          trim(str(int(10000*abs(E_f/delta_0-int(E_f/delta_0)))))//& 
          ":U="//trim(str(int(U/delta_0)))//"."//&
          trim(str(int(10000*abs(U/delta_0-int(U/delta_0)))))//& 
          ".txt")
          x1=E_f
          x2=E_q
          a_2_p=-(3*x1+3*x2+U)
          V=delta_0
          delta=dsqrt((x1-x2)**2+4*V**2)
          delta_p=dsqrt((x1+U-x2)**2+4*V**2)
          Q=-(12*V**2+(x1+x2)**2+(2*x1+U)**2+(2*x2)**2-&
          (x2+x1)*(2*x1+U)-(x2+x1)*(2*x2)-(2*x1+U)*(2*x2))/9.0
          R=(-3*((x2+x1)**2*(2*x1+U)+(x2+x1)**2*(2*x2)+&
          (2*x1+U)**2*(x2+x1)+(2*x1+U)**2*(2*x2)+(2*x2)**2*&
          (x2+x1)+(2*x2)**2*(2*x1+U))+12*(x2+x1)*(2*x1+U)*&
          (2*x2)+18*V**2*(2*(x2+x1)-(2*x1+U)-2*x2)+2*&
          ((x2+x1)**3+(2*x1+U)**3+(2*x2)**3))/54.0
          theta_1=acos(R/sqrt(-Q**3))
          E_1=0.0
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
          E(1)=E_1
          E(2)=E_2
          E(3)=E_3
          E(4)=E_4
          E(5)=E_5
          E(6)=E_6
          E(7)=E_7
          E(8)=E_8
          E(9)=E_9
          E(10)=E_10
          E(11)=E_11
          E(12)=E_12
          E(13)=E_13
          E(14)=E_14
          E(15)=E_15
          E(16)=E_16
          E_min=minval(E)
          phi=atan(2.0*V/(E_q-E_f+delta))
          theta=atan(2.0*V/(E_f+U-E_q-delta_p))
          a_9=a_0(E_9,E_f,E_q,U,V)
          a_10=a_0(E_10,E_f,E_q,U,V)
          a_11=a_0(E_11,E_f,E_q,U,V)
          a_(1)=a_9
          a_(2)=a_10
          a_(3)=a_11
          b_9=b_0(E_9,E_f,E_q,U,V)
          b_10=b_0(E_10,E_f,E_q,U,V)
          b_11=b_0(E_11,E_f,E_q,U,V)
          b_(1)=b_9
          b_(2)=b_10
          b_(3)=b_11
          c_9=c_0(E_9,E_f,E_q,U,V,E_9)
          c_10=c_0(E_10,E_f,E_q,U,V,E_10)
          c_11=c_0(E_11,E_f,E_q,U,V,E_11)
          c_(1)=c_9
          c_(2)=c_10
          c_(3)=c_11
          !!!!!!!! Funcion de particion !!!!
          Z_=0.0d0
          do i=1,16
              Z_=Z_+dexp(-beta*(E(i)-E_min))
          end do
          !!!!!!!! Numeros de ocupacion !!!!!!!!
          n_f_s=0.0d0
          n_cc_up=0.0d0
          n_f_00_1=0.0d0
          n_f_pp_1=0.0d0
          n_f_mm_1=0.0d0
          n_f_dd_1=0.0d0
          n_f_00_2=0.0d0
          n_f_mm_2=0.0d0
          n_f_pp_2=0.0d0
          n_f_dd_2=0.0d0
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!! Conductancia !!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          G_=0.0d0
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
          do i=1,16
              E(i)=E(i)-E_min
          end do
          !!!! g11_mu g11(w=mu)
          g11_mu=(sin(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_4))&
          /(w_mu+E_1-E_4)+&
          1.5*(dexp(-beta*E_2)+dexp(-beta*E_6))/(w_mu+E_2-E_6))+&
          cos(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_2))/(w_mu+E_1-E_2)&
          +1.5*(dexp(-beta*E_4)+dexp(-beta*E_6))/(w_mu+E_4-E_6)))    
          do i=9,11
          g11_mu=g11_mu+(((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_mu+E_3-E(i)))*&
          (a_(i-8)*sin(phi))**2+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w_mu+E_5-E(i)))*(a_(i-8)*cos(phi))**2+((dexp(-beta*E(i))+&
          dexp(-beta*E_14))/(w_mu+E(i)-E_14))*(c_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w_mu+E(i)-E_12))*&
          (c_(i-8)*cos(theta))**2)
          end do
          g11_mu=g11_mu/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          phi_sigma_mu=cdlog((w_mu-B+mu)/(w_mu-A+mu))/(2*D)
          phi_o_sigma_mu=-1/(w_mu-E_q-mu)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g33_mu g33(w=mu)
          g33_mu=(sin(theta)**2*((dexp(-beta*E_15)+dexp(-beta*E_16))/&
          (w_mu+E_15-E_16)+1.5*(dexp(-beta*E_8)+dexp(-beta*E_12))/&
          (w_mu+E_8-E_12))+cos(theta)**2*((dexp(-beta*E_13)+&
          dexp(-beta*E_16))/(w_mu+E_13-E_16)+1.5*(dexp(-beta*E_8)+&
          dexp(-beta*E_14))/(w_mu+E_8-E_14)))
          do i=9,11
          g33_mu=g33_mu+((dexp(-beta*E_4)+dexp(-beta*E(i)))/&
          (w_mu+E_4-E(i)))*&
          (b_(i-8)*sin(phi))**2+((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_mu+E_3-E(i)))*(b_(i-8)*cos(phi))**2+((dexp(-beta*E(i))&
          +dexp(-beta*E_12))/(w_mu+E(i)-E_12))*(a_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w_mu+E(i)-E_14))*&
          (a_(i-8)*cos(theta))**2
          end do
          g33_mu=g33_mu/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g13
          g13_mu=0.0d0
          do i=9,11
          g13_mu=g13_mu+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w_mu+E_5-E(i))-(dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_mu+E_3-E(i)))*(a_(i-8)*b_(i-8)*sin(phi)*cos(phi))+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w_mu+E(i)-E_14)-&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w_mu+E(i)-E_12)))*&
          (a_(i-8)*c_(i-8)*sin(theta)*cos(theta))
          end do
          g13_mu=g13_mu/Z_
          g31_mu=g13_mu
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          g11_mu=-g11_mu
          g33_mu=-g33_mu
          g13_mu=-g13_mu
          g31_mu=-g31_mu
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!! bucle para el calculo de la densidad de estados !!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do j=0,n_puntos-1
          w=DCMPLX((A_n+j*(B_n-A_n)/(n_puntos-1)),eta)
          w_n=DCMPLX((A_I+j*(B_I-A_I)/(n_puntos-1)),eta)
          n_F=1.0d0/(1+exp(beta*DREAL(w)))
          diff_n_F=-beta*exp(beta*DREAL(w_n))/&
                   ((1+exp(beta*DREAL(w_n)))**2)
          !!!! g11
          g11=(sin(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_4))&
          /(w+E_1-E_4)+&
          1.5*(dexp(-beta*E_2)+dexp(-beta*E_6))/(w+E_2-E_6))+&
          cos(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_2))/(w+E_1-E_2)+&
          1.5*(dexp(-beta*E_4)+dexp(-beta*E_6))/(w+E_4-E_6)))    
          do i=9,11
          g11=g11+(((dexp(-beta*E_3)+dexp(-beta*E(i)))/(w+E_3-E(i)))*&
          (a_(i-8)*sin(phi))**2+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w+E_5-E(i)))*(a_(i-8)*cos(phi))**2+((dexp(-beta*E(i))+&
          dexp(-beta*E_14))/(w+E(i)-E_14))*(c_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w+E(i)-E_12))*(c_(i-8)*&
          cos(theta))**2)
          end do
          g11=g11/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!! g11_n
          g11_n=(sin(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_4))&
          /(w_n+E_1-E_4)+&
          1.5*(dexp(-beta*E_2)+dexp(-beta*E_6))/(w_n+E_2-E_6))+&
          cos(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_2))/(w_n+E_1-E_2)+&
          1.5*(dexp(-beta*E_4)+dexp(-beta*E_6))/(w_n+E_4-E_6)))    
          do i=9,11
          g11_n=g11_n+(((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_n+E_3-E(i)))*&
          (a_(i-8)*sin(phi))**2+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w_n+E_5-E(i)))*(a_(i-8)*cos(phi))**2+((dexp(-beta*E(i))+&
          dexp(-beta*E_14))/(w_n+E(i)-E_14))*(c_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w_n+E(i)-E_12))*&
          (c_(i-8)*cos(theta))**2)
          end do
          g11_n=g11_n/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g33
          g33=(sin(theta)**2*((dexp(-beta*E_15)+dexp(-beta*E_16))/&
          (w+E_15-E_16)+1.5*(dexp(-beta*E_8)+dexp(-beta*E_12))/&
          (w+E_8-E_12))+cos(theta)**2*((dexp(-beta*E_13)+&
          dexp(-beta*E_16))/(w+E_13-E_16)+1.5*(dexp(-beta*E_8)+&
          dexp(-beta*E_14))/(w+E_8-E_14)))
          do i=9,11
          g33=g33+((dexp(-beta*E_4)+dexp(-beta*E(i)))/(w+E_4-E(i)))*&
          (b_(i-8)*sin(phi))**2+((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w+E_3-E(i)))*(b_(i-8)*cos(phi))**2+((dexp(-beta*E(i))&
          +dexp(-beta*E_12))/(w+E(i)-E_12))*(a_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w+E(i)-E_14))*&
          (a_(i-8)*cos(theta))**2
          end do
          g33=g33/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g33_n
          g33_n=(sin(theta)**2*((dexp(-beta*E_15)+dexp(-beta*E_16))/&
          (w_n+E_15-E_16)+1.5*(dexp(-beta*E_8)+dexp(-beta*E_12))/&
          (w_n+E_8-E_12))+cos(theta)**2*((dexp(-beta*E_13)+&
          dexp(-beta*E_16))/(w_n+E_13-E_16)+1.5*(dexp(-beta*E_8)+&
          dexp(-beta*E_14))/(w_n+E_8-E_14)))
          do i=9,11
          g33_n=g33_n+((dexp(-beta*E_4)+dexp(-beta*E(i)))/&
          (w_n+E_4-E(i)))*&
          (b_(i-8)*sin(phi))**2+((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_n+E_3-E(i)))*(b_(i-8)*cos(phi))**2+((dexp(-beta*E(i))&
          +dexp(-beta*E_12))/(w_n+E(i)-E_12))*(a_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w_n+E(i)-E_14))*&
          (a_(i-8)*cos(theta))**2
          end do
          g33_n=g33_n/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g13
          g13=0.0d0
          do i=9,11
          g13=g13+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w+E_5-E(i))-(dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w+E_3-E(i)))*(a_(i-8)*b_(i-8)*sin(phi)*cos(phi))+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w+E(i)-E_14)-&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w+E(i)-E_12)))*&
          (a_(i-8)*c_(i-8)*sin(theta)*cos(theta))
          end do
          g13=g13/Z_
          g31=g13
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g13_n
          g13_n=0.0d0
          do i=9,11
          g13_n=g13_n+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w_n+E_5-E(i))-(dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_n+E_3-E(i)))*(a_(i-8)*b_(i-8)*sin(phi)*cos(phi))+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w_n+E(i)-E_14)-&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w_n+E(i)-E_12)))*&
          (a_(i-8)*c_(i-8)*sin(theta)*cos(theta))
          end do
          g13_n=g13_n/Z_
          g31_n=g13_n
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g22
          g22=(sin(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_5))/&
          (w+E_1-E_5)+1.5*(dexp(-beta*E_3)+dexp(-beta*E_7))/&
          (w+E_3-E_7))+cos(phi)**2*((dexp(-beta*E_1)+&
          dexp(-beta*E_3))/(w+E_1-E_3)+1.5*(dexp(-beta*E_5)+&
          dexp(-beta*E_7))/(w+E_5-E_7)))
          do i=9,11
          g22=g22+((dexp(-beta*E_2)+dexp(-beta*E(i)))/&
          (w+E_2-E(i)))*(a_(i-8)*sin(phi))**2+&
          ((dexp(-beta*E_4)+dexp(-beta*E(i)))/(w+E_4-E(i)))*&
          (a_(i-8)*cos(phi))**2+((dexp(-beta*E(i))+&
          dexp(-beta*E_15))/(w+E(i)-E_15))*(c_(i-8)*&
          sin(theta))**2+((dexp(-beta*E(i))+dexp(-beta*E_13))/&
          (w+E(i)-E_13))*(c_(i-8)*cos(theta))**2
          end do
          g22=g22/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g44
          g44=(sin(theta)**2*((dexp(-beta*E_14)+dexp(-beta*E_16))/&
          (w+E_14-E_16)+1.5*(dexp(-beta*E_6)+dexp(-beta*E_12))/&
          (w+E_6-E_12))+cos(theta)**2*((dexp(-beta*E_12)+&
          dexp(-beta*E_16))/(w+E_12-E_16)+1.5*(dexp(-beta*E_6)+&
          dexp(-beta*E_14))/(w+E_6-E_14)))
          g44_n=0.0d0
          do i=9,11
          g44_n=g44_n+((dexp(-beta*E_4)+dexp(-beta*E(i)))/&
          (w+E_4-E(i)))*(b_(i-8)*sin(phi))**2+((dexp(-beta*E_2)+&
          dexp(-beta*E(i)))/(w+E_2-E(i)))*(b_(i-8)*cos(phi))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_13))/(w+E(i)-E_13))*&
          (a_(i-8)*sin(theta))**2+((dexp(-beta*E(i))+&
          dexp(-beta*E_15))/(w+E(i)-E_15))*(a_(i-8)*cos(theta))**2
          end do
          g44=g44+g44_n
          g44=g44/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g24
          g24=0.0
          do i=9,11
          g24=g24+((dexp(-beta*E_2)+dexp(-beta*E(i)))/(w+E_2-E(i))-&
          (dexp(-beta*E_4)+dexp(-beta*E(i)))/(w+E_4-E(i)))*(a_(i-8)*&
          b_(i-8)*sin(phi)*cos(phi))+((dexp(-beta*E(i))+&
          dexp(-beta*E_13))/&
          (w+E(i)-E_13)-((dexp(-beta*E(i))+dexp(-beta*E_15))/&
          (w+E(i)-E_15)))*&
          (a_(i-8)*c_(i-8)*sin(theta)*cos(theta))
          end do
          g24=g24/Z_
          g42=g24
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          g11=-g11
          g33=-g33
          g13=-g13
          g31=-g31
          g11_n=-g11_n
          g33_n=-g33_n
          g13_n=-g13_n
          g31_n=-g31_n
          g22=-g22
          g44=-g44
          g24=-g24
          g42=-g42
          V=dsqrt(delta_0*2.0d0*D/pi)  ! Hibridacion  
          phi_o_sigma=-1/(w-E_q-mu)
          phi_o_sigma_n=-1/(w_n-E_q-mu)
          den_=1+delta_0**2*phi_o_sigma*(g11+g33+g13+g31)
          den_n=1+delta_0**2*phi_o_sigma_n*(g11_n+g33_n+g13_n+g31_n)
          m_11_up=(g11+delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
          m_11_up_n=(g11_n+delta_0**2*phi_o_sigma_n*(g11_n&
          *g33_n-g13_n*g31_n))/den_n
          m_13_up=(g13-delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
          m_13_up_n=(g13_n+delta_0**2*phi_o_sigma_n*(g11_n&
          *g33_n-g13_n*g31_n))/den_n
          m_31_up=(g31-delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
          m_31_up_n=(g31_n+delta_0**2*phi_o_sigma_n*(g11_n&
          *g33_n-g13_n*g31_n))/den_n
          m_33_up=(g33+delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
          m_33_up_n=(g33_n+delta_0**2*phi_o_sigma_n*(g11_n&
          *g33_n-g13_n*g31_n))/den_n
          den_0=1+delta_0**2*phi_o_sigma*(g22+g44-g24-g42)
          m_22_down=(g22+delta_0**2*phi_o_sigma*(g22*g44-g24*g42))/den_0
          m_24_down=(g24+delta_0**2*phi_o_sigma*(g22*g44-g24*g42))/den_0
          m_42_down=(g42+delta_0**2*phi_o_sigma*(g22*g44-g24*g42))/den_0
          m_44_down=(g44+delta_0**2*phi_o_sigma*(g22*g44-g24*g42))/den_0
          phi_sigma=cdlog((w-B+mu)/(w-A+mu))/(2*D)
          phi_sigma_n=cdlog((w_n-B+mu)/(w_n-A+mu))/(2*D)
          den_2=1-V**2*phi_sigma*(m_11_up+m_13_up+m_31_up+m_33_up)
          den_2_n=1-V**2*phi_sigma_n*(m_11_up_n+m_13_up_n&
                  +m_31_up_n+m_33_up_n)
          !!!! funciones de green del estado localizado
          G_11_up=(m_11_up-V**2*phi_sigma*(m_11_up*m_33_up-&
                  m_13_up*m_31_up))/den_2
          G_11_up_n=(m_11_up_n-V**2*phi_sigma_n*(m_11_up_n*m_33_up_n-&
                  m_13_up_n*m_31_up_n))/den_2_n
          G_13_up=(m_13_up+V**2*phi_sigma*(m_11_up*m_33_up-&
                  m_13_up*m_31_up))/den_2
          G_13_up_n=(m_13_up_n+V**2*phi_sigma_n*(m_11_up_n*m_33_up_n-&
                  m_13_up_n*m_31_up_n))/den_2_n
          G_31_up=(m_31_up+V**2*phi_sigma*(m_11_up*m_33_up-&
                  m_13_up*m_31_up))/den_2
          G_31_up_n=(m_31_up_n+V**2*phi_sigma_n*(m_11_up_n*m_33_up_n-&
                  m_13_up_n*m_31_up_n))/den_2_n
          G_33_up=(m_33_up-V**2*phi_sigma*(m_11_up*m_33_up-&
                  m_13_up*m_31_up))/den_2
          G_33_up_n=(m_33_up_n-V**2*phi_sigma_n*(m_11_up_n*m_33_up_n-&
                  m_13_up_n*m_31_up_n))/den_2_n
          den_3=1-V**2*phi_sigma*(m_22_down+m_44_down-&
                m_24_down-m_42_down)
          G_22_down=(m_22_down-V**2*phi_sigma*(m_22_down*m_44_down-&
                    m_24_down*m_42_down))/den_3
          G_24_down=(m_24_down-V**2*phi_sigma*(m_22_down*m_44_down-&
                    m_24_down*m_42_down))/den_3
          G_42_down=(m_42_down-V**2*phi_sigma*(m_22_down*m_44_down-&
                    m_24_down*m_42_down))/den_3
          G_44_down=(m_44_down-V**2*phi_sigma*(m_22_down*m_44_down-&
                    m_24_down*m_42_down))/den_3
          !!!! funciones de green de la banda de conduccion
          G_cc_up=phi_sigma+(V**2/N_s)*phi_sigma**2*(m_11_up+&
                  m_33_up+m_13_up+m_31_up)/den_2
          G_cf_0=-(V/dsqrt(N_s))*phi_sigma**2*(m_11_up+&
                  m_31_up)/den_2
          G_cf_s=-(V/dsqrt(N_s))*phi_sigma**2*(m_13_up+&
                  m_33_up)/den_2
          G_F=-(G_11_up+G_33_up+G_13_up+G_31_up)
          !!!!!!!! numeros de ocupacion !!!!!!!!!
          G_F_n_2=-(G_22_down+G_44_down-G_24_down-G_42_down)
          sg_n=g11_n+g33_n+g13_n+g31_n
          G_F_G=-1/(1/sg_n+delta_0**2*phi_o_sigma_n-V**2*phi_sigma_n)
          G_00_s=(phi_sigma_n/(1+phi_sigma_n*V**2*&
                 (m_11_up_n+m_33_up_n)))
          S_=(gamma_**2)*(abs(G_00_s)**2)
          !if (j==0 .or. j==n_puntos-1) then
          G_=G_+(-diff_n_F*S_)*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_s=n_f_s+(-1.0d0/pi)*DIMAG(G_F)*n_F*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_cc_up=n_cc_up+(-1.0d0/pi)*DIMAG(-G_cc_up)*n_F*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_00_1=n_f_00_1+(-1.0d0/pi)*DIMAG(G_11_up)*(1.0d0-n_F)*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_pp_1=n_f_pp_1+(-1.0d0/pi)*DIMAG(G_11_up)*n_F*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_mm_1=n_f_mm_1+(-1.0d0/pi)*DIMAG(G_33_up)*(1.0d0-n_F)*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_dd_1=n_f_dd_1+(-1.0d0/pi)*DIMAG(G_33_up)*n_F*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_00_2=n_f_00_2+(-1.0d0/pi)*DIMAG(G_22_down)*(1.0d0-n_F)*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_mm_2=n_f_mm_2+(-1.0d0/pi)*DIMAG(G_22_down)*n_F*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_pp_2=n_f_pp_2+(-1.0d0/pi)*DIMAG(G_44_down)*(1.0d0-n_F)*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          n_f_dd_2=n_f_dd_2+(-1.0d0/pi)*DIMAG(G_44_down)*n_F*&
             abs(abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)-1)
          !else if (mod(j,2)==0 .and. j/=0 .and. j/=n_puntos-1) then
          G_=G_+2.0d0*(-diff_n_F*S_)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_s=n_f_s+2.0d0*(-1.0d0/pi)*DIMAG(G_F)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_cc_up=n_cc_up+2.0d0*(-1.0d0/pi)*DIMAG(-G_cc_up)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_00_1=n_f_00_1+2.0d0*(-1.0d0/pi)*DIMAG(G_11_up)*(1.0d0-n_F)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_pp_1=n_f_pp_1+2.0d0*(-1.0d0/pi)*DIMAG(G_11_up)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_mm_1=n_f_mm_1+2.0d0*(-1.0d0/pi)*DIMAG(G_33_up)*(1.0d0-n_F)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_dd_1=n_f_dd_1+2.0d0*(-1.0d0/pi)*DIMAG(G_33_up)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_00_2=n_f_00_2+2.0d0*(-1.0d0/pi)*DIMAG(G_22_down)*(1.0d0-n_F)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_mm_2=n_f_mm_2+2.0d0*(-1.0d0/pi)*DIMAG(G_22_down)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_pp_2=n_f_pp_2+2.0d0*(-1.0d0/pi)*DIMAG(G_44_down)*(1.0d0-n_F)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_dd_2=n_f_dd_2+2.0d0*(-1.0d0/pi)*DIMAG(G_44_down)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          !else if (mod(j,2)==1 .and. j/=0 .and. j/=n_puntos-1) then
          G_=G_+4.0d0*(-diff_n_F*S_)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_s=n_f_s+4.0d0*(-1.0d0/pi)*DIMAG(G_F)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_cc_up=n_cc_up+4.0d0*(-1.0d0/pi)*DIMAG(-G_cc_up)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_00_1=n_f_00_1+4.0d0*(-1.0d0/pi)*DIMAG(G_11_up)*(1.0d0-n_F)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_pp_1=n_f_pp_1+4.0d0*(-1.0d0/pi)*DIMAG(G_11_up)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_mm_1=n_f_mm_1+4.0d0*(-1.0d0/pi)*DIMAG(G_33_up)*(1.0d0-n_F)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_dd_1=n_f_dd_1+4.0d0*(-1.0d0/pi)*DIMAG(G_33_up)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_00_2=n_f_00_2+4.0d0*(-1.0d0/pi)*DIMAG(G_22_down)*(1.0d0-n_F)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_mm_2=n_f_mm_2+4.0d0*(-1.0d0/pi)*DIMAG(G_22_down)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_pp_2=n_f_pp_2+4.0d0*(-1.0d0/pi)*DIMAG(G_44_down)*(1.0d0-n_F)*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_dd_2=n_f_dd_2+4.0d0*(-1.0d0/pi)*DIMAG(G_44_down)*n_F*&
             abs(int(n_puntos/(j+n_puntos)-1))*&
             abs(int(j/n_puntos)-1)*mod(j,2)
          !end if
          sg=g11+g33+g13+g31
          G_cc_up=-G_cc_up
          G_cf_0=-G_cf_0
          G_cf_s=-G_cf_s
          G_F=-1/(1/sg+delta_0**2*phi_o_sigma-V**2*phi_sigma)
          G_F_n=-(G_11_up+G_33_up+G_13_up+G_31_up)
          write(1,*) DREAL(DCMPLX(w)), (-1/pi)*DIMAG(DCMPLX(G_F_n)),&
          (-1/pi)*DIMAG(DCMPLX(G_cc_up)),(-1/pi)*DIMAG(DCMPLX(G_cf_0)),&
          (-1/pi)*DIMAG(DCMPLX(G_cf_s))
          end do
          sg_mu=g11_mu+g33_mu+g13_mu+g31_mu
          G_F_mu=-1/(1/sg_mu+delta_0**2*phi_o_sigma_mu-V**2*phi_sigma_mu)
          den_mu=1+delta_0**2*phi_o_sigma_mu*&
           (g11_mu+g33_mu+g13_mu+g31_mu)
          m_11_up_mu=(g11_mu+delta_0**2*phi_o_sigma_mu*&
           (g11_mu*g33_mu-g13_mu*g31_mu))/den_mu
          m_13_up_mu=(g13_mu-delta_0**2*phi_o_sigma_mu*&
           (g11_mu*g33_mu-g13_mu*g31_mu))/den_mu
          m_31_up_mu=(g31_mu-delta_0**2*phi_o_sigma_mu*&
           (g11_mu*g33_mu-g13_mu*g31_mu))/den_mu
          m_33_up_mu=(g33_mu+delta_0**2*phi_o_sigma_mu*&
           (g11_mu*g33_mu-g13_mu*g31_mu))/den_mu
          den_2_mu=1-V**2*phi_sigma_mu*&
           (m_11_up_mu+m_13_up_mu+m_31_up_mu+m_33_up_mu)
          G_cc_up_mu=phi_sigma_mu+(V**2/N_s)*phi_sigma_mu**2*&
           (m_11_up_mu+m_33_up_mu+m_13_up_mu+m_31_up_mu)/den_2_mu
          !!!! funcion de green G_F_mu(w=mu)
          G_=G_*(B_I-A_I)/(3.0d0*n_puntos)
          n_f_s=n_f_s*(B_n-A_n)/(3.0d0*n_puntos)
          n_cc_up=n_cc_up*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_00_1=n_f_00_1*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_pp_1=n_f_pp_1*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_mm_1=n_f_mm_1*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_dd_1=n_f_dd_1*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_00_2=n_f_00_2*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_mm_2=n_f_mm_2*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_pp_2=n_f_pp_2*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_dd_2=n_f_dd_2*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_tot_1=n_f_00_1+n_f_pp_1+n_f_mm_1+n_f_dd_1
          n_f_tot_2=n_f_00_2+n_f_pp_2+n_f_mm_2+n_f_dd_2
          pho_fsmu=(-1.0d0/pi)*DIMAG(DCMPLX(G_F_mu))    
          pho_csmu=(-1.0d0/pi)*DIMAG(DCMPLX(-G_cc_up_mu))    
          !F_Eq=(sin(pi*(n_f_s))**2/(delta_0*pi))-pho_fsmu!-abs(E_q)    
          F_Eq=F_Eq_n(E_q)
          print *, "FSR",E_q/delta_0,F_Eq,E_f/delta_0,U/delta_0,G_!,pho_csmu
          if (EU_==1) then
              write(2,*) U/delta_0,2.0*n_f_s,2.0*n_cc_up
              write(3,*) U/delta_0,pho_fsmu,pho_csmu
              write(10,*) U/delta_0,G_
          end if
          if (EU_==0) then
              write(4,*) E_f/delta_0,2.0*n_f_s,2.0*n_cc_up
              write(5,*) E_f/delta_0,pho_fsmu,pho_csmu
              write(7,*) E_f/delta_0,G_
              end if
          close(1)
          end do
          end do
          if (EU_==1) then
              close(2)
              close(3)
              close(10)
          end if
          if (EU_==0) then
              close(4)
              close(5)
              close(7)
          end if
          end
          !!!!!!!!!! funciones !!!!!!!!!
          function f_E(po_,EU_)
              implicit none
              real(8) f_E,E_fA,E_fB,E_f,U,pi,V,x1,x2,delta_0,A,B,mu
              integer po_,E_f_puntos,EU_
              common/VARS1/pi,V,E_f,U,x1,x2,delta_0,A,B,mu
              common/VARS5/E_fA,E_fB,E_f_puntos
              f_E=EU_*E_f+abs(EU_-1)*(E_fA+po_*(E_fB-E_fA)&
              /(E_f_puntos+EU_-1))
              return
              end
          function f_U(po_,EU_)
              implicit none
              real(8) f_U,U_A,U_B,E_f,U,pi,V,x1,x2,delta_0,A,B,mu
              integer po_,U_puntos,EU_
              common/VARS1/pi,V,E_f,U,x1,x2,delta_0,A,B,mu
              common/VARS6/U_A,U_B,U_puntos
              f_U=abs(EU_-1)*U+EU_*(U_A+po_*(U_B-U_A)/&
              (U_puntos+abs(EU_-1)-1))
              return
              end
          function V_E_lim(T) !! valor esperado de la energia
          implicit none
          real(8) pi,delta_0,E_f,E_q,U,D,V,T,T_0,beta,V_E_lim
          real(8) x1,x2,a_2_p,delta,delta_p,Q,R,theta_1   
          real(8) E_1,E_2,E_3,E_4,E_5,E_6,E_7,E_8,E_9,E_10,E_11,E_12
          real(8) E_13,E_14,E_15,E_16,E_min,Z_
          real(8), dimension(16) :: E
          real(8) A,B,A_n,B_n,eta,mu
          complex(8) w_mu
          integer i,j,n_puntos
          common/VARS1/pi,V,E_f,U,x1,x2,delta_0,A,B,mu
          common/VARS2/w_mu,beta,D,A_n,B_n,eta,n_puntos
          common/VARS3/E_q,T_0
          beta=1.0d0/T
          x1=E_f
          x2=E_q
          a_2_p=-(3*x1+3*x2+U)
          delta=dsqrt((x1-x2)**2+4*V**2)
          delta_p=dsqrt((x1+U-x2)**2+4*V**2)
          Q=-(12*V**2+(x1+x2)**2+(2*x1+U)**2+(2*x2)**2-&
          (x2+x1)*(2*x1+U)-(x2+x1)*(2*x2)-(2*x1+U)*(2*x2))/9.0
          R=(-3*((x2+x1)**2*(2*x1+U)+(x2+x1)**2*(2*x2)+&
          (2*x1+U)**2*(x2+x1)+(2*x1+U)**2*(2*x2)+(2*x2)**2*&
          (x2+x1)+(2*x2)**2*(2*x1+U))+12*(x2+x1)*(2*x1+U)*&
          (2*x2)+18*V**2*(2*(x2+x1)-(2*x1+U)-2*x2)+2*&
          ((x2+x1)**3+(2*x1+U)**3+(2*x2)**3))/54.0
          theta_1=acos(R/sqrt(-Q**3))
          E_1=0.0
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
          E(1)=E_1
          E(2)=E_2
          E(3)=E_3
          E(4)=E_4
          E(5)=E_5
          E(6)=E_6
          E(7)=E_7
          E(8)=E_8
          E(9)=E_9
          E(10)=E_10
          E(11)=E_11
          E(12)=E_12
          E(13)=E_13
          E(14)=E_14
          E(15)=E_15
          E(16)=E_16
          E_min=minval(E)
          Z_=0.0d0
          do i=1,16
              Z_=Z_+dexp(-beta*(E(i)-E_min))
          end do
          V_E_lim=0.0d0
          do j=1,16
              V_E_lim=V_E_lim+E(j)*dexp(-beta*(E(j)-E_min))
          end do
          V_E_lim=V_E_lim/Z_
          T=T_0
          beta=1/T_0
          return
          end
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          function C_v_lim(T) !! capacidad calorifica a volumen constante
          implicit none
          real(8) pi,delta_0,E_f,E_q,U,D,V,T,T_0,beta,C_v_lim
          real(8) x1,x2,a_2_p,delta,delta_p,Q,R,theta_1,E2_prom,E_prom
          real(8) E_1,E_2,E_3,E_4,E_5,E_6,E_7,E_8,E_9,E_10,E_11,E_12
          real(8) E_13,E_14,E_15,E_16,E_min,Z_
          real(8), dimension(16) :: E
          real(8) A,B,A_n,B_n,eta,mu
          complex(8) w_mu
          integer i,j,n_puntos
          common/VARS1/pi,V,E_f,U,x1,x2,delta_0,A,B,mu
          common/VARS2/w_mu,beta,D,A_n,B_n,eta,n_puntos
          common/VARS3/E_q,T_0
          beta=1.0d0/T
          x1=E_f
          x2=E_q
          a_2_p=-(3*x1+3*x2+U)
          delta=dsqrt((x1-x2)**2+4*V**2)
          delta_p=dsqrt((x1+U-x2)**2+4*V**2)
          Q=-(12*V**2+(x1+x2)**2+(2*x1+U)**2+(2*x2)**2-&
          (x2+x1)*(2*x1+U)-(x2+x1)*(2*x2)-(2*x1+U)*(2*x2))/9.0
          R=(-3*((x2+x1)**2*(2*x1+U)+(x2+x1)**2*(2*x2)+&
          (2*x1+U)**2*(x2+x1)+(2*x1+U)**2*(2*x2)+(2*x2)**2*&
          (x2+x1)+(2*x2)**2*(2*x1+U))+12*(x2+x1)*(2*x1+U)*&
          (2*x2)+18*V**2*(2*(x2+x1)-(2*x1+U)-2*x2)+2*&
          ((x2+x1)**3+(2*x1+U)**3+(2*x2)**3))/54.0
          theta_1=acos(R/sqrt(-Q**3))
          E_1=0.0
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
          E(1)=E_1
          E(2)=E_2
          E(3)=E_3
          E(4)=E_4
          E(5)=E_5
          E(6)=E_6
          E(7)=E_7
          E(8)=E_8
          E(9)=E_9
          E(10)=E_10
          E(11)=E_11
          E(12)=E_12
          E(13)=E_13
          E(14)=E_14
          E(15)=E_15
          E(16)=E_16
          E_min=minval(E)
          Z_=0.0d0
          E_prom=0.0d0
          E2_prom=0.0d0
          do i=1,16
              Z_=Z_+dexp(-beta*(E(i)-E_min))
              E_prom=E_prom+E(i)*dexp(-beta*(E(i)-E_min))
              E2_prom=E2_prom+(E(i)**2)*dexp(-beta*(E(i)-E_min))
          end do
          E_prom=E_prom/Z_
          E2_prom=E2_prom/Z_
          !do j=1,16
          !    C_v=C_v+E(j)*dexp(-beta*(E(j)-E_min))
          !end do
          !C_v=-C_v**2/(Z_**2)
          !do j=1,16
          !    C_v=C_v+E(j)**2*dexp(-beta*(E(j)-E_min))/Z_
          !end do
          !C_v=C_v/T**2
          C_v_lim=(E2_prom-E_prom**2)/T**2
          T=T_0
          beta=1/T_0
          return
          end
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          function V_E(T) !! capacidad calorifica a volumen constante,
          !!! en el caso limite y el valor esperado de la energia
          !!! usando los numeros de ocupacion
          implicit none
          real(8) pi,delta_0,E_f,E_q,U,D,V,T,T_0,beta,C_v,V_E
          real(8) x1,x2,a_2_p,delta,delta_p,Q,R,theta_1,E2_prom,E_prom
          real(8) E_1,E_2,E_3,E_4,E_5,E_6,E_7,E_8,E_9,E_10,E_11,E_12
          real(8) E_13,E_14,E_15,E_16,E_min,Z_,N_s
          real(8) phi,theta,a_9,a_10,a_11,b_9,b_10,b_11,c_9,c_10,c_11
          real(8) a_0,b_0,c_0
          real(8), dimension(3) :: a_,b_,c_
          complex(8) g11,g33,g13,g31,w
          complex(8) m_11_up,m_13_up,m_31_up,m_33_up
          complex(8) phi_sigma,phi_o_sigma,den_,den_0,den_2
          complex(8) G_cc_up,G_cf_0,G_cf_s
          complex(8) G_11_up,G_13_up,G_31_up,G_33_up
          real(8) n_cc_up,n_cf_0,n_cf_s,n_F
          real(8) n_f_00_1,n_f_pp_1,n_f_mm_1,n_f_dd_1
          real(8), dimension(16) :: E
          real(8) A,B,A_n,B_n,eta,mu
          complex(8) w_mu
          integer i,j,n_puntos
          common/VARS1/pi,V,E_f,U,x1,x2,delta_0,A,B,mu
          common/VARS2/w_mu,beta,D,A_n,B_n,eta,n_puntos
          common/VARS3/E_q,T_0
          common/VARS4/N_s
          external a_0,b_0,c_0
          beta=1.0d0/T
          x1=E_f
          x2=E_q
          V=delta_0
          a_2_p=-(3*x1+3*x2+U)
          delta=dsqrt((x1-x2)**2+4*V**2)
          delta_p=dsqrt((x1+U-x2)**2+4*V**2)
          Q=-(12*V**2+(x1+x2)**2+(2*x1+U)**2+(2*x2)**2-&
          (x2+x1)*(2*x1+U)-(x2+x1)*(2*x2)-(2*x1+U)*(2*x2))/9.0
          R=(-3*((x2+x1)**2*(2*x1+U)+(x2+x1)**2*(2*x2)+&
          (2*x1+U)**2*(x2+x1)+(2*x1+U)**2*(2*x2)+(2*x2)**2*&
          (x2+x1)+(2*x2)**2*(2*x1+U))+12*(x2+x1)*(2*x1+U)*&
          (2*x2)+18*V**2*(2*(x2+x1)-(2*x1+U)-2*x2)+2*&
          ((x2+x1)**3+(2*x1+U)**3+(2*x2)**3))/54.0
          theta_1=acos(R/sqrt(-Q**3))
          E_1=0.0
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
          E(1)=E_1
          E(2)=E_2
          E(3)=E_3
          E(4)=E_4
          E(5)=E_5
          E(6)=E_6
          E(7)=E_7
          E(8)=E_8
          E(9)=E_9
          E(10)=E_10
          E(11)=E_11
          E(12)=E_12
          E(13)=E_13
          E(14)=E_14
          E(15)=E_15
          E(16)=E_16
          E_min=minval(E)
          phi=atan(2.0*V/(E_q-E_f+delta))
          theta=atan(2.0*V/(E_f+U-E_q-delta_p))
          a_9=a_0(E_9,E_f,E_q,U,V)
          a_10=a_0(E_10,E_f,E_q,U,V)
          a_11=a_0(E_11,E_f,E_q,U,V)
          a_(1)=a_9
          a_(2)=a_10
          a_(3)=a_11
          b_9=b_0(E_9,E_f,E_q,U,V)
          b_10=b_0(E_10,E_f,E_q,U,V)
          b_11=b_0(E_11,E_f,E_q,U,V)
          b_(1)=b_9
          b_(2)=b_10
          b_(3)=b_11
          c_9=c_0(E_9,E_f,E_q,U,V,E_9)
          c_10=c_0(E_10,E_f,E_q,U,V,E_10)
          c_11=c_0(E_11,E_f,E_q,U,V,E_11)
          c_(1)=c_9
          c_(2)=c_10
          c_(3)=c_11
          Z_=0.0d0
!          E_prom=0.0d0
!          E2_prom=0.0d0
          !!!!!!!! Numeros de ocupacion !!!!!!!!
          n_cc_up=0.0d0
          n_cf_0=0.0d0
          n_cf_s=0.0d0
          n_f_00_1=0.0d0
          n_f_pp_1=0.0d0
          n_f_mm_1=0.0d0
          n_f_dd_1=0.0d0
         do i=1,16
             Z_=Z_+dexp(-beta*(E(i)-E_min))
!             E_prom=E_prom+E(i)*dexp(-beta*(E(i)-E_min))
!             E2_prom=E2_prom+(E(i)**2)*dexp(-beta*(E(i)-E_min))
         end do
!         E_prom=E_prom/Z_
!         E2_prom=E2_prom/Z_
!         C_v=(E2_prom-E_prom**2)/T**2 !! capacidad calorifica en el
          !!caso limite
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
          do i=1,16
              E(i)=E(i)-E_min
          end do
          do j=0,n_puntos-1
          w=DCMPLX((A_n+j*(B_n-A_n)/(n_puntos-1)),eta)
          n_F=1.0d0/(1+exp(beta*DREAL(w)))
          !!!! g11
          g11=(sin(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_4))&
          /(w+E_1-E_4)+&
          1.5*(dexp(-beta*E_2)+dexp(-beta*E_6))/(w+E_2-E_6))+&
          cos(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_2))/(w+E_1-E_2)+&
          1.5*(dexp(-beta*E_4)+dexp(-beta*E_6))/(w+E_4-E_6)))    
          do i=9,11
          g11=g11+(((dexp(-beta*E_3)+dexp(-beta*E(i)))/(w+E_3-E(i)))*&
          (a_(i-8)*sin(phi))**2+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w+E_5-E(i)))*(a_(i-8)*cos(phi))**2+((dexp(-beta*E(i))+&
          dexp(-beta*E_14))/(w+E(i)-E_14))*(c_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w+E(i)-E_12))*(c_(i-8)*&
          cos(theta))**2)
          end do
          g11=g11/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g33
          g33=(sin(theta)**2*((dexp(-beta*E_15)+dexp(-beta*E_16))/&
          (w+E_15-E_16)+1.5*(dexp(-beta*E_8)+dexp(-beta*E_12))/&
          (w+E_8-E_12))+cos(theta)**2*((dexp(-beta*E_13)+&
          dexp(-beta*E_16))/(w+E_13-E_16)+1.5*(dexp(-beta*E_8)+&
          dexp(-beta*E_14))/(w+E_8-E_14)))
          do i=9,11
          g33=g33+((dexp(-beta*E_4)+dexp(-beta*E(i)))/(w+E_4-E(i)))*&
          (b_(i-8)*sin(phi))**2+((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w+E_3-E(i)))*(b_(i-8)*cos(phi))**2+((dexp(-beta*E(i))&
          +dexp(-beta*E_12))/(w+E(i)-E_12))*(a_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w+E(i)-E_14))*&
          (a_(i-8)*cos(theta))**2
          end do
          g33=g33/Z_
          !!! g13
          g13=0.0d0
          do i=9,11
          g13=g13+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w+E_5-E(i))-(dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w+E_3-E(i)))*(a_(i-8)*b_(i-8)*sin(phi)*cos(phi))+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w+E(i)-E_14)-&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w+E(i)-E_12)))*&
          (a_(i-8)*c_(i-8)*sin(theta)*cos(theta))
          end do
          g13=g13/Z_
          g31=g13
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          g11=-g11
          g33=-g33
          g13=-g13
          g31=-g31
          V=dsqrt(delta_0*2.0d0*D/pi)  ! Hibridacion  
          phi_o_sigma=-1/(w-E_q-mu)
          den_=1+delta_0**2*phi_o_sigma*(g11+g33+g13+g31)
          m_11_up=(g11+delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
          m_13_up=(g13-delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
          m_31_up=(g31-delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
          m_33_up=(g33+delta_0**2*phi_o_sigma*(g11*g33-g13*g31))/den_
          phi_sigma=cdlog((w-B+mu)/(w-A+mu))/(2*D)
          den_2=1-V**2*phi_sigma*(m_11_up+m_13_up+m_31_up+m_33_up)
          !!! funciones de green de la banda de conduccion y el esto
          !!! cruzado
          G_cc_up=phi_sigma+(V**2/N_s)*phi_sigma**2*(m_11_up+&
                  m_33_up+m_13_up+m_31_up)/den_2
          G_cc_up=-G_cc_up
          G_cf_0=-(V/dsqrt(N_s))*phi_sigma**2*(m_11_up+&
                  m_31_up)/den_2
          G_cf_s=-(V/dsqrt(N_s))*phi_sigma**2*(m_13_up+&
                  m_33_up)/den_2
          !!! funciones de green del estado localizado
          G_11_up=(m_11_up-V**2*phi_sigma*(m_11_up*m_33_up-&
                  m_13_up*m_31_up))/den_2
          G_13_up=(m_13_up+V**2*phi_sigma*(m_11_up*m_33_up-&
                  m_13_up*m_31_up))/den_2
          G_31_up=(m_31_up+V**2*phi_sigma*(m_11_up*m_33_up-&
                  m_13_up*m_31_up))/den_2
          G_33_up=(m_33_up-V**2*phi_sigma*(m_11_up*m_33_up-&
                  m_13_up*m_31_up))/den_2
          !!!!!!!! numeros de ocupacion !!!!!!!!!
          !if (j==0 .or. j==n_puntos-1) then
          n_cc_up=n_cc_up+(-1.0d0/pi)*DIMAG(G_cc_up)*(n_F)*&
            abs(abs(int(n_puntos/(j+n_puntos))-1)*&
            abs(int(j/n_puntos)-1)-1)
          n_cf_0=n_cf_0+(-1.0d0/pi)*DIMAG(G_cf_0)*n_F*&
            abs(abs(int(n_puntos/(j+n_puntos))-1)*&
            abs(int(j/n_puntos)-1)-1)
          n_cf_s=n_cf_s+(-1.0d0/pi)*DIMAG(G_cf_s)*&
            abs(abs(int(n_puntos/(j+n_puntos))-1)*&
            abs(int(j/n_puntos)-1)-1)
          n_f_00_1=n_f_00_1+(-1.0d0/pi)*DIMAG(G_11_up)*(1.0d0-n_F)*&
            abs(abs(int(n_puntos/(j+n_puntos))-1)*&
            abs(int(j/n_puntos)-1)-1)
          n_f_pp_1=n_f_pp_1+(-1.0d0/pi)*DIMAG(G_11_up)*n_F*&
            abs(abs(int(n_puntos/(j+n_puntos))-1)*&
            abs(int(j/n_puntos)-1)-1)
          n_f_mm_1=n_f_mm_1+(-1.0d0/pi)*DIMAG(G_33_up)*(1.0d0-n_F)*&
            abs(abs(int(n_puntos/(j+n_puntos))-1)*&
            abs(int(j/n_puntos)-1)-1)
          n_f_dd_1=n_f_dd_1+(-1.0d0/pi)*DIMAG(G_33_up)*n_F*&
            abs(abs(int(n_puntos/(j+n_puntos))-1)*&
            abs(int(j/n_puntos)-1)-1)
          !else if (mod(j,2)==0 .and. j/=0 .and. j/=n_puntos-1) then
          n_cc_up=n_cc_up+2.0d0*(-1.0d0/pi)*DIMAG(G_cc_up)*(n_F)*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_cf_0=n_cf_0+2.0d0*(-1.0d0/pi)*DIMAG(G_cf_0)*n_F*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_cf_s=n_cf_s+2.0d0*(-1.0d0/pi)*DIMAG(G_cf_s)*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_00_1=n_f_00_1+2.0d0*(-1.0d0/pi)*DIMAG(G_11_up)*(1.0d0-n_F)*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_pp_1=n_f_pp_1+2.0d0*(-1.0d0/pi)*DIMAG(G_11_up)*n_F*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_mm_1=n_f_mm_1+2.0d0*(-1.0d0/pi)*DIMAG(G_33_up)*(1.0d0-n_F)*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          n_f_dd_1=n_f_dd_1+2.0d0*(-1.0d0/pi)*DIMAG(G_33_up)*n_F*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          !else if (mod(j,2)==1 .and. j/=0 .and. j/=n_puntos-1) then
          n_cc_up=n_cc_up+4.0d0*(-1.0d0/pi)*DIMAG(G_cc_up)*(n_F)*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*mod(j,2)
          n_cf_0=n_cf_0+4.0d0*(-1.0d0/pi)*DIMAG(G_cf_0)*n_F*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*mod(j,2)
          n_cf_s=n_cf_s+4.0d0*(-1.0d0/pi)*DIMAG(G_cf_s)*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_00_1=n_f_00_1+4.0d0*(-1.0d0/pi)*DIMAG(G_11_up)*(1.0d0-n_F)*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_pp_1=n_f_pp_1+4.0d0*(-1.0d0/pi)*DIMAG(G_11_up)*n_F*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_mm_1=n_f_mm_1+4.0d0*(-1.0d0/pi)*DIMAG(G_33_up)*(1.0d0-n_F)*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*mod(j,2)
          n_f_dd_1=n_f_dd_1+4.0d0*(-1.0d0/pi)*DIMAG(G_33_up)*n_F*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*mod(j,2)
          !end if
          end do
          !!!! funcion de green G_F_mu(w=mu)
          n_cc_up=n_cc_up*(B_n-A_n)/(3.0d0*n_puntos)
          n_cf_0=n_cf_0*(B_n-A_n)/(3.0d0*n_puntos)
          n_cf_s=n_cf_s*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_00_1=-n_f_00_1*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_pp_1=-n_f_pp_1*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_mm_1=-n_f_mm_1*(B_n-A_n)/(3.0d0*n_puntos)
          n_f_dd_1=-n_f_dd_1*(B_n-A_n)/(3.0d0*n_puntos)
          V_E=E_q*(2.0d0*n_cc_up)+E_f*(n_f_00_1+n_f_pp_1+n_f_mm_1)+&
              (E_f+U)*n_f_dd_1+2*V*(n_cf_0+n_cf_s)
          T=T_0
          beta=1/T_0
          return
          end
          function F_Eq_n(E_q)
          implicit none
          real(8) pi,eta,delta_0,E_f,E_q,U,A,B,A_n,B_n,D,V,T,beta,mu,n_F
          real(8) x1,x2,a_2_p,delta,delta_p,Q,R,theta_1   
          real(8) E_1,E_2,E_3,E_4,E_5,E_6,E_7,E_8,E_9,E_10,E_11,E_12
          real(8) E_13,E_14,E_15,E_16,E_min,Z_
          real(8) phi,theta,a_9,a_10,a_11,b_9,b_10,b_11,c_9,c_10,c_11
          real(8) a_0,b_0,c_0
          complex(8) g11,g33,g13,g31,g11_mu,g33_mu,g13_mu,g31_mu,&
          phi_sigma,phi_o_sigma,phi_sigma_mu,phi_o_sigma_mu,w,w_mu,&
          G_F,G_F_mu,sg,sg_mu
          real(8) F_Eq_n,n_f_s,pho_fsmu
          real(8), dimension(16) :: E
          real(8), dimension(3) :: a_,b_,c_
          integer i,j,k,l,n_puntos
          common/VARS1/pi,V,E_f,U,x1,x2,delta_0,A,B,mu
          common/VARS2/w_mu,beta,D,A_n,B_n,eta,n_puntos
          common/T_/T
          external a_0,b_0,c_0
          !!!!!!!!!!!!
          !!!!!!!!!!!!
          x1=E_f
          x2=E_q
          a_2_p=-(3*x1+3*x2+U)
          V=delta_0
          delta=dsqrt((x1-x2)**2+4*V**2)
          delta_p=dsqrt((x1+U-x2)**2+4*V**2)
          Q=-(12*V**2+(x1+x2)**2+(2*x1+U)**2+(2*x2)**2-&
          (x2+x1)*(2*x1+U)-(x2+x1)*(2*x2)-(2*x1+U)*(2*x2))/9.0
          R=(-3*((x2+x1)**2*(2*x1+U)+(x2+x1)**2*(2*x2)+&
          (2*x1+U)**2*(x2+x1)+(2*x1+U)**2*(2*x2)+(2*x2)**2*&
          (x2+x1)+(2*x2)**2*(2*x1+U))+12*(x2+x1)*(2*x1+U)*&
          (2*x2)+18*V**2*(2*(x2+x1)-(2*x1+U)-2*x2)+2*&
          ((x2+x1)**3+(2*x1+U)**3+(2*x2)**3))/54.0
          theta_1=acos(R/sqrt(-Q**3))
          E_1=0.0
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
          E(1)=E_1
          E(2)=E_2
          E(3)=E_3
          E(4)=E_4
          E(5)=E_5
          E(6)=E_6
          E(7)=E_7
          E(8)=E_8
          E(9)=E_9
          E(10)=E_10
          E(11)=E_11
          E(12)=E_12
          E(13)=E_13
          E(14)=E_14
          E(15)=E_15
          E(16)=E_16
          E_min=minval(E)
          phi=atan(2.0*V/(E_q-E_f+delta))
          theta=atan(2.0*V/(E_f+U-E_q-delta_p))
          a_9=a_0(E_9,E_f,E_q,U,V)
          a_10=a_0(E_10,E_f,E_q,U,V)
          a_11=a_0(E_11,E_f,E_q,U,V)
          a_(1)=a_9
          a_(2)=a_10
          a_(3)=a_11
          b_9=b_0(E_9,E_f,E_q,U,V)
          b_10=b_0(E_10,E_f,E_q,U,V)
          b_11=b_0(E_11,E_f,E_q,U,V)
          b_(1)=b_9
          b_(2)=b_10
          b_(3)=b_11
          c_9=c_0(E_9,E_f,E_q,U,V,E_9)
          c_10=c_0(E_10,E_f,E_q,U,V,E_10)
          c_11=c_0(E_11,E_f,E_q,U,V,E_11)
          c_(1)=c_9
          c_(2)=c_10
          c_(3)=c_11
          !!!!!!!! Funcion de particion !!!!
          Z_=0.0d0
          do i=1,16
              Z_=Z_+dexp(-beta*(E(i)-E_min))
          end do
          !!!!!!!! Numeros de ocupacion !!!!!!!!
          n_f_s=0.0d0
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!! Conductancia !!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
          do i=1,16
              E(i)=E(i)-E_min
          end do
          !!!! g11_mu g11(w=mu)
          g11_mu=(sin(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_4))&
          /(w_mu+E_1-E_4)+&
          1.5*(dexp(-beta*E_2)+dexp(-beta*E_6))/(w_mu+E_2-E_6))+&
          cos(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_2))/(w_mu+E_1-E_2)&
          +1.5*(dexp(-beta*E_4)+dexp(-beta*E_6))/(w_mu+E_4-E_6)))    
          do i=9,11
          g11_mu=g11_mu+(((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_mu+E_3-E(i)))*&
          (a_(i-8)*sin(phi))**2+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w_mu+E_5-E(i)))*(a_(i-8)*cos(phi))**2+((dexp(-beta*E(i))+&
          dexp(-beta*E_14))/(w_mu+E(i)-E_14))*(c_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w_mu+E(i)-E_12))*&
          (c_(i-8)*cos(theta))**2)
          end do
          g11_mu=g11_mu/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          phi_sigma_mu=cdlog((w_mu-B+mu)/(w_mu-A+mu))/(2*D)
          phi_o_sigma_mu=-1/(w_mu-E_q-mu)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g33_mu g33(w=mu)
          g33_mu=(sin(theta)**2*((dexp(-beta*E_15)+dexp(-beta*E_16))/&
          (w_mu+E_15-E_16)+1.5*(dexp(-beta*E_8)+dexp(-beta*E_12))/&
          (w_mu+E_8-E_12))+cos(theta)**2*((dexp(-beta*E_13)+&
          dexp(-beta*E_16))/(w_mu+E_13-E_16)+1.5*(dexp(-beta*E_8)+&
          dexp(-beta*E_14))/(w_mu+E_8-E_14)))
          do i=9,11
          g33_mu=g33_mu+((dexp(-beta*E_4)+dexp(-beta*E(i)))/&
          (w_mu+E_4-E(i)))*&
          (b_(i-8)*sin(phi))**2+((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_mu+E_3-E(i)))*(b_(i-8)*cos(phi))**2+((dexp(-beta*E(i))&
          +dexp(-beta*E_12))/(w_mu+E(i)-E_12))*(a_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w_mu+E(i)-E_14))*&
          (a_(i-8)*cos(theta))**2
          end do
          g33_mu=g33_mu/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g13
          g13_mu=0.0d0
          do i=9,11
          g13_mu=g13_mu+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w_mu+E_5-E(i))-(dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w_mu+E_3-E(i)))*(a_(i-8)*b_(i-8)*sin(phi)*cos(phi))+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w_mu+E(i)-E_14)-&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w_mu+E(i)-E_12)))*&
          (a_(i-8)*c_(i-8)*sin(theta)*cos(theta))
          end do
          g13_mu=g13_mu/Z_
          g31_mu=g13_mu
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          g11_mu=-g11_mu
          g33_mu=-g33_mu
          g13_mu=-g13_mu
          g31_mu=-g31_mu
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!! bucle para el calculo de la densidad de estados !!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do j=0,n_puntos-1
          w=DCMPLX((A_n+j*(B_n-A_n)/(n_puntos-1)),eta)    
          n_F=1.0d0/(1+exp(beta*DREAL(w)))
          !!!! g11
          g11=(sin(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_4))&
          /(w+E_1-E_4)+&
          1.5*(dexp(-beta*E_2)+dexp(-beta*E_6))/(w+E_2-E_6))+&
          cos(phi)**2*((dexp(-beta*E_1)+dexp(-beta*E_2))/(w+E_1-E_2)+&
          1.5*(dexp(-beta*E_4)+dexp(-beta*E_6))/(w+E_4-E_6)))    
          do i=9,11
          g11=g11+(((dexp(-beta*E_3)+dexp(-beta*E(i)))/(w+E_3-E(i)))*&
          (a_(i-8)*sin(phi))**2+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w+E_5-E(i)))*(a_(i-8)*cos(phi))**2+((dexp(-beta*E(i))+&
          dexp(-beta*E_14))/(w+E(i)-E_14))*(c_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w+E(i)-E_12))*(c_(i-8)*&
          cos(theta))**2)
          end do
          g11=g11/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g33
          g33=(sin(theta)**2*((dexp(-beta*E_15)+dexp(-beta*E_16))/&
          (w+E_15-E_16)+1.5*(dexp(-beta*E_8)+dexp(-beta*E_12))/&
          (w+E_8-E_12))+cos(theta)**2*((dexp(-beta*E_13)+&
          dexp(-beta*E_16))/(w+E_13-E_16)+1.5*(dexp(-beta*E_8)+&
          dexp(-beta*E_14))/(w+E_8-E_14)))
          do i=9,11
          g33=g33+((dexp(-beta*E_4)+dexp(-beta*E(i)))/(w+E_4-E(i)))*&
          (b_(i-8)*sin(phi))**2+((dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w+E_3-E(i)))*(b_(i-8)*cos(phi))**2+((dexp(-beta*E(i))&
          +dexp(-beta*E_12))/(w+E(i)-E_12))*(a_(i-8)*sin(theta))**2+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w+E(i)-E_14))*&
          (a_(i-8)*cos(theta))**2
          end do
          g33=g33/Z_
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!! g13
          g13=0.0d0
          do i=9,11
          g13=g13+((dexp(-beta*E_5)+dexp(-beta*E(i)))/&
          (w+E_5-E(i))-(dexp(-beta*E_3)+dexp(-beta*E(i)))/&
          (w+E_3-E(i)))*(a_(i-8)*b_(i-8)*sin(phi)*cos(phi))+&
          ((dexp(-beta*E(i))+dexp(-beta*E_14))/(w+E(i)-E_14)-&
          ((dexp(-beta*E(i))+dexp(-beta*E_12))/(w+E(i)-E_12)))*&
          (a_(i-8)*c_(i-8)*sin(theta)*cos(theta))
          end do
          g13=g13/Z_
          g31=g13
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          g11=-g11
          g33=-g33
          g13=-g13
          g31=-g31
          V=dsqrt(delta_0*2.0d0*D/pi)  ! Hibridacion  
          phi_o_sigma=-1/(w-E_q-mu)
          phi_sigma=cdlog((w-B+mu)/(w-A+mu))/(2*D)
          sg=g11+g33+g13+g31
          G_F=-1/(1/sg+delta_0**2*phi_o_sigma-V**2*phi_sigma)
          !!!!!!!! numeros de ocupacion !!!!!!!!!!
          !if (j==0 .or. j==n_puntos-1) then
          n_f_s=n_f_s+(-1.0d0/pi)*DIMAG(G_F)*n_F*&
                abs(abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)-1)
          !else if (mod(j,2)==0 .and. j/=0 .and. j/=n_puntos-1) then
          n_f_s=n_f_s+2.0d0*(-1.0d0/pi)*DIMAG(G_F)*n_F*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*abs(mod(j,2)-1)
          !else if (mod(j,2)==1 .and. j/=0 .and. j/=n_puntos-1) then
          n_f_s=n_f_s+4.0d0*(-1.0d0/pi)*DIMAG(G_F)*n_F*&
                abs(int(n_puntos/(j+n_puntos))-1)*&
                abs(int(j/n_puntos)-1)*mod(j,2)
          !end if
          end do
          sg_mu=g11_mu+g33_mu+g13_mu+g31_mu
          G_F_mu=-1/(1/sg_mu+delta_0**2*phi_o_sigma_mu-V**2*phi_sigma_mu)
          !!!! funcion de green G_F_mu(w=mu)
          n_f_s=n_f_s*(B_n-A_n)/(3.0d0*n_puntos)
          pho_fsmu=(-1.0d0/pi)*DIMAG(DCMPLX(G_F_mu))    
          F_Eq_n=abs((sin(pi*(n_f_s))**2/(delta_0*pi))-pho_fsmu)!-abs(E_q))
          !!! descomentar para evitar el ruido y discontinuidades
          !!! para calcular conductancia dejarlo comentado
          return
          end
          function a_0(E_i,E_f,E_q,U,V)
          implicit none
          real(8) a_0,E_i,E_f,E_q,U,V
          a_0=1.0/sqrt(2.0+4.0*V**2*(1.0/ &
          (E_i-2*E_f-U)**2+1.0/(E_i-2*E_q)**2))
          return
          end
          !!!!!!!!!!!
          function b_0(E_i,E_f,E_q,U,V)
          implicit none
          real(8) a,b_0,E_i,E_f,E_q,U,V
          a=1.0/sqrt(2.0+4.0*V**2*(1.0/ &
          (E_i-2*E_f-U)**2+1.0/(E_i-2*E_q)**2))
          b_0=2.0*V*a/(E_i-2*E_f-U)
          return 
          end
          !!!!!!!!!!!!
          function c_0(E_i,E_f,E_q,U,V,E_9)
          implicit none
          real(8) a,c_0,E_i,E_f,E_q,U,V,E_9
          a=1.0/sqrt(2.0+4.0*V**2*(1.0/ &
          (E_i-2*E_f-U)**2+1.0/(E_i-2*E_q)**2))
          c_0=2.0*V*a/(E_9-2*E_q)
          return
          end
          !!!!!!!!!!!!! convert integer to string
          character(len=20) function str(k) 
          ! "Convert an integer to string." 
              integer, intent(in) :: k
              write (str, *) k
              str = adjustl(str)
              return
          end
!*********************************************************************************************
         Function Fmin_(ax,bx,f,Tol)
!      Para funcion f(x) unimodular dentro del intervalo (ax,bx)
!      se encuentra el punto x=Fmin, con incertidumbre Tol,
!      donde la funcion acepta el valor minimo
!
!  Parametros de la entrada:
!        ax, bx -  los extremos del intervalo analizado;
!        f      -  subprograma-funcion, que calcula el valor
!                 de la funcion f(x) dentro (ax,bx)
!        Tol    - la precision requerida
!  Parametros de la salida:
!        Fmin   - la abscisa, que aproxima el punto, donde
!                 f(x) tiene el minimo
!  Metodo:
!     Se utiliza la combinacion del metodo de Fibonacci y
!     de interpolacion cuadratica sucesiva. El algoritmo
!     no puede ser menos rapido que el metodo de Fibonacci
!     y es mucho mas rapido si la funcion f(x) tiene la segunda
!     derivada continua.
!     Programa utiliza el minimo posible accesos al calculo
!     de la funcion, en cado paso se calcula el nuevo valor
!     no mas que una vez y nunca se calcula el nuevo valor,
!     si la distancia entre dos puntos sucesivos es menor que
!     macheps*abs(x)+tol/3, donde macheps es la precision
!     aritmetica del computador y x la abscisa del punto anterior.
!     Se utilizo el algoritmo, que esta descrito en el libro:
!     RICHARD BRENT. "ALGORITHMS FOR MINIMIZATION
!     WITHOUT DERIVATVES", Prentice -Hall,1973
            Implicit double precision (a-h,o-z)
            External f
!    El parametro C es el cuadrato del valor inversa
!            a seccion de oro
            C=0.5d0*(3.0d0-dsqrt(5.0d0))
!   El parametro Eps es igual a la raiz cuadrada de
!   la precision aritmetica del computador
            Eps=1.0d0
10          Eps=Eps/2.0d0
            Tol1=1.0d0+Eps
            If(tol1.gt.1.0d0) go to 10
            Eps=dsqrt(Eps)
!  Asignacion de los valores iniciales
            A=ax
            B=bx
            V=A+C*(B-A)
            W=V
            X=V
            E=0.0d0
            FX=f(X)
            FV=FX
            FW=FX
!  Aqui se inicia el ciclo principal
20          XM=0.5d0*(A+B)
            TOL1=Eps*Dabs(X)+tol/3.0d0
            TOL2=2.0d0*TOL1
!  Chequeo de terminacion del ciclo principal
            If(Dabs(X-XM).le.(TOL2-0.5d0*(B-A))) go to 90
!   Es necesario la seccion de oro ?
            If(Abs(E).le.Tol1) go to 40
!   Construir la parabola
            R=(X-W)*(FX-FV)
            Q=(X-V)*(FX-FW)
            P=(X-V)*Q-(X-W)*R
            Q=2.0d0*(Q-R)
            If(Q.gt.0.0d0) P=-P
            Q=Dabs(Q)
            R=E
            E=D
!     Parabola es aceptable ?
30          If(dabs(P).ge.dabs((0.5d0*Q*R))) go to 40
            If(P.ge.Q*(A-X)) go to 40
            If(P.ge.Q*(B-X)) go to 40
!   El paso de la interpolacion parabolica
            D=P/Q
            U=X+D
!   F no se permite calcular demasiado cerca al AX o BX
            If((U-A).lt.Tol2) D=Dsign(Tol1,XM-X)
            If((B-U).lt.Tol2) D=Dsign(Tol1,XM-X)
            Go to 50
!   El paso de seccion de oro
40          If(X.ge.XM) E=A-X
            If(X.lt.XM) E=B-X
            D=C*E
!   F no se permite calcular demasiado cerca al X
50          If(Dabs(D).ge.Tol1) U=X+D
            If(Dabs(D).lt.Tol1) U=X+Dsign(Tol1,D)
            FU=f(u)
!  Asignar valores nuevos a los parametros A,B,V,W y X
            If(FU.gt.FX) go to 60
            If(U.ge.X) A=X
            If(U.lt.X) B=X
            V=W
            FV=FW
            W=X
            FW=FX
            X=U
            FX=FU
            go to 20
60          if(U.lt.X) A=U
            If(U.ge.X) B=U
            If(FU.le.FW) go to 70
            If(W.eq.X) go to 70
            If(FU.le.FV) go to 80
            If(V.eq.X) go to 80
            If(V.eq.W) go to 80
            go to 20
70          V=W
            FV=FW
            W=U
            FW=FU
            go to 20
80          V=U
            FV=FU
            go to 20
!   Fin del ciclo principal
90          Fmin=X
            Return
            End
       
! ************************************************************************************************
