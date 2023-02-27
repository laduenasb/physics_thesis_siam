from sympy import *
x1,x2,U,V,lam=symbols("xi_f xi_q U V lambda")
f=(2*x1+U-lam)*(x1+x2-lam)*(2*x2-lam)-2*V**2*(2*(x1+x2)+U-2*lam)
f_c=collect(expand(-f),lam)
print "Polinomio de grado 3 a solucionar por formula de cardano"
pprint(f_c)
a2=-(3*x1+3*x2+U)
a1=x1*U+3*x2*U-4*V**2+2*(x1+x2)**2+4*x1*x2
a0=2*U*V**2-2*U*x1*x2-2*U*x2**2+4*V**2*x1+4*V**2*x2-4*x1**2*x2-4*x1*x2**2
Qp=-(3*a1-a2**2)
Qp_articulo=12*V**2+(x1+x2)**2+(2*x1+U)**2+(2*x2)**2-(x2+x1)*(2*x1+U)-(x2+x1)*(2*x2)-(2*x1+U)*(2*x2)
qp=27*(2*a2**3/27-a1*a2/3+a0)
R=(-3*((x2+x1)**2*(2*x1+U)+(x2+x1)**2*(2*x2)+(2*x1+U)**2*(x2+x1)+(2*x1+U)**2*(2*x2)+(2*x2)**2*(x2+x1)+(2*x2)**2*(2*x1+U))+12*(x2+x1)*(2*x1+U)*(2*x2)+18*V**2*(2*(x2+x1)-(2*x1+U)-2*x2)+2*((x2+x1)**3+(2*x1+U)**3+(2*x2)**3))/54
qp_articulo=-54*R
print "Qp"
pprint(simplify(expand(Qp)))
print "Qp_articulo"
pprint(simplify(expand(Qp_articulo)))
print "qp"
pprint(simplify(expand(qp)))
print "qp_articulo"
pprint(simplify(expand(qp_articulo)))
