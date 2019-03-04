from sympy import symbols, exp, log, diff, Rational
from sympy.printing import ccode

V, m, h, n, Ca, K, Na = symbols("x[0] x[1] x[2] x[3] x[4] x[5] x[6]")

Cm, GNa, GK, GAHP, GKL, GNaL, GClL, GCa, Gglia, Koinf, gamma, tau, control = symbols(
    "Cm GNa GK GAHP GKL GNaL GClL GCa Gglia Koinf gamma1 tau control"
)
rho, eps0, beta0, ECa, Cli, Clo, ECl, phi = symbols("rho eps0 beta0 ECa Cli Clo ECl phi")


Ip1, Ip2, Ip3 = symbols("25.0/3.0 5.5 1.0/3.0")
Ipump = rho*(1.0/(1.0 + exp(Ip1 - Na*Ip3)))*(1.0/(1.0 + exp(Ip2 - K)))

Ig1, Ig2 = symbols("7.2 1.0/2.5")
IGlia = Gglia/(1.0 + exp(Ig1 - K*Ig2))

Idiff = eps0*(K - Koinf)

ea1 = symbols("26.64")
Ki = 140.0 + (18.0 - Na)
Nao = 144.0 - beta0*(Na - 18.0)
ENa = ea1*log(Nao/Na)
EK = ea1*log(control*K/Ki)

am1, am2, bm1, bm2 = symbols("3.0 0.1 55.0/18.0 1.0/18.0")
a_m = (3.0 + am2*V)/(1.0 - exp(-am1 - am2*V))
b_m = 4*exp(-bm1 - bm2*V)

ah1, ah2, ah3, bh1, bh2 = symbols("2.2 0.07 0.05 1.4 0.1")
ah = ah2*exp(-ah1 - ah3*V)
bh = 1.0/(1.0 + exp(-bh1 - bh2*V))

an1, an2, bn1, bn2, an3, bn3, bn4 = symbols("3.4 0.1 0.55 0.01 0.34 0.125 0.0125")
an = (an3 + bn2*V)/(1.0 - exp(-an1 - an2*V))
bn = bn3*exp(-bn1 - bn4*V)

minf = a_m/(a_m + b_m)
taum = 1.0/(a_m + b_m)
h_inf = ah/(bh + ah)
tauh = 1.0/(bh + ah)
ninf = an/(an + bn)
taun = 1.0/(an + bn)

INa = GNa*m**3*h*(V - ENa) + GNaL*(V - ENa)
IK = (GK*n**4 + GAHP*Ca/(1.0 + Ca) + GKL)*(V - EK)
ICl = GClL*(V - ECl)


F0 = -(INa + IK + ICl)/Cm
F1 = phi*(minf - m)/taum
F2 = phi*(h_inf - h)/tauh
F3 = phi*(ninf - n)/taun

F41 = symbols("2.5")
F4 = -Ca/80.0 - 0.002*GCa*(V - ECa)/(1.0 + exp(-V/F41 + 10))
F5 = (gamma*beta0*IK - 2*beta0*Ipump - IGlia - Idiff)/tau
F6 = -(gamma*INa + 3*Ipump)/tau

current_func = F6

print(ccode(diff(current_func, V)))
print(ccode(diff(current_func, m)))
print(ccode(diff(current_func, h)))
print(ccode(diff(current_func, n)))
print(ccode(diff(current_func, Ca)))
print(ccode(diff(current_func, K)))
print(ccode(diff(current_func, Na)))

# (V, m, h, n, Ca, K, Na)
