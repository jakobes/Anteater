from sympy import symbols, exp, log, diff, Rational
from sympy.printing import ccode


def cprint(arg):
    print(ccode(arg) + ";")


V, m, n, h, Ca, K, Na = symbols("x[0] x[1] x[2] x[3] x[4] x[5] x[6]")

Cm, GNa, GK, GAHP, GKL, GNaL, GClL, GCa, Gglia, Koinf, gamma1, tau, control = symbols(
    "Cm GNa GK GAHP GKL GNaL GClL GCa Gglia Koinf gamma1 tau control"
)

rho, eps0, beta0, ECa, Cli, Clo, ECl, phi = symbols("rho eps0 beta0 ECa Cli Clo ECl phi")


Ip1, Ip2, Ip3 = symbols("25. 3. 5.5")
Ipump = rho*(1/(1 + exp(Ip1/Ip2 - Na/Ip2)))*(1/(1 + exp(Ip3 - K)))

Ig1, Ig2 = symbols("18./2.5 2.5")
IGlia = Gglia/(1 + exp(Ig1 - K/Ig2))

Idiff = eps0*(K - Koinf)        # NB! Koinf varying in space

Ep1 = symbols("26.64")
Ki = 140 + (18 - Na)
Nao = 144 - beta0*(Na - 18)
ENa = Ep1*log(Nao/Na)
EK = Ep1*log(control*K/Ki)
ECl = Ep1*log(Cli/Clo)

am1, am2, bm1, bm2 = symbols("0.1 3 55./18. 1./18.")
am = (3.0 + am1*V)*(1 - exp(-am2 - am1*V))**(-1)
bm = 4*exp(-bm1 - bm2*V)

ah1, ah2, ah3, bh1, bh2 = symbols("0.07 2.2 0.05 1.4 0.1")
ah = ah1*exp(-ah2 - ah3*V)
bh = (1 + exp(-bh1 - bh2*V))**(-1)

an1, an2, an3, an4, bn1, bn2, bn3 = symbols("0.34 0.01 3.4 0.1 0.125 0.55 0.0125")
an = (an1 + an2*V)*(1 - exp(-an3 - an4*V))**(-1)
bn = bn1*exp(-bn2 - bn3*V)

taum = 1./(am + bm)
minf = 1./(am + bm)*am
tauh = 1./(bh + ah)
hinf = 1./(bh + ah)*ah
taun = 1./(an + bn)
ninf = 1./(an + bn)*an

INa = GNa*m**3*h*(V - ENa) + GNaL*(V - ENa)
IK = (GK*n**4 + GAHP*Ca/(1 + Ca) + GKL)*(V - EK)
ICl = GClL*(V - ECl)

F0 = -(INa + IK + ICl)/Cm
F1 = phi*(minf - m)/taum
F2 = phi*(ninf - n)/taun
F3 = phi*(hinf - h)/tauh

F41, F42, F43, F44 = symbols("80. 0.002 25. 2.5")
F4 = Ca/F41 - F42*GCa*(V - ECa)/(1 + exp(-(V + F43)/F44))
F5 = (gamma1*beta0*IK - 2*beta0*Ipump - IGlia - Idiff)/tau
F6 = -(gamma1*INa + 3*Ipump)/tau

variables = (V, m, n, h, Ca, K, Na)

for var in variables:
    cprint(diff(F6, var))
