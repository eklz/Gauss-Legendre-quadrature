# %%
import numpy as np
import numpy.polynomial.polynomial as nppol
import matplotlib.pyplot as plt
import cmath
# %%


def Legendre(n, t):
    Pk0 = np.poly1d([1])
    Pk1 = np.poly1d([1, -1/len(t)*np.sum(t)])

    for _ in range(2, n+1):
        Bk = np.sum([i*Pk1(i)**2 for i in t])
        Bk = Bk/np.sum([Pk1(i)**2 for i in t])

        Gk = np.sum([i*Pk1(i)*Pk0(i) for i in t])
        Gk = Gk/np.sum([Pk0(i)**2 for i in t])

        Pk2 = Pk1*np.poly1d([1, -Bk]) - Gk*Pk0

        Pk0, Pk1 = Pk1, Pk2

    return Pk2
# %%


def Laguerre(m, Pm, sigma=0.001):
    k, Pk = m, Pm
    res = []

    while k > 0:
        tau0 = complex(0)
        eps = 1
        Pk1 = np.polyder(Pk)
        Pk2 = np.polyder(Pk, 2)
        Hk = (k-1)*((k-1)*Pk1**2 - k*Pk*Pk2)
        print(k)
        while eps >= sigma:
            Pk2tau = Pk2(tau0)
            Pk1tau = Pk1(tau0)
            PkTau = Pk(tau0)

            Hktau = Hk(tau0)

            if np.abs((Pk1tau+cmath.sqrt(Hktau))) >= np.abs((Pk1tau-cmath.sqrt(Hktau))):
                tau1 = tau0 - (k*PkTau/(Pk1tau+cmath.sqrt(Hktau)))

            else:
                tau1 = tau0 - (k*PkTau/(Pk1tau-cmath.sqrt(Hktau)))

            eps = np.abs(tau1-tau0)
            tau0 = tau1
            print(k, eps)
        res.append(tau1)
        Pk = (Pk/np.poly1d([1, -tau1]))[0]
        k = k-1

    return np.sort(res)


# %%
if __name__ == "__main__":
    Pm = Legendre(10, [i/1000 for i in range(1000)])
    Laguerre(10, Pm)
# %%
