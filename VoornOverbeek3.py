import numpy as np
import matplotlib.pyplot as plt
from scipy.special import genlaguerre

# Parameters
num_terms = 100
laguerre_polynomials = [genlaguerre(n, 1) for n in range(num_terms)]
N = 1000
y = 1-np.logspace(-5, np.log10(0.999), 200)
sigma = 0.1

# Calculates z given A and B
def hvoinv(A,B, laguerre_polynomials, num_terms=100):
    A = np.asarray(A)
    B = np.asarray(B)

    c1=1/(A+B)
    c2=(B-A)/(A+B)

    x=-c1
    for k in range(num_terms):
        L_k_minus_1_1 = laguerre_polynomials[k]  # Use precomputed Laguerre polynomial
        term = -(c2) ** (k+1) * ((c1-c1/c2)* np.exp(-(k+1)*c1) /(k+1)) * L_k_minus_1_1((k+1)*c1-(k+1)*c1/c2)

        x += term
    z=-A*x/(1+B*x)
    z=np.array(z)
    # if z.size>1:
    #     cond1 = (abs(z)+abs(B*z/A))<0.5
    #     A = np.complex128(A)
    #     B = np.complex128(B)

    #     A=A[cond1]
    #     B=B[cond1]
    #     z[cond1]=(-3 * B + np.sqrt(3) * np.sqrt(2 * A - 4 * A**2 + 3 * B**2)) / (2 * A)
    #     z[cond1] = ((-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3) / (6 * 2 ** (1 / 3) * B) -(36 * B ** 2 - 4 * A ** 2) / (3 * 2 ** (2 / 3) * B * (-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3)) -A / (3 * B))
    #     z=z.real
    # else:
    #     if (abs(z)+abs(B*z/A))<0.5:
    #         A=np.complex128(A)
    #         B=np.complex128(B)
    #         z = ((-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3) / (6 * 2 ** (1 / 3) * B) - (36 * B ** 2 - 4 * A ** 2) / (3 * 2 ** (2 / 3) * B * (-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3)) - A / (3 * B))
    #         z=z.real

    return z


# Solves for all variables given a range of y1 and tanh(y1, w2)
def master(N, sigma, y, num_terms):
    h = -3*(1/N-1) + 1/y*(1+np.sqrt((1+y)*(1-y))/2)*(1/N+sigma)*np.log((1+y)/(1-y))
    Ah = 3/2*(1+sigma)
    Bh = 3/(2*y)*(1+sigma) - (1+sigma)/y*(1+np.sqrt((1+y)*(1-y))/2)
    A = Ah/h
    B = Bh/h

    z=abs(hvoinv(A,B, laguerre_polynomials, num_terms))
    cond1 = (abs(z)+abs(B*z/A))<0.5
    A = np.complex128(A)
    B = np.complex128(B)

    A=A[cond1]
    B=B[cond1]
    z[cond1]=(-3 * B + np.sqrt(3) * np.sqrt(2 * A - 4 * A**2 + 3 * B**2)) / (2 * A)
    z[cond1] = ((-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3) / (6 * 2 ** (1 / 3) * B) -(36 * B ** 2 - 4 * A ** 2) / (3 * 2 ** (2 / 3) * B * (-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3)) -A / (3 * B))
    z=z.real

    beta = 2*z/((1+sigma)*(y+z))
    phiA = 1/2*beta*(1+y)
    phiB = 1/2*beta*(1-y)

    alpha=((1/N-1)*(phiA-phiB) - np.log((1-(1+sigma)*phiA)/(1-(1+sigma)*phiB)))/(np.sqrt(2)*(sigma**(3/2))*((phiA)**(3/2) - (phiB)**(3/2)))
    return y, z, phiA, phiB, alpha

y, z, phiA, phiB, alpha = master(N, sigma, y, num_terms)

# Plots a phase diagram of phi and alpha
def plot1(phiA, phiB, alpha):
    plt.figure(figsize=(8,6))
    plt.plot(phiA, alpha)
    plt.plot(phiB, alpha)
    index = np.argmin(phiA)
    plt.scatter([phiA[index]], [alpha[index]], color='black')
    plt.xlim(0, 1)
    plt.ylim(0, 20)
    plt.title('Phase diagram')
    plt.xlabel('phi')
    plt.ylabel('alpha')
    plt.show()
    return

plot1(phiA, phiB, alpha)

# Plots a phase diagram for multiple values of N
def multiNplot(sigma, y, num_terms):
    plt.figure(figsize=(8,6))

    y, z, phiA, phiB, alpha = master(10, sigma, y, num_terms)
    plt.plot(phiA, alpha, color='red', label='N=10')
    plt.plot(phiB, alpha, color='red')
    index = np.argmin(phiA)
    plt.scatter([phiA[index]], [alpha[index]], color='red')

    y, z, phiA, phiB, alpha = master(100, sigma, y, num_terms)
    plt.plot(phiA, alpha, color='green', label='N=100')
    plt.plot(phiB, alpha, color='green')
    index = np.argmin(phiA)
    plt.scatter([phiA[index]], [alpha[index]], color='green')

    y, z, phiA, phiB, alpha = master(1000, sigma, y, num_terms)
    plt.plot(phiA, alpha, color='blue', label='N=1000')
    plt.plot(phiB, alpha, color='blue')
    index = np.argmin(phiA)
    plt.scatter([phiA[index]], [alpha[index]], color='blue')

    plt.xlim(0, 1)
    plt.ylim(0, 20)
    plt.title('Phase diagram')
    plt.xlabel('phi')
    plt.ylabel('alpha')
    plt.legend()
    plt.show()
    return

multiNplot(0.2, y, num_terms)

# Plots a phase diagram for multiple values of N
def multisigmaplot(y, num_terms):
    plt.figure(figsize=(8,6))

    y, z, phiA, phiB, alpha = master(1000, 0.05, y, num_terms)
    plt.plot(phiA, alpha, color='red', label='sigma=0.05')
    plt.plot(phiB, alpha, color='red')
    index = np.argmin(phiA)
    plt.scatter([phiA[index]], [alpha[index]], color='red')

    y, z, phiA, phiB, alpha = master(1000, 0.1, y, num_terms)
    plt.plot(phiA, alpha, color='green', label='sigma=0.1')
    plt.plot(phiB, alpha, color='green')
    index = np.argmin(phiA)
    plt.scatter([phiA[index]], [alpha[index]], color='green')

    y, z, phiA, phiB, alpha = master(1000, 0.2, y, num_terms)
    plt.plot(phiA, alpha, color='blue', label='sigma=0.2')
    plt.plot(phiB, alpha, color='blue')
    index = np.argmin(phiA)
    plt.scatter([phiA[index]], [alpha[index]], color='blue')

    plt.xlim(0, 1)
    plt.ylim(0, 40)
    plt.xlabel('phi')
    plt.ylabel('alpha')
    plt.title('Phase diagram')
    plt.legend()
    plt.show()
    return

multisigmaplot(y, num_terms)

# Checks the validity of the solution
def check(N, sigma, phiA, phiB, alpha):
    muA = (1/N + sigma)*np.log(phiA) - (1+sigma)*np.log(1-(1+sigma)*phiA)-3*alpha*sigma*np.sqrt(2*sigma*phiA)
    muB = (1/N + sigma)*np.log(phiB) - (1+sigma)*np.log(1-(1+sigma)*phiB)-3*alpha*sigma*np.sqrt(2*sigma*phiB)
    piA = (1/N - 1)*phiA - np.log(1-(1+sigma)*phiA) - np.sqrt(2)*alpha*np.power(sigma*phiA, 1.5)
    piB = (1/N - 1)*phiB - np.log(1-(1+sigma)*phiB) - np.sqrt(2)*alpha*np.power(sigma*phiB, 1.5)

    plt.figure(figsize=(8,6))
    plt.scatter(muA, muA-muB)
    plt.xlabel('muA')
    plt.ylabel('mudiff')
    plt.title('mu')
    plt.show()

    plt.figure(figsize=(8,6))
    plt.scatter(piA, piA-piB)
    plt.xlabel('PiA')
    plt.ylabel('Pidiff')
    plt.title('Pi')
    plt.show()
    return

check(N, sigma, phiA, phiB, alpha)

def discriminant(phi, alpha, N, sigma):
    return (1/(N*phi) + 1/(1-(1+sigma)*phi) - 3/4*alpha*np.power(sigma, 2)*np.power(2*sigma*phi, -0.5)) * (1/(sigma*phi) + 1/(1-(1+sigma)*phi) - 3/4*alpha*np.power(2*sigma*phi, -0.5)) - np.power(1/(1-(1+sigma)*phi) - 3/4*alpha*sigma*np.power(2*sigma*phi, -0.5), 2)
    
def spinodal(N, sigma):
    x = np.linspace(1e-10, 0.5, 1000)
    y = np.linspace(1e-10, 20, 1000)
    X, Y = np.meshgrid(x, y)
    Z = discriminant(X, Y, N, sigma)

    plt.contour(X, Y, Z, levels=[0])
    plt.xlabel('phi')
    plt.ylabel('alpha')
    plt.show()

spinodal(1000, 0.1)