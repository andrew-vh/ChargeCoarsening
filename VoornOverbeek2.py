import numpy as np
import matplotlib.pyplot as plt
from scipy.special import genlaguerre

# Parameters
num_terms = 100
laguerre_polynomials = [genlaguerre(n, 1) for n in range(num_terms)]
N1 = 100
y1_values = 1-np.logspace(-2, np.log10(0.999), 200)
tanhy1w2_values = np.sort(np.concatenate([np.logspace(-8, -0.0001, 200),-np.logspace(-8, -0.0001, 200)]))
alpha_crit = 3

# Calculates z given A and B
def hvoinv(A,B, laguerre_polynomials, num_terms):
    A = np.asarray(A)
    B = np.asarray(B)
    c1 = 1/(A+B)
    c2 = (B-A)/(A+B)

    x = -c1

    for k in range(num_terms):
        L_k_minus_1_1 = laguerre_polynomials[k]
        term = -(c2) ** (k+1) * ((c1-c1/c2)* np.exp(-(k+1)*c1) /(k+1)) * L_k_minus_1_1((k+1)*c1-(k+1)*c1/c2)
        x += term
        print(k)
    
    z=-A*x/(1+B*x)

    cond1 = abs(z)<0.01

    A=A[cond1]
    B=B[cond1]
    z[cond1] = ((-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3) / (6 * 2 ** (1 / 3) * B) -(36 * B ** 2 - 4 * A ** 2) / (3 * 2 ** (2 / 3) * B * (-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3)) -A / (3 * B))
    return z

# Solves for all variables given a range of y1 and tanh(y1, w2)
def master(N1, y1_values, tanhy1w2_values, num_terms):
    y1, tanhy1w2 = np.meshgrid(y1_values, tanhy1w2_values)
    w2 = np.arctanh(tanhy1w2)/y1

    y2=np.tanh(np.arctanh(y1)/N1)
    y3=(1 + w2) / (1/y1 + w2/y2)
    w3 = ((1 + y1) * y3) / (y1 * (1 + y3)) + (w2 * (1 + y2) * y3) / (y2 * (1 + y3))

    term25 = (1 / y1) + (w2 / y2) + 0.5 * np.sqrt((1 + y1) / y1 + w2 * (1 + y2) / y2) * np.sqrt((1 - y1) / y1 + w2 * (1 - y2) / y2)
    term26= (1+y2)/(1-y2)
    term27= ((1 + y1) / y1 + w2 * (1 + y2) / y2) / ((1 - y1) / y1 + w2 * (1 - y2) / y2)

    h = (3 * (1 / N1 - 1) - term25 * (np.log(term26) + np.log(term27)))

    A=-((3 * (1 + w2 + w3)) / 2)/h
    B=((-3 * (1/y1 + w2 / y2 + w3 / y3) /2)+2*term25)/h

    z=abs(hvoinv(A,B, laguerre_polynomials, num_terms))

    beta1=2*z/(y1*(1+w2+w3+z/y1+w2*z/y2+w3*z/y3))
    beta2=2*z*w2/(y2*(1+w2+w3+z/y1+w2*z/y2+w3*z/y3))
    beta3=2*z*w3/(y3*(1+w2+w3+z/y1+w2*z/y2+w3*z/y3))

    phi1A=0.5*beta1*(1+y1)
    phi1B=0.5*beta1*(1-y1)
    phi2A=0.5*beta2*(1+y2)
    phi2B=0.5*beta2*(1-y2)
    phi3A=0.5*beta3*(1+y3)
    phi3B=0.5*beta3*(1-y3)

    alpha=((1/N1 - 1) * (phi1A - phi1B) - np.log((1 - 2*phi1A - 2*phi2A) / (1 - 2*phi1B - 2*phi2B)))/(np.sqrt(2) * ((phi1A + phi2A)**(3/2) - (phi1B + phi2B)**(3/2)))
    return y1, tanhy1w2, y2, w2, y3, w3, z, phi1A, phi1B, phi2A, phi2B, phi3A, phi3B, alpha

y1, tanhy1w2, y2, w2, y3, w3, z, phi1A, phi1B, phi2A, phi2B, phi3A, phi3B, alpha = master(N1, y1_values, tanhy1w2_values, num_terms)

# Plots a phase diagram of y1 and tanh(y1w2)
def plot1(y1, tanhy1w2, alpha, alpha_crit):
    plt.figure(figsize=(8,6))
    plt.contourf(y1, tanhy1w2, alpha, cmap='viridis', levels=100)
    highlight = plt.contour(y1, tanhy1w2, alpha, levels=[0.99*alpha_crit, alpha_crit, 1.01*alpha_crit], colors='red', linewidths=0.5)
    plt.colorbar(label='alpha')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.title('Phase diagram')
    plt.show()
    return

plot1(y1, tanhy1w2, alpha, alpha_crit)

# Selects variables that give the critical alpha
def selectcrit(y1, tanhy1w2, y2, w2, y3, w3, z, phi1A, phi1B, phi2A, phi2B, phi3A, phi3B, alpha, alpha_crit):
    rows, columns = y1.shape
    y1crit = []
    tanhy1w2crit = []
    y2crit = []
    w2crit = []
    y3crit = []
    w3crit = []
    zcrit = []
    phi1Acrit = []
    phi1Bcrit = []
    phi2Acrit = []
    phi2Bcrit = []
    phi3Acrit = []
    phi3Bcrit = []
    for x in range(rows):
        for y in range(columns):
            if (alpha[x, y]>0.995*alpha_crit and alpha[x, y]<1.005*alpha_crit and phi1A[x, y]>=0 and phi1B[x, y]>=0 and phi2A[x, y]>=0 and phi2B[x, y]>=0):
                y1crit.append(y1[x, y])
                tanhy1w2crit.append(tanhy1w2[x, y])
                y2crit.append(y2[x, y])
                w2crit.append(w2[x, y])
                y3crit.append(y3[x, y])
                w3crit.append(w3[x, y])
                zcrit.append(z[x, y])
                phi1Acrit.append(phi1A[x, y])
                phi1Bcrit.append(phi1B[x, y])
                phi2Acrit.append(phi2A[x, y])
                phi2Bcrit.append(phi2B[x, y])
                phi3Acrit.append(phi3A[x, y])
                phi3Bcrit.append(phi3B[x, y])
    return y1crit, tanhy1w2crit, y2crit, w2crit, y3crit, w3crit, zcrit, phi1Acrit, phi1Bcrit, phi2Acrit, phi2Bcrit, phi3Acrit, phi3Bcrit

y1crit, tanhy1w2crit, y2crit, w2crit, y3crit, w3crit, zcrit, phi1Acrit, phi1Bcrit, phi2Acrit, phi2Bcrit, phi3Acrit, phi3Bcrit = selectcrit(y1, tanhy1w2, y2, w2, y3, w3, z, phi1A, phi1B, phi2A, phi2B, phi3A, phi3B, alpha, alpha_crit)

# Plots a phase diagram of phi1 and phi_salt
def plot2(phi1Acrit, phi1Bcrit, phi2Acrit, phi2Bcrit):
    np1A = np.asarray(phi1Acrit)
    np1B = np.asarray(phi1Bcrit)
    np2A = np.asarray(phi2Acrit)
    np2B = np.asarray(phi2Bcrit)
    plt.figure(figsize=(8,6))
    plt.scatter(np1A, np2A)
    plt.scatter(np1B, np2B)
    for i in range(len(np1A)):
        plt.plot([np1A[i], np1B[i]], [np2A[i], np2B[i]], linestyle='--', linewidth=0.1)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.title('Phase diagram')
    plt.xlabel('phi1')
    plt.ylabel('phi_salt')
    plt.show()
    return

plot2(phi1Acrit, phi1Bcrit, phi2Acrit, phi2Bcrit)

# Checks the validity of the solution
def check(N1, phi1A, phi1B, phi2A, phi2B, alpha):
    mudiffA = np.log(phi1A) / N1 - np.log(phi2A)
    mudiffB = np.log(phi1B) / N1 - np.log(phi2B)
    mu2A = np.log(phi2A) + np.log(phi1A+phi2A) - 2*np.log(1-2*(phi1A+phi2A)) - 3*np.sqrt(2)*alpha*np.sqrt(phi1A+phi2A)
    mu2B = np.log(phi2B) + np.log(phi1B+phi2B) - 2*np.log(1-2*(phi1B+phi2B)) - 3*np.sqrt(2)*alpha*np.sqrt(phi1B+phi2B)
    piA = -np.log(1-2*(phi1A+phi2A)) - np.sqrt(2)*alpha*np.power(phi1A+phi2A, 1.5) + (1/N1 - 1)*phi1A
    piB = -np.log(1-2*(phi1B+phi2B)) - np.sqrt(2)*alpha*np.power(phi1B+phi2B, 1.5) + (1/N1 - 1)*phi1B
    
    plt.figure(figsize=(8,6))
    plt.plot(mudiffA, mudiffA, color='red')
    plt.scatter(mudiffA, mudiffB)
    plt.xlabel('mu1A - mu2A')
    plt.ylabel('mu1B - mu2B')
    plt.title('mu1 - mu2')
    plt.show()

    plt.figure(figsize=(8,6))
    plt.plot(mu2A, mu2A, color='red')
    plt.scatter(mu2A, mu2B)
    plt.xlabel('mu2A')
    plt.ylabel('mu2B')
    plt.title('mu1')
    plt.show()

    plt.figure(figsize=(8,6))
    plt.plot(piA, piA, color='red')
    plt.scatter(piA, piB)
    plt.xlabel('PiA')
    plt.ylabel('PiB')
    plt.title('Pi')
    plt.show()
    return

check(N1, phi1A, phi1B, phi2A, phi2B, alpha)