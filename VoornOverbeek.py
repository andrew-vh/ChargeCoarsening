import numpy as np
import matplotlib.pyplot as plt
from scipy.special import genlaguerre

def generate_laguerre_polynomials(num_terms):
    # Generate Laguerre polynomials up to the specified number of terms
    laguerre_polynomials = [genlaguerre(n - 1, 1) for n in range(1, num_terms + 1)]
    return laguerre_polynomials

def hinv(hx, laguerre_polynomials, num_terms=100):
    logeps = np.zeros_like(hx)  # Initialize eps as an array of zeros with the same shape as hx

    # Conditions for different ranges of hx
    cond1 = hx < 1.076
    cond2 = hx > 5
    cond3 = ~cond1 & ~cond2  # hx between 1.06 and 5

    # Apply the first condition
    logeps[cond1] = np.log(1 - np.sqrt(((637875 * hx[cond1] + np.sqrt((637875 * hx[cond1] - 557172)**2 + 5833096416) - 557172)**(1/3) / (45 * 2**(1/3)) - (126 * 2**(1/3)) / (5 * (637875 * hx[cond1] + np.sqrt((637875 * hx[cond1] - 557172)**2 + 5833096416) - 557172)**(1/3)) - 7/15)))

    # Apply the second condition
    logeps[cond2] = np.log(2) - 2 * hx[cond2]

    # Apply the third condition using a vectorized approach
    hx_cond3 = hx[cond3]
    logeps_cond3 = np.zeros_like(hx_cond3)
    for n in range(num_terms):
        L_n_minus_1_1 = laguerre_polynomials[n]  # Use precomputed Laguerre polynomial
        term = -(-1) ** (n + 1) * (2 * np.exp(-2 * hx_cond3 * (n + 1)) / (n + 1)) * L_n_minus_1_1(4 * hx_cond3 * (n + 1))
        logeps_cond3 += term
    logeps_cond3 = np.log(logeps_cond3)
    logeps[cond3] = logeps_cond3

    return logeps

def hvoinv(A,B, laguerre_polynomials, num_terms=100):
    A = np.asarray(A)
    B = np.asarray(B)
    c1=1/(A-B)
    c2=-(A+B)/(A-B)
    c10=-c1/c2
    c20=1/c2

    c10=1/(A+B)
    c20=(B-A)/(A+B)

    x=-c1
    x0=-c10
    for k in range(num_terms):
        L_k_minus_1_1 = laguerre_polynomials[k]  # Use precomputed Laguerre polynomial
        term = -(c2) ** (k+1) * ((c1-c1/c2)* np.exp(-(k+1)*c1) /(k+1)) * L_k_minus_1_1((k+1)*c1-(k+1)*c1/c2)
        term0 = -(c20) ** (k+1) * ((c10-c10/c20)* np.exp(-(k+1)*c10) /(k+1)) * L_k_minus_1_1((k+1)*c10-(k+1)*c10/c20)

        x += term
        x0 +=term0
    x=-x0
    z=A*x/(1-B*x)

    cond1 = abs(z)<0.01

    A=A[cond1]
    B=B[cond1]
    z[cond1] = ((-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3) / (6 * 2 ** (1 / 3) * B) -(36 * B ** 2 - 4 * A ** 2) / (3 * 2 ** (2 / 3) * B * (-16 * A ** 3 + np.sqrt((-16 * A ** 3 - 432 * A * B ** 2 + 324 * B ** 2) ** 2 + 4 * (36 * B ** 2 - 4 * A ** 2) ** 3) - 432 * A * B ** 2 + 324 * B ** 2) ** (1 / 3)) -A / (3 * B))
    return z

# Generate Laguerre polynomials
num_terms = 100
laguerre_polynomials = generate_laguerre_polynomials(num_terms)
z=-0.01
A=0.2
B=-A/z+1/np.log((1+z)/(1-z))
# print(hvoinv(A,B, laguerre_polynomials))


x=0.1
hx=np.arctanh(x)/x
# print(1-np.exp(hinv(hx, laguerre_polynomials)))


#Specify parameters
N1=10
# Define the range and number of points for y1 and w2
y1_values = 1-np.logspace(-8, np.log10(0.99), 200)
tanhy1w2_values = np.sort(np.concatenate([np.logspace(-6, -0.0001, 100),-np.logspace(-6, -0.0001, 100)]))

# Create the grid
Y1, tanhY1W2 = np.meshgrid(y1_values, tanhy1w2_values)

y1=Y1
w2=np.arctanh(tanhY1W2)/ Y1

y2=np.tanh(np.arctanh(y1)/N1)
y3=(1 + w2) / (1/y1 + w2/y2)
w3 = ((1 + y1) * y3) / (y1 * (1 + y3)) + (w2 * (1 + y2) * y3) / (y2 * (1 + y3))

term25 = (1 / y1) + (w2 / y2) + 0.5 * np.sqrt((1 + y1) / y1 + w2 * (1 + y2) / y2) * np.sqrt((1 - y1) / y1 + w2 * (1 - y2) / y2)
term26= (1+y2)/(1-y2)
term27= ((1 + y1) / y1 + w2 * (1 + y2) / y2) / ((1 - y1) / y1 + w2 * (1 - y2) / y2)


h = (3 * (1 / N1 - 1) - term25 * (np.log(term26) + np.log(term27)))


A=-((3 * (1 + w2 + w3)) / 2)/h
B=((-3 * (1/y1 + w2 / y2 + w3 / y3) /2)+2*term25)/h


z=abs(hvoinv(A,B, laguerre_polynomials))

# print(z)
#

print(((A+B*z)*np.log((1+z)/(1-z))-z)/z)



beta1=2*z/(y1*(1+w2+w3+z/y1+w2*z/y2+w3*z/y3))
beta2=2*z*w2/(y2*(1+w2+w3+z/y1+w2*z/y2+w3*z/y3))
beta3=2*z*w3/(y3*(1+w2+w3+z/y1+w2*z/y2+w3*z/y3))

# SIGN ERROR FIXED
phi1A=0.5*beta1*(1+y1)
phi1B=0.5*beta1*(1-y1)

phi2A=0.5*beta2*(1+y2)
phi2B=0.5*beta2*(1-y2)

phi3A=0.5*beta3*(1+y3)
phi3B=0.5*beta3*(1-y3)

alpha=((1/N1 - 1) * (phi1A - phi1B) - np.log((1 - 2*phi1A - 2*phi2A) / (1 - 2*phi1B - 2*phi2B)))/(np.sqrt(2) * ((phi1A + phi2A)**(3/2) - (phi1B + phi2B)**(3/2)))

# Define the contour level to highlight
highlight_level = 4

# Plot the grid using a heatmap
plt.figure(figsize=(8, 6))
plt.contourf(Y1, tanhY1W2, alpha, levels=100, cmap='viridis')
levels = [5, 10, 20]
colors=['blue','green','red']
# Highlight the specific contour level
highlight_contour = plt.contour(Y1, tanhY1W2, alpha, levels=levels, colors=colors, linewidths=2)

plt.colorbar(label='alpha value')
plt.xlabel('y1')
plt.ylabel('tanh(w2y1)')
plt.ylim(0, 1)
plt.title('Phase diagram')
plt.show()