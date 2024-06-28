import numpy as np
import matplotlib.pyplot as plt

# Generate an array of x values from 0 to 2*pi
x = np.linspace(0, 2 * np.pi, 1000)

# Compute the y values for sin(x) and cos(x)
y_sin = np.sin(x)
y_cos = np.cos(x)

# Create the plot
plt.figure(figsize=(10, 5))

# Plot sin(x)
plt.plot(x, y_sin, label='sin(x)', color='blue')

# Plot cos(x)
plt.plot(x, y_cos, label='cos(x)', color='red')

# Add title and labels
plt.title('Plot of sin(x) and cos(x)')
plt.xlabel('x')
plt.ylabel('y')

# Add a legend
plt.legend()

# Display the plot
plt.grid(True)
plt.show()