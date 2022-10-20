
import math


# Function that calculate cross product
def cross(a, b):
    result = [a[1]*b[2] - a[2]*b[1],
              a[2]*b[0] - a[0]*b[2],
              a[0]*b[1] - a[1]*b[0]]

    return result

# Normalize Vector
def norm(a):
    module = math.sqrt(a[0]**2+a[1]**2+a[2]**2)
    V = []
    V.append(a[0]/module)
    V.append(a[1]/module)
    V.append(a[2]/module)
    return V

# Product of vectors by value
def prod(a, val):
    products = []
    for num1 in a:
        products.append(num1*val)
    return products

# Division of vector by value
def div(a, val):
    products = []
    for num1 in a:
        products.append(num1/val)
    return products

# Sum of vectors
def sum(a, b):
    result = []
    for num1, num2 in zip(a, b):
        result.append(num1+num2)
    return result

# Substraction of vectors
def sub(a, b):
    result = []
    for num1, num2 in zip(a, b):
        result.append(num1-num2)
    return result

# Modulo of a vector
def modulo(a):
    result = math.sqrt(a[0]**2+a[1]**2+a[2]**2)
    return result

# Dot product of two vectors
def dot(a, b):
    result = 0.0
    for n in range(len(a)):
        result = result+a[n]*b[n]
    return result

