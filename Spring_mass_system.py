import numpy as np
n = 4
spring_k = 89

spring_vec = np.zeros((1, 4*n))

spring_vec[0][n - 1] = spring_k

print(spring_vec)