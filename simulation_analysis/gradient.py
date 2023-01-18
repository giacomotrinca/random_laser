import tensorflow as tf
import numpy as np


def funct(x):
    #x^4 +y^4 + x^2 + y^2 - x -y
    x2 = tf.multiply(x,x)
    x4 = tf.multiply(x2, x2)
    #print(x)
    return tf.reduce_sum(x4) + tf.reduce_sum(x2) - tf.reduce_sum(x)
    

def grad(x):
    x2 = tf.multiply(x, x)
    x3 = tf.multiply(x, x2)
    
    return 4 * x3 + 2 * x2 -1


learning_rate = 0.0001
max_iter = 10000000
tol = 0.000001
file = open('minimization.dat', "w")

for k in range(50):
    start = tf.random.uniform(shape=[2], minval=-5., maxval=5.)

    x = start
    steps = [start]

    for i in range(max_iter):
        diff = learning_rate * grad(x)
        print(f'[{k+1}][{i+1}/{max_iter}] f({x}): {funct(x)}\t\t grad: {grad(x)}')
        file.write(f'{x[0]}\t{x[1]}\t{funct(x)}\t{k}\n')
        if tf.math.abs(tf.reduce_max(diff))<tol:
            print('Critical Point reached!')
            break
        
        x_old = x
        x = x - diff
    
        if tf.math.reduce_all(tf.equal(x_old,x)):
            print('Try to tune parameters...')
            break
    file.write("\n\n")
    #steps.append(x)
file.close()


