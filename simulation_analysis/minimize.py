import tensorflow as tf

# Define the function f(x) and its constraints g(x)
def f(x):
    return x**3 - 6*x**2 + 4*x + 12

def g(x):
    return x

#define boundaries
min_bound = 2.
max_bound = 6.

#tuning parameters
learning_rate = 1e-4
tol = 1e-7

# Define the optimizer
optimizer = tf.optimizers.Adam(learning_rate=learning_rate)

# Define a variable for the input to the function
x = tf.Variable(3.)

# Define the constraint
constraint = tf.math.logical_and(g(x) >= min_bound, g(x) <= max_bound)

# Define the constraint as a function of x
def constraint_fn(x):
    return tf.math.logical_and(g(x) >= min_bound, g(x) <= max_bound)

# Define the optimization objective
@tf.function
def objective():
    with tf.GradientTape() as tape:
        tape.watch(x)
        loss = f(x)
    grads, = tape.gradient(loss, [x])
    optimizer.apply_gradients([(grads, x)])
    return loss, grads, x

# Optimize the function
for step in range(10000000):
    ob = objective()
    loss_value = ob[0]
    point = ob[2]
    gr = ob[1]
    tf.debugging.check_numerics(loss_value, message='Loss is inf or nan.')
    tf.Assert(constraint_fn(x), [x])
    if step % 10 == 0:
        print(f'[{step}] -> Loss = {loss_value:.4e} \tPoint = {point:.4e}\tgrad = {gr:.4e}')
    if abs(gr) <= tol:
        break
print("x: {}, grad: {}, loss: {}".format(point, gr, loss_value))