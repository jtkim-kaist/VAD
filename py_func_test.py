import tensorflow as tf
import numpy as np


def test_func(xx):
    return np.mean(xx[1:3])

x = tf.constant(np.arange(10), tf.float32)
y = tf.py_func(test_func, [x], tf.float32)

sess = tf.InteractiveSession()
yy = y.eval()

print("aa")