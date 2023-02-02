import numpy as np

from network import Network
from fc_layer import FCLayer
from activation_layer import ActivationLayer
from activations import tanh, tanh_prime
from losses import mse, mse_prime
from common import ChNbr

# training data
x_train = np.array([[[ChNbr(0.01, follow=True),ChNbr(0.01, follow=True)]],
    [[ChNbr(0.01, follow=True),ChNbr(1, follow=True)]],
    [[ChNbr(1, follow=True),ChNbr(0.01, follow=True)]],
    [[ChNbr(1, follow=True),ChNbr(1, follow=True)]]])
y_train = np.array([[[ChNbr(0.01, follow=True)]], [[ChNbr(1, follow=True)]], [[ChNbr(1, follow=True)]], [[ChNbr(0.01, follow=True)]]])

# network
net = Network()
net.add(FCLayer(2, 3))
net.add(ActivationLayer(tanh, tanh_prime))
net.add(FCLayer(3, 1))
net.add(ActivationLayer(tanh, tanh_prime))

# train
net.use(mse, mse_prime)
net.fit(x_train, y_train, epochs=1000, learning_rate=ChNbr(0.1, follow=True))

# test
out = net.predict(x_train)
print (out)
# print([[[item.val for item in arr1] for arr1 in arr] for arr in out])
