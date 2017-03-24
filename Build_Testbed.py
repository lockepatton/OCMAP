import numpy as np

filters = ['v','b']

print np.zeros(len(filters))

shifts = []
for filter_ in filters:
    shifts.append([0.,0.])

print shifts

x_,y_ = [[0,1],[2,3]]
print x_,y_

for it_,(x,y) in enumerate(zip(x_,y_)):
    print it_,x,y

x = np.array([[0, 1],[2, 3]]);
np.transpose(x)

T = np.array([np.array([0,10]),np.array([2,12]),np.array([4,14])])
print T
print np.transpose(T)
