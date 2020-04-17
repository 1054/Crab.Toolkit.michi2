#!/usr/bin/env python

import numpy as np
from copy import copy, deepcopy

aa = []
a = {'var':np.array([1,1]),}
aa.append(a)

b = copy(a)
b['var'] += np.array([1,1])
print(a, b) # a value changed





aa = []
a = {'var':np.array([1,1]),}
aa.append(a)

b = copy(aa[0])
b['var'] += np.array([2,2])
print(a, b) # a value changed





aa = []
a = {'var':np.array([1,1]),}
aa.append(a)

b = copy(aa[0])
b['var'] = np.array([1,1]) + np.array([2,2])
print(a, b) # a value unchanged





aa = []
a = {'var':np.array([1,1]),}
aa.append(a)

b = deepcopy(aa[0])
b['var'] += np.array([3,3])
print(a, b) # a value unchanged






