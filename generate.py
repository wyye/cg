#!/usr/bin/env python

import sys
import random

with open(sys.argv[1], 'w') as f:
    for i in range(int(sys.argv[2])):
        for j in range(4):
            f.write('%f ' % random.uniform(-10, 10))
        f.write('\n')
        
