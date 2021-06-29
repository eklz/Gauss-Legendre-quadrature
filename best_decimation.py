from CACTUS import *
from multiprocessing import Pool
import numpy as np
from itertools import product
import pandas as pd

data = Cn2.read_csv('/scratchm/eklotz/Cn2_Tenerife_2020_fromSoundings.csv')
data.set_ground_level()
data.rm_zeros()
data.filtre(50)
data.rm_incomplete(20000)

def dec(i, j):
    tmp = data.decimate(i,j)._data
    tmp['nbsegements'] = i
    tmp['nbpoints_per_segment'] = j 
    print(f'{i},{j} ok')
    return tmp

with Pool(88) as p : 
    res = p.starmap(dec, product([i for i in range(1,11)], repeat = 2))

res1 = pd.concat(res)
res1.to_csv('/scratchm/eklotz/every_possible_decimation2.csv')
