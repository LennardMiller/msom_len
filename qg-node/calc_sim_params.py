''' script to calculate simulation parameters '''

import numpy as np

N = 2048 # grid size
N_proc = 256 # number of processors
t_phys = 24 # runtime on dahu in hours
size_out = 10 # size of output file in GB

safety = 2.5 # safety factor for computation length variations


core_time = N_proc*t_phys
print(f'Core time on dahu (max = 6400h): {core_time}h')

ppp = int(np.sqrt(N**2/N_proc)) # grid points per processor
print(f'Grid size per processor (sweet spot = 256): {ppp}')

ind = int(np.log2(ppp) - 5)

pspspp_list = [62500, 8.5e5, 1.4e6, 1.4e6] # list of performances vs ppp ( at 32, 64, ...)

pspspp = pspspp_list[ind]

i_end = pspspp*t_phys*3600*N_proc/(N**2*safety)
print(f'iend should be: {i_end:.2E}')

di_out = i_end/(size_out/(7.5e-9*N**2))
print(f'diout should be: {di_out:.2E}')