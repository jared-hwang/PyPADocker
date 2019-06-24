"""
This is a restart script for Warp.
It can restart a simulation that was run with the script lpa_script.py

Usage:
-----
- If warp_script.py was run in serial (python warp_script.py), then type
  python restart_after_lpa_script.py

- It warp_script.py was run in parallel (mpirun -np N python warp_script.py):
  mpirun -np N python restart_after_lpa_script.py

  (where `N` should be replaced by the number of cores used, and should
  be the same for the initial simulation and for the restart simulation.)
"""
from warp import *

# Reload the simulation from the dump files
dump_name = 'lpa_script000500'
# - Single proc case
if npes == 1:
    restart( dump_name + '.dump')
# - Parallel case
else:
    restart( dump_name )

# Proceed with 500 steps
N_steps = 500
n_stepped=0
while n_stepped < N_steps:
    step(10)
    n_stepped = n_stepped + 10
printtimers()
