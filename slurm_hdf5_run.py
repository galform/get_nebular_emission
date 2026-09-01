''' Generate and submit SLURM jobs for running GNE (get_nebular_emission) '''
from gne.gne_slurm import create_slurm_script, submit_slurm_job

verbose = True
nvol = 2
submit_jobs = False  # False to only generate scripts

# Parameter file to use as base (will modify subvols and root/snapshot)
param_file = 'run_gne_SU1.py'

# Optional: user-defined suffix for job names
# If None, suffix is derived from cutcols/mincuts/maxcuts in param_file
job_suffix = None  # e.g., 'lbol45' or None

# Simulations to process: list of (snapshot_list, subvols_list) tuples
# The simulation name is already defined in the param_file's outpath/root
taurus_runs = [
    ([87, 128], list(range(nvol))),
    #([109, 104, 98, 90, 87, 128, 96, 78], list(range(nvol))),
]

# Galform in cosma - example
cosma_runs = [
    ([39, 61], list(range(64)))
]

# Select which runs to process
runs = taurus_runs
hpc = 'taurus'

job_count = 0
for snaps, subvols in runs:
    for snap in snaps:
        # Generate SLURM script
        script_path = create_slurm_script(
            hpc, param_file, snap, subvols,
            verbose=verbose, job_suffix=job_suffix
        )
        if verbose: 
            print(f'  Created script: {script_path}')
            
        # Submit the job
        if submit_jobs:
            submit_slurm_job(script_path)
            job_count += 1

if submit_jobs and verbose:
    print(f'Total jobs submitted: {job_count}')
