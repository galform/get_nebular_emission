''' Auxiliary functions for SLURM submission - HDF5 input processing '''
import os, sys
import re, glob
import subprocess
import gne.gne_const as c
from datetime import datetime

def generate_job_name(model, snap, job_suffix=None):
    """
    Generate an informative prefix for job names
    
    Parameters
    ----------
    model : string
        Identifier of the model
    snap : int
        Snapshot number
    job_suffix : string or None
        User-defined suffix. 
        
    Returns
    -------
    job_name : string
        Informative job name
    """
    # Get suffix from cutcols if not provided
    if job_suffix is None:
        job_name = f'gne_{model}_iz{snap}'
    else:
        job_name = f'gne_{model}_iz{snap}_{job_suffix}'
    
    return job_name


def get_slurm_template(hpc):
    """Read the SLURM template file for the specified HPC."""
    fnom = f'slurm_{hpc}_template.sh'
    template_file = os.path.join(c.slurm_temp_dir, fnom)
    # Check if template file exists
    if not os.path.exists(template_file):
        print(f'ERROR: Template file {template_file} not found')
        sys.exit()

    # Read template content
    with open(template_file, 'r') as f:
        slurm_template = f.read()
    return slurm_template


def modify_param_file(param_file, simpath, snap):
    """
    Read parameter file and modify the subvols and root lines.
    
    Parameters
    ----------
    param_file : string
        Path to the parameter file
    simpath : string
        Path to model catalogues
    snap : int
        Snapshot number to set in root path
        
    Returns
    -------
    modified_content : string
        Modified content of the parameter file
    """
    with open(param_file, 'r') as f:
        content = f.read()

    # Modify outpat line: 
    content = re.sub(
        r"^(outpath\s*=\s*).*$",
        rf"\g<1>'{simpath}'",
        content,
        flags=re.MULTILINE
    )
        
    # Modify root line: replace iz<number> with iz<snap>
    # This handles patterns like: root = ...'iz87/ivol'
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if re.match(r'\s*root\s*=', line):
            lines[i] = re.sub(r'(iz)\d+', rf'\g<1>{snap}', line)
    content = '\n'.join(lines)
    return content


def create_slurm_script(hpc, param_file, simpath, model, snap, subvols,
                        logdir=None, job_suffix=None, verbose=True):
    """
    Create a SLURM script that runs the modified parameter file.

    Parameters
    ----------
    hpc : string
        HPC machine to submit jobs
    param_file : string
        Path to the parameter file (e.g., 'run_gne_SU1.py')
    model : string
        Name of the used model
    snap : int
        Simulation snapshot number 
    subvols : string
        String with numbers indicating the range or list of subvolumes
    logdir : string
        Name of log directory
    verbose : bool
        Verbose output flag
    job_suffix : string or None
        User-defined suffix for job name.
    
    Returns
    -------
    script_path : string
        Path to the generated SLURM script
    job_name : string
        Prefix for slurm job
    """
    # Check parameter file exists
    if not os.path.exists(param_file):
        print(f'ERROR: Parameter file {param_file} not found')
        sys.exit()
    
    job_name = generate_job_name(model,snap,job_suffix=job_suffix)

    # Read the SLURM template
    slurm_template = get_slurm_template(hpc)

    # Get modified parameter file content
    modified_params=modify_param_file(param_file, simpath, snap)

    # Replace placeholders in template
    script_content = slurm_template
    script_content = script_content.replace('__GNE_LOG_DIR__', logdir)
    script_content = script_content.replace('__GNE_JOB_NAME__', job_name)
    script_content = script_content.replace('__GNE_VOLS__', subvols)
    script_content = script_content.replace('__GNE_PARAM_CONTENT__',
                                            modified_params)

    # Output directory
    if logdir is None:
        output_dir = 'logs'
    else:
        output_dir = logdir    
    os.makedirs(output_dir, exist_ok=True)

    # Write script to file
    now = datetime.now()
    date_time = now.strftime("%d_%m_%Y_%Hh_%Mm_%Ss%f")
    snom = f'submit_{job_name}_{date_time}.sh'
    script_path = os.path.join(output_dir, snom)
    with open(script_path, 'w') as f:
        f.write(script_content)

    return script_path


def submit_slurm_job(script_path,verbose=False):
    """Submit a SLURM job and return the job ID."""
    job_id = None
    
    # Check if the slurm script exists
    if not os.path.exists(script_path):
        print(f'ERROR: Script file {script_path} not found')
        return None

    try:
        # Submit the job
        process = subprocess.Popen(
            ['sbatch', script_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    except FileNotFoundError:
        print("Warning: sbatch not found. Job submission skipped.")
        return None
    
    stdout, stderr = process.communicate()
    if process.returncode == 0:
        # Extract job ID from output (format: "Submitted batch job XXXXX")
        output = stdout.decode('utf-8').strip()
        job_id = output.split()[-1] if output else 'unknown'
        if verbose:
            print(f'  Submitted {script_path}, with job ID {job_id}')
        return job_id
    else:
        print(f'  ERROR submitting {script_path}: {stderr.decode("utf-8")}')
        return None


def check_job_status(err_file, success_string='SUCCESS',verbose=True):
    """
    Check the status of a completed SLURM job.

    Parameters
    ----------
    err_file : string
        Name of the error file
    success_string : string
        String to search for in .out file to confirm success.
    verbose : bool
        If True, print detailed status messages

    Returns
    -------
    status : string
        'success' - *.err empty, .out contains success_string
        'error' - *.err file is not empty
        'incomplete' - *.out missing success_string
        'not_found' - output files not found
    """
    out_file = err_file.replace(".err",".out")
    
    # Check if files exist
    out_exists = os.path.exists(out_file)
    err_exists = os.path.exists(err_file)
    
    if not out_exists or not err_exists:
        if verbose and not out_exists:
            print(f' Log not found: {out_file}')
        elif verbose and not err_exists:
            print(f' Log not found: {err_file}')
        return 'not_found'
    
    # Check .err file (should be empty)
    has_errors = False
    error_content = None
    if err_exists:
        with open(err_file, 'r') as f:
            error_content = f.read().strip()
        if error_content:
            has_errors = True
            if verbose:
                print(f' ERROR message in {err_file}')
    
    # Check .out file for success string
    has_success = False
    if out_exists:
        with open(out_file, 'r') as f:
            out_content = f.read()
        if success_string in out_content:
            has_success = True
    
    # Determine overall status
    if has_errors:
        return 'error'
    elif has_success:
        if verbose:
            print(f'  SUCCESS for {out_file}')
        return 'success'
    else:
        if verbose:
            print(f'  No {success_string} found, incomplete run {out_file}')
        return 'incomplete'


def check_all_jobs(model, snap, logdir, job_suffix=None,
                   success_string='SUCCESS', verbose=True):
    """
    Check the status of all jobs for a list of simulations.

    Parameters
    ----------
    model : string
        Model name (e.g., 'Galform')
    snap : string
        Sanpshot number
    logdir : string
        Directory containing output files
    job_suffix : string or None
        User-defined suffix    
    success_string : string
        String to search for in .out file to confirm success
    verbose : bool
        If True, print detailed status messages

    Returns
    -------
    results : dict
        Dictionary with keys 'success', 'error', 'incomplete', 'not_found',
        each containing a list of job names in that status
    """
    results = {
        'success': [],
        'error': [],
        'incomplete': [],
        'not_found': []
    }

    if job_suffix is None:
        contains = f'{logdir}/*{model}*{snap}*'
    else:
        contains = f'{logdir}/*{model}*{snap}*{job_suffix}*'

    fnames = [f for f in glob.glob(contains) if f.endswith('.err')]
    for iname in fnames:
        status = check_job_status(iname,success_string='SUCCESS',
                                     verbose=verbose)
        results[status].append(iname)
    
    # Print summary
    if verbose:
        print('\n--- Summary ---')
        print(f'  Success:    {len(results["success"])}')
        print(f'  Error:      {len(results["error"])}')
        print(f'  Incomplete: {len(results["incomplete"])}')
        print(f'  Not found:  {len(results["not_found"])}')
    
    return results


def clean_all_jobs(model, snap, logdir, job_suffix=None,
                   only_show=True, verbose=True):
    """
    Remove .out, .err, and .sh files for all jobs in a simulation list.

    Parameters
    ----------
    model : string
        Model name (e.g., 'Galform')
    snap : string
        Snapshot name
    logdir : string
        Directory containing log files
    job_suffix : string or None
        User-defined suffix.
    only_show : bool
        If True, only list files that would be deleted.
    verbose : bool
        If True, print information about deleted files

    Returns
    -------
    deleted_files : list
        List of all files that would be or were deleted
    """
    all_deleted = []

    if job_suffix is None:
        contains = f'{logdir}/*{model}*{snap}*'
    else:
        contains = f'{logdir}/*{model}*{snap}*{job_suffix}*'

    suffixes = ('err', 'out', 'sh')
    fnames = [f for f in glob.glob(contains)
              if any(f.endswith(suf) for suf in suffixes)]
    for iname in fnames:
        if os.path.exists(iname):
            all_deleted.append(iname)
            if not only_show:
                os.remove(iname)

    # Print summary
    if verbose:
        action = 'Would delete' if only_show else 'Deleted'
        if all_deleted:
            print(f'{action} {len(all_deleted)} file(s):')
            for f in all_deleted:
                print(f'  {f}')
        else:
            print('No files to delete')
        
        if only_show and all_deleted:
            print('\n(Set only_show=False to delete.)')
    
    return all_deleted
