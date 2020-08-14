#!/bin/bash
#
#SBATCH --job-name=test_run
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --mail-user=moritz.oberrauch@student.uibk.ac.at
#SBATCH --mail-type=ALL

# Abort whenever a single step fails.
# Without this, bash will just continue on errors.
# set -e

# On every node, when slurm starts a job, it will make sure the directory
# /work/username exists and is writable by the jobs user.
# We create a subdirectory there for this job to store its runtime data at.
WORKDIR="/work/$SLURM_JOB_USER/$SLURM_JOB_ID/my_run"
mkdir -p "$WORKDIR"
echo "Workdir for this run: $WORKDIR"

# Export the WORKDIR as environment variable
# so our script can use it to find its working directory.
export WORKDIR

# Use the local data download cache
export OGGM_DOWNLOAD_CACHE=/home/data/download
export OGGM_DOWNLOAD_CACHE_RO=1
export OGGM_EXTRACT_DIR="/work/$SLURM_JOB_USER/$SLURM_JOB_ID/oggm_tmp"

# All commands in the EOF block run inside of the container. Adjust container
# version to your needs, they are guaranteed to never change after their
# respective day has passed.
srun -n 1 -c "${SLURM_JOB_CPUS_PER_NODE}" singularity exec docker://oggm/oggm:20200802 bash -s <<EOF
  set -e
  # Setup a fake home dir inside of our workdir, so we don't clutter the actual
  # shared homedir with potentially incompatible stuff.
  export HOME="$WORKDIR/fake_home"
  mkdir "\$HOME"
  # Create a venv that _does_ use system-site-packages, since everything is
  # already installed on the container. We cannot work on the container itself,
  # as the base system is immutable.
  python3 -m venv --system-site-packages "$WORKDIR/oggm_env"
  source "$WORKDIR/oggm_env/bin/activate"
  # Make sure latest pip is installed
  pip install --upgrade pip setuptools
  # OPTIONAL: install OGGM latest and OGGM-VAS
  pip install --no-deps "git+https://github.com/OGGM/oggm.git@5d61f9dd48cd4e389ca092d4b1865ff28e0499c0"
  pip install --no-deps "git+https://github.com/OGGM/oggm-vas.git@3327888b24f86dd047481a3d9c2f8a5d6e405b65"
  # Finally, the run
  python3 ./my_eq_run.py
EOF

# Copy files
cp "$WORKDIR"/*.nc /home/users/moberrauch/run_output/
# cp "$WORKDIR"/*.csv /home/users/moberrauch/run_output/

# Print a final message so you can actually see it being done in the output log.
echo "SLURM DONE"