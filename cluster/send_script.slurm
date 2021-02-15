#!/bin/bash
#
# Add description
#
#SBATCH --job-name=my_run
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --mail-user=moritz.oberrauch@student.uibk.ac.at
#SBATCH --mail-type=ALL

# Abort whenever a single step fails.
# Without this, bash will just continue on errors.
set -e

# On every node, when slurm starts a job, it will make sure the directory
# /work/username exists and is writable by the jobs user.
# We create a sub-directory there for this job to store its runtime data at.
WORKDIR="/work/$SLURM_JOB_USER/$SLURM_JOB_ID/my_run/"
mkdir -p "$WORKDIR"
echo "Workdir for this run: $WORKDIR"

# Similarly, let's create a dedicated output directory
OUTDIR="/home/users/moberrauch/run_output/my_run/"
mkdir -p "$OUTDIR"
echo "Outdir for this run: $OUTDIR"

# Export the WORKDIR and OUTDIR as environment variable
# so our script can use it to find its working directory.
export WORKDIR
export OUTDIR

# Use the local data download cache
export OGGM_DOWNLOAD_CACHE=/home/data/download
export OGGM_DOWNLOAD_CACHE_RO=1
export OGGM_EXTRACT_DIR="/work/$SLURM_JOB_USER/$SLURM_JOB_ID/oggm_tmp"

# All commands in the EOF block run inside of the container. Adjust container
# version to your needs, they are guaranteed to never change after their
# respective day has passed.
srun -n 1 -c "${SLURM_JOB_CPUS_PER_NODE}" singularity exec /home/users/moberrauch/oggm_20210122.sif bash -s <<EOF
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
  # Install needed dependencies
  pip install sklearn
  # OPTIONAL: install OGGM latest and OGGM-VAS
  pip install --no-deps "git+https://github.com/OGGM/oggm.git@c37ba50058c4ad7f66f8332bbf49a786c9f46506
  "
  pip install --no-deps "git+https://github.com/OGGM/oggm-vas.git@3a08963344af4cfb17398947676ccb9b51d61cd9
  "
  # Finally, the run
  python3 my_script.py
EOF

# OPTIONAL: Copy files from working directory to home directory
cp "$WORKDIR"/*.csv /home/users/moberrauch/run_output/ 2>/dev/null || :
cp "$WORKDIR"/*.nc /home/users/moberrauch/run_output/ 2>/dev/null || :

# Print a final message so you can actually see it being done in the output log.
echo "SLURM DONE"
