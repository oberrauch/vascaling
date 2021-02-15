#!/bin/bash
#
# Running RGI region 11 with CESM Data
#
#SBATCH --job-name=cesm_run
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --mail-user=moritz.oberrauch@student.uibk.ac.at
#SBATCH --mail-type=ALL

# Abort whenever a single step fails. Without this, bash will just continue on errors.
set -e

# Current RGI region
OGGM_RGI_REG=`printf "%02d" $SLURM_ARRAY_TASK_ID`
export OGGM_RGI_REG
echo "RGI Region: $OGGM_RGI_REG"

# On every node, when slurm starts a job, it will make sure the directory
# /work/username exists and is writable by the jobs user.
# We create a sub-directory there for this job to store its runtime data at.
OGGM_WORKDIR="/work/$SLURM_JOB_USER/$SLURM_JOB_ID/cesm_run"
mkdir -p "$OGGM_WORKDIR"
export OGGM_WORKDIR
echo "Workdir for this run: $OGGM_WORKDIR"

# Similarly, let's create a dedicated output directory
OUTDIR="/home/users/moberrauch/run_output/cesm_run/"
mkdir -p "$OUTDIR"
export OUTDIR
echo "Outdir for this run: $OUTDIR"

# Use the local data download cache
export OGGM_DOWNLOAD_CACHE=/home/data/download
export OGGM_DOWNLOAD_CACHE_RO=1
export OGGM_EXTRACT_DIR="/work/$SLURM_JOB_USER/$SLURM_JOB_ID/oggm_tmp"

# Try to make mp better
# export OGGM_USE_MP_SPAWN=1

# Link www fmaussion data here to avoid useless downloads
mkdir -p "$OGGM_WORKDIR/cache/cluster.klima.uni-bremen.de"
ln -s /home/www/fmaussion "$OGGM_WORKDIR/cache/cluster.klima.uni-bremen.de/~fmaussion"

# Add other useful defaults
export LRU_MAXSIZE=1000

# All commands in the EOF block run inside of the container
# Adjust container version to your needs, they are guaranteed to never change after their respective day has passed.
srun -n 1 -c "${SLURM_JOB_CPUS_PER_NODE}" singularity exec /home/users/moberrauch/oggm_20210122.sif bash -s <<EOF
    set -e
  # Setup a fake home dir inside of our workdir, so we don't clutter the actual shared homedir with potentially incompatible stuff.
  export HOME="$OGGM_WORKDIR/fake_home"
  mkdir "\$HOME"
  # Create a venv that _does_ use system-site-packages, since everything is already installed on the container.
  # We cannot work on the container itself, as the base system is immutable.
  python3 -m venv --system-site-packages "$OGGM_WORKDIR/oggm_env"
  source "$OGGM_WORKDIR/oggm_env/bin/activate"
  # Make sure latest pip is installed
  pip install --upgrade pip setuptools
  # OPTIONAL: install OGGM latest
  pip install --no-deps "git+https://github.com/OGGM/oggm.git@a91f577f32654787710fcf3ffc8d7056fdc2ee9d"
  # Increase number of allowed open file descriptors
  ulimit -n 65000
  # Finally, the run
  python my_cmip6.py "$1"
EOF

# Write out
echo "Copying files..."
rsync -avzh "$OUTDIR/" cmip6_output

# Print a final message so you can actually see it being done in the output log.
echo "SLURM DONE"