#!/bin/bash
#SBATCH --output=trackmate_%j.out
#SBATCH --error=trackmate_%j.err
#SBATCH --time=02:00:00
#SBATCH --account=cavitation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --partition=l40s_normal_q
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kbsanjayvasanth@vt.edu
#SBATCH --job-name=trackanalysis

set -euo pipefail

module reset
module load MATLAB

WORKDIR="/home/kbsanjayvasanth/Tracking dataanlaysis/Code"
ENTRY_SCRIPT=${ENTRY_SCRIPT:-"$WORKDIR/main_batch_trackmate_arc.m"}

cd "$WORKDIR" || { echo "Failed to cd to WORKDIR"; exit 1; }
if [[ ! -f "$ENTRY_SCRIPT" ]]; then
    echo "MATLAB entry script not found: $ENTRY_SCRIPT"
    exit 1
fi

echo "============================================================"
echo " SLURM Job ID:   $SLURM_JOB_ID"
echo " Node:           $(hostname)"
echo " Workdir:        $WORKDIR"
echo " Entry script:   $ENTRY_SCRIPT"
echo " Start time:     $(date)"
echo " Cores:          $SLURM_CPUS_PER_TASK"
echo "============================================================"

export TMPDIR="$WORKDIR/tmp_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

export MATLAB_LOG_DIR="$TMPDIR"
export MATLAB_PREFDIR="$TMPDIR/matlab_prefs"
mkdir -p "$MATLAB_PREFDIR"

JAVA_HEAP_GB=${JAVA_HEAP_GB:-12}
JAVA_OPTS_FILE_WORKDIR="$WORKDIR/java.opts"
JAVA_OPTS_FILE_PREFDIR="$MATLAB_PREFDIR/java.opts"
for JAVA_OPTS_FILE in "$JAVA_OPTS_FILE_WORKDIR" "$JAVA_OPTS_FILE_PREFDIR"; do
cat > "$JAVA_OPTS_FILE" <<EOF
-Xms512m
-Xmx${JAVA_HEAP_GB}g
-XX:+UseG1GC
EOF
done

export _JAVA_OPTIONS="-Xms512m -Xmx${JAVA_HEAP_GB}g -XX:+UseG1GC"

echo "Java heap cap: ${JAVA_HEAP_GB}G (via java.opts + _JAVA_OPTIONS)"

cleanup() {
    rm -f "$JAVA_OPTS_FILE_WORKDIR" "$JAVA_OPTS_FILE_PREFDIR"
}
trap cleanup EXIT

srun --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK \
matlab -nodisplay -nosplash -nodesktop -batch "\
fprintf('Java max heap at startup: %.3f GB\n', java.lang.Runtime.getRuntime.maxMemory/1024^3); \
try, run('$ENTRY_SCRIPT'); \
catch ME, \
fprintf(2,'MATLAB run failed. Identifier: %s\n', ME.identifier); \
try, fprintf(2,'Message: %s\n', ME.message); catch, fprintf(2,'Message: <unavailable>\n'); end; \
for si = 1:min(numel(ME.stack),10), fprintf(2,'  at %s:%d\n', ME.stack(si).name, ME.stack(si).line); end; \
exit(1); end; exit(0);"

echo "============================================================"
echo " Job finished at: $(date)"
echo "============================================================"
