qsub -l select=1:ncpus=64,walltime=30:00,filesystems=grand:home -q debug -A remote_offloading run.sh
