sbatch -p yejin --nodelist=haic-hgx-7 --output=logs/fold.log bin/run_fold.sh data/structure/fold

sbatch -p yejin --nodelist=haic-hgx-7 --output=logs/inpaint.log bin/run_inpaint.sh data/structure/inpaint

sbatch -p yejin --nodelist=haic-hgx-7 --output=logs/mask.log bin/run_inpaint.sh data/structure/mask
