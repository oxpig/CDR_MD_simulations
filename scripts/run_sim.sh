# shell script to run a single antibody simulation
# CONDA_PATH=$(which conda)
echo $CONDA_EXE
eval "$($CONDA_EXE shell.bash hook)"

# Check for pdb NAME and chain arguments
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]; then
    echo "Error: Missing arguments. Usage: $0 <NAME> <INSTRUC> <CHAINH> <CHAINL> <OUTDIR>"
    exit 1
fi

NAME=$1
INSTRUC=$2
CHAINH=$3
CHAINL=$4
OUTDIR=$5

echo "[INFO] Starting simulation for $NAME"

LOGFILE=${OUTDIR}/simulation.log


# Create log file and output directory
mkdir -p $OUTDIR
exec > >(tee -a "$LOGFILE") 2>&1

echo "[INFO] Simulating $NAME with H chain $CHAINH and L chain $CHAINL, saving to $OUTDIR"

echo "[INFO] Copying imgt numbered pdb"
cp "$INSTRUC" "$OUTDIR/${NAME}_${CHAINH}${CHAINL}.pdb"


echo "[INFO] Renumbering and relaxing pdb"

conda activate prep_env
python CALVADOS3_Fv/CV3_Fv_pipeline/imgt2relax2calvados3.py --pdb-file ${OUTDIR}/${NAME}_${CHAINH}${CHAINL}.pdb --chains $CHAINH,$CHAINL --output-dir $OUTDIR  --job-id $NAME --overwrite --relax


echo "[INFO] Creating dssp file"

python CALVADOS3_Fv/CV3_Fv_pipeline/create_dssp.py --infile ${OUTDIR}/${NAME}_Fv_${CHAINH}${CHAINL}.pdb --outdir $OUTDIR
conda deactivate


conda activate calvados

echo "[INFO] Creating folded constraints"
mkdir ${OUTDIR}/constraints_folded
python CALVADOS3_Fv/CV3_Fv_pipeline/relaxed_dssp2constraints.py --dssp-file ${OUTDIR}/${NAME}_Fv_${CHAINH}${CHAINL}.dssp --constraints-dir ${OUTDIR}/constraints_folded/ --map-file ${OUTDIR}/${NAME}_dict_Fv_${CHAINH}${CHAINL}.txt --chain ${CHAINH}${CHAINL} --energy-threshold -0.01

echo "[INFO] Updating folded constraints for CDR1"
python CALVADOS3_Fv/CV3_Fv_pipeline/add_cdr1_constraints.py --dssp-file ${OUTDIR}/${NAME}_Fv_${CHAINH}${CHAINL}.dssp --map-file ${OUTDIR}/${NAME}_dict_Fv_${CHAINH}${CHAINL}.txt --constraints-file ${OUTDIR}/constraints_folded/${NAME}_Fv_${CHAINH}${CHAINL}_constraints.yaml --pdb ${NAME} --chain ${CHAINH}${CHAINL}

echo "[INFO] Creating CDR constraints"
python CALVADOS3_Fv/CV3_Fv_pipeline/dssp2ij_energy_relaxed.py --dssp-file ${OUTDIR}/${NAME}_Fv_${CHAINH}${CHAINL}.dssp --out-dir $OUTDIR --map-file ${OUTDIR}/${NAME}_dict_Fv_${CHAINH}${CHAINL}.txt --chains ${CHAINH}${CHAINL} --energy-threshold -0.01

echo "[INFO] Creating constraints for contacts"
python CALVADOS3_Fv/CV3_Fv_pipeline/pdb2ij_energy_relaxed.py --pdb-file ${OUTDIR}/${NAME}_Fv_${CHAINH}${CHAINL}.pdb --energy-HB-file ${OUTDIR}/${NAME}_Fv_${CHAINH}${CHAINL}_ij_energy.csv --out-dir ${OUTDIR} --map-file ${OUTDIR}/${NAME}_dict_Fv_${CHAINH}${CHAINL}.txt --chains ${CHAINH}${CHAINL}

echo "[INFO] Preparing simulation"
python CALVADOS3_Fv/prepare_ABrun.py --input-pdb ${OUTDIR}/${NAME}_Fv_${CHAINH}${CHAINL}.pdb --fdomains-file ${OUTDIR}/constraints_folded/${NAME}_Fv_${CHAINH}${CHAINL}_constraints.yaml --custom-restraints-file ${OUTDIR}/${NAME}_Fv_${CHAINH}${CHAINL}_ij_energy_contact.csv --out-dir $OUTDIR --fresidues-file /vols/opig/projects/cagiada-CDRsimulations/analysis_files/CDRsimulations/CALVADOS/residues_C3.csv #--steps 3010

echo "[INFO] Running simulation"
cd $OUTDIR
python run.py
conda deactivate

echo "[INFO] Done"
