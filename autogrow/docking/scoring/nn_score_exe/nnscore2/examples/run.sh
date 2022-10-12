# You must set your own path to the Vina 1.1.2 executable
export VINA_EXEC=/Applications/autodock_vina_1_1_2_mac/bin/vina

# Run NNScore2
python ../NNScore2.py -receptor myreceptor.pdbqt -ligand myligand.pdbqt -vina_executable $VINA_EXEC

# If it works, you should get this output:
echo "====================================="
echo " CORRECT OUTPUT (TO VERIFY IT WORKS) "
echo "====================================="
echo

cat << EOF
When the poses were ranked by the best of the 20 network scores
associated with each pose, the best-scoring pose was MODEL 3 (Score =
7.9 = 12.59 nM)

When the poses were ranked by the average of the 20 network scores
associated with each pose, the best-scoring pose was MODEL 1 (Score =
1.765 +/- 2.274 = 17.18 mM). This is the recommended ranking/scoring
metric.
EOF
