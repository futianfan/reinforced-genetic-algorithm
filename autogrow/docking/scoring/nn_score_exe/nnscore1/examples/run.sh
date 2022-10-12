# Run NNScore
python ../NNScore.py -receptor myreceptor.pdbqt -vina_output myligand.pdbqt -networks_dir ../networks/top_3_networks/

# If it works, you should get this output:
echo "====================================="
echo " CORRECT OUTPUT (TO VERIFY IT WORKS) "
echo "====================================="
echo

echo "Best score: -1.2194509130711613 (myligand.pdbqt, MODEL 1)"

echo
