rm -rf Outputfolder*
rm -rf temp_user_files*
rm -rf output_and_log_dir*
cd ../

export rt="QuickVina2Docking"
sudo python autogrow_in_docker.py -j tests/${rt}.json
mv output_and_log_dir tests/output_and_log_dir.${rt}
mv temp_user_files tests/temp_user_files.${rt}
mv tests/Outputfolder tests/Outputfolder.${rt}

export rt="VinaDocking"
sudo python autogrow_in_docker.py -j tests/${rt}.json
mv output_and_log_dir tests/output_and_log_dir.${rt}
mv temp_user_files tests/temp_user_files.${rt}
mv tests/Outputfolder tests/Outputfolder.${rt}

export rt="VinaNN1Docking"
sudo python autogrow_in_docker.py -j tests/${rt}.json
mv output_and_log_dir tests/output_and_log_dir.${rt}
mv temp_user_files tests/temp_user_files.${rt}
mv tests/Outputfolder tests/Outputfolder.${rt}

export rt="VinaNN2Docking"
sudo python autogrow_in_docker.py -j tests/${rt}.json
mv output_and_log_dir tests/output_and_log_dir.${rt}
mv temp_user_files tests/temp_user_files.${rt}
mv tests/Outputfolder tests/Outputfolder.${rt}
