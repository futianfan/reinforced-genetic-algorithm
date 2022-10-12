#!/bin/bash
"""
This script will take a file path to a tab-delineated
.smi file. It will then filter it for redundancies in the 
1st and 2nd columns of the file. 
The output file is the input file + '_no_dup.smi'
"""

filename=$1
tmp_str=_tmp_dummy.smi
no_dup_str=_no_dup.smi
tmp_filename=$filename$tmp_str
new_filename=$filename$no_dup_str
echo $tmp_filename

awk '!seen[$1]++' $filename >> $tmp_filename
echo $new_filename
awk '!seen[$2]++' $tmp_filename >> $new_filename

rm $tmp_filename
