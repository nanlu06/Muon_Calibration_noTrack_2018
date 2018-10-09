declare -a arr=("SingleMuon_2018Dv2_ieta27_Oct2_test")
## now loop through the above array
for i in "${arr[@]}"
do
 count=0
 while read p; do
   echo "$p"
   count=$((count+1))
   prefix=${i}_${count}
   echo "$p" > MuonCali_test/${prefix}_input_list.txt
   input=${prefix}_input_list.txt
   path=/eos/cms/store/group/dpg_hcal/comm_hcal/nlu/ntuples/Sept12018/HCALSample/
   echo $input $count $path
   bsub -R "pool>200000" -q 1nh run_test.sh $input $prefix $path
 done < datalist/${i}.txt
done
