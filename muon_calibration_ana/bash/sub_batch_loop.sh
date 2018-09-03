for index in {1..9}
  do
  echo $index
  bsub -R "pool>80000" -q 1nd run_test.sh $index
done
