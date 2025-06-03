
# 循环 1 到 35，提交后台任务
for i in {1..35}; do
  # nohup zsh batch_run_haps.sh "sample_dirs_${i}.txt" > "log.batch_run_haps.hisat2.${i}.log" 2>&1 &
  # nohup zsh batch_gen_seqs.sh sample_dirs_${i}.txt > log.gen_seq.${i}.log 2>&1 &
  nohup zsh batch_eva_map.sh ../sample_dirs_${i}.txt > log.batch_eva_map_sct.${i}.log 2>&1 & 
  # nohup zsh batch_eva_map_invas.sh ../sample_dirs_${i}.txt > log.batch_eva_map_invas.${i}.log 2>&1 & 
done