bsub -I  -b -q q_sw_expr -n 1 -cgsp 64 -share_size 4096 -host_stack 128 ./result
