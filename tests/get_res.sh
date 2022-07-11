#!/bin/bash

# get_iops filename dbname
function microbench_get_iops()
{
    echo $1 $2
    tail $1 -n 100 | grep "Metic-Operate" | awk '{print $9/1e+03}'
}

function scalability_get_write_iops()
{
    echo $1 $2
    cat $1 | grep "Metic-write" | grep "iops" | awk '{print $6/1e+03}'
}

function scalability_get_read_iops()
{
    echo $1 $2
    cat $1 | grep "Metic-read"  | grep "iops" | awk '{print $6/1e+03}'
}

# function scalability_get_nvm_write()
# {
#     echo $1 $2
#     cat $1 | grep "NVM WRITE :" | awk '{print $4/(16*10000000)}'
# }

function get_write_ampli() 
{
    echo $1
    cat $1 | grep "write ampli:" | awk '{print $3*4}'
}

function get_dram_size()
{
    echo $1
    cat $1 | grep "dram space use:" | awk '$4>"0.1"{print $4}'
}

function scalability_get_clevel_nvm_size()
{
    echo $1 $2
    cat $1 | grep "/mnt/AEP0/combotree-clevel-0" | awk '{print $3/(1024*1024*1024)}'
}

function get_multi_put()
{
    cat $1 | grep "Metic-Put" | grep "iops" | awk '{print $7/1e+06}'
}

function get_multi_get()
{
    cat $1 | grep "Metic-Get" | grep "iops" | awk '{print $7/1e+06}'
}

dbname=nnbtree
# workload=ycsb-400m
workload=longlat-400m
# workload=longtitudes-200m
# workload=lognormal-150m

# logfile="microbench-$dbname-$workload.txt"

for thread in 4 8 16 24 32
do
    logfile="multi-$dbname-$workload-zipf0.9-th$thread.txt"
    get_multi_put $logfile $dbname
    # get_multi_get $logfile $dbname
done

# scalability_get_write_iops $logfile $dbname

# scalability_get_read_iops $logfile $dbname
# microbench_get_iops $logfile $dbname

# dbs="letree fastfair pgm xindex lipp alex"

# for dbname in $dbs; do
#     echo "$dbname"
#     logfile="nvm-write-$dbname-$workload.txt"
#     scalability_get_nvm_write $logfile $dbname
# done

# if [ $# -ge 1 ]; then
#     dbname=$1
# fi

#scalability_get_write_iops $logfile $dbname

#scalability_get_nvm_write $logfile $dbname
#scalability_get_nvm_size $logfile $dbname
#scalability_get_clevel_nvm_size $logfile $dbname

# logfile="microbench-$dbname-$workload.txt"
# microbench_get_iops $logfile $dbname