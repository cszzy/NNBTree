#!/bin/bash
BUILDDIR=$(dirname "$0")/../build/
WorkLoad="/home/zzy/dataset/asia-latest.csv"
Loadname="longlat-400m"
function Run() {
    dbname=$1
    loadnum=$2
    opnum=$3
    scansize=$4
    thread=$5

    rm -rf /mnt/AEP0/*
    rm -rf /mnt/AEP1/*
    Loadname="ycsb-400m-zipf0.9"
    # Loadname="ycsb-400m"
    date | tee multi-${dbname}-${Loadname}-th${thread}.txt
    # gdb --args \
    LD_PRELOAD=libhugetlbfs.so HUGETLB_MORECORE=yes timeout 660  ${BUILDDIR}/multibench --dbname ${dbname} \
        --loadstype 3 --load-size ${loadnum} --put-size ${opnum} --get-size ${opnum} \
        -t $thread | tee -a multi-${dbname}-${Loadname}-th${thread}.txt
    echo "----------------"
    sleep 60

    # rm -rf /mnt/AEP0/*
    # rm -rf /mnt/AEP1/*
    # Loadname="longlat-400m-zipf0.9"
    # # Loadname="longlat-400m"
    # date | tee multi-${dbname}-${Loadname}-th${thread}.txt
    # # gdb --args \
    # LD_PRELOAD=libhugetlbfs.so HUGETLB_MORECORE=yes timeout 660 ${BUILDDIR}/multibench --dbname ${dbname} \
    #     --loadstype 4 --load-size ${loadnum} --put-size ${opnum} --get-size ${opnum} \
    #     -t $thread | tee -a multi-${dbname}-${Loadname}-th${thread}.txt
    # echo "----------------"
    # sleep 60

    # rm -rf /mnt/AEP0/*
    # # Loadname="longtitudes-200m"
    # Loadname="longtitudes-200m-zipf0.9"
    # loadnum=200000000
    # date | tee multi-${dbname}-${Loadname}-th${thread}.txt
    # # gdb --args \
    # timeout 660 numactl --cpubind=0 --membind=0 ${BUILDDIR}/multibench --dbname ${dbname} \
    #     --loadstype 5 --load-size ${loadnum} --put-size ${opnum} --get-size ${opnum} \
    #     -t $thread | tee -a multi-${dbname}-${Loadname}-th${thread}.txt
    # echo "----------------"
    # sleep 60

    # rm -rf /mnt/AEP0/*
    # # Loadname="lognormal-150m"
    # Loadname="lognormal-150m-zipf0.9"
    # loadnum=130000000
    # date | tee multi-${dbname}-${Loadname}-th${thread}.txt
    # # gdb --args \
    # timeout 660 numactl --cpubind=0 --membind=0 ${BUILDDIR}/multibench --dbname ${dbname} \
    #     --loadstype 6 --load-size ${loadnum} --put-size ${opnum} --get-size ${opnum} \
    #     -t $thread | tee -a multi-${dbname}-${Loadname}-th${thread}.txt
    # echo "----------------"
    # sleep 60
}

# DBName: combotree fastfair pgm xindex alex
function run_all() {
    dbs="nnbtree fastfair"
    for dbname in $dbs; do
        echo "Run: " $dbname
        Run $dbname $1 $2 $3 1
        # sleep 100
    done
}

dbname="nnbtree"
loadnum=200000000
opnum=10000000
scansize=4000000

for thread in 16
do
    Run $dbname $loadnum $opnum $scansize $thread
done

# dbname="fastfair"
# for thread in 16
# do
#     Run $dbname $loadnum $opnum $scansize $thread
# done

# if [ $# -ge 1 ]; then
#     dbname=$1
# fi
# if [ $# -ge 2 ]; then
#     loadnum=$2
# fi
# if [ $# -ge 3 ]; then
#     opnum=$3
# fi
# if [ $# -ge 4 ]; then
#     scansize=$4
# fi
# if [ $# -ge 5 ]; then
#     thread=$5
# fi
# if [ $dbname == "all" ]; then
#     run_all $loadnum $opnum $scansize $thread
# else
#     Run $dbname $loadnum $opnum $scansize $thread
# fi 

# Run combotree 4000000 100000 4000000 1