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

    #scalability
    # Loadname="longlat-400m"
    # date | tee scalability-tmp-${dbname}-${Loadname}.txt
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 3 -t $thread | tee -a scalability-tmp-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 3 -t $thread"

    #operation bench
    # Loadname="ycsb-400m"
    # date | tee scan-${dbname}-${Loadname}.txt
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 6 -t $thread | tee -a scan-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 6 -t $thread"

    # Loadname="longtitude-200m"
    # date | tee scan-${dbname}-${Loadname}.txt
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 2 -t $thread | tee -a scan-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 2 -t $thread"

    # Loadname="lognormal-150m"
    # date | tee scan-${dbname}-${Loadname}.txt
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 5 -t $thread | tee -a scan-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 5 -t $thread"


    # nvm write nvm size
    # Loadname="longlat-400m"
    # date | tee nvm-write-${dbname}-${Loadname}.txt
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 3 -t $thread | tee -a nvm-write-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 3 -t $thread"

    # diff entry count
    # Loadname="ycsb-400m"
    # date | tee entrycount-16-${dbname}-${Loadname}.txt
    # # gdb --args \
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 3 -t $thread | tee -a entrycount-16-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 3 -t $thread"

    # Loadname="lognormal-150m"
    # date | tee nvm-write-${dbname}-${Loadname}.txt
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 5 -t $thread | tee -a nvm-write-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 5 -t $thread"

    # Loadname="longtitude-200m"
    # date | tee nvm-write-${dbname}-${Loadname}.txt
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 2 -t $thread | tee -a nvm-write-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 2 -t $thread"

    # microbench
    rm -rf /mnt/AEP0/*
    Loadname="longlat-400m"
    date | tee microbench-${dbname}-${Loadname}.txt
    # gdb --args \
    numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    --loadstype 3 -t $thread | tee -a microbench-${dbname}-${Loadname}.txt

    echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 3 -t $thread"
    sleep 30

    rm -rf /mnt/AEP0/*
    Loadname="ycsb-400m"
    date | tee microbench-${dbname}-${Loadname}.txt
    LD_PRELOAD=libhugetlbfs.so HUGETLB_MORECORE=yes numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    --loadstype 6 -t $thread | tee -a microbench-${dbname}-${Loadname}.txt

    echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 6 -t $thread"
    sleep 30

    rm -rf /mnt/AEP0/*
    Loadname="lognormal-150m"
    loadnum=150000000
    date | tee microbench-${dbname}-${Loadname}.txt
    # gdb --args \
    LD_PRELOAD=libhugetlbfs.so HUGETLB_MORECORE=yes numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    --loadstype 5 -t $thread | tee -a microbench-${dbname}-${Loadname}.txt

    echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 5 -t $thread"
    sleep 30

    rm -rf /mnt/AEP0/*
    Loadname="longtitude-200m"
    loadnum=200000000
    date | tee microbench-${dbname}-${Loadname}.txt
    LD_PRELOAD=libhugetlbfs.so HUGETLB_MORECORE=yes numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    --loadstype 2 -t $thread | tee -a microbench-${dbname}-${Loadname}.txt

    echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 2 -t $thread"
    sleep 30

    # expandtest
    # date | tee microbench-expand-times-${dbname}-${Loadname}.txt
    # numactl --cpubind=0 --membind=0 ${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} \
    # --put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} \
    # --loadstype 3 -t $thread | tee -a microbench-expand-times-${dbname}-${Loadname}.txt

    # echo "${BUILDDIR}/microbench --dbname ${dbname} --load-size ${loadnum} "\
    # "--put-size ${opnum} --get-size ${opnum} --workload ${WorkLoad} --loadstype 3 -t $thread"
}

# DBName: combotree fastfair pgm xindex alex
function run_all() {
    dbs="letree fastfair pgm xindex"
    for dbname in $dbs; do
        echo "Run: " $dbname
        Run $dbname $1 $2 $3 1
        sleep 100
    done
}

function main() {
    dbname="combotree"
    loadnum=4000000
    opnum=1000000
    scansize=4000000
    thread=1
    if [ $# -ge 1 ]; then
        dbname=$1
    fi
    if [ $# -ge 2 ]; then
        loadnum=$2
    fi
    if [ $# -ge 3 ]; then
        opnum=$3
    fi
    if [ $# -ge 4 ]; then
        scansize=$4
    fi
    if [ $# -ge 5 ]; then
        thread=$5
    fi
    if [ $dbname == "all" ]; then
        run_all $loadnum $opnum $scansize $thread
    else
        echo "Run $dbname $loadnum $opnum $scansize $thread"
        Run $dbname $loadnum $opnum $scansize $thread
    fi 
}
# main fastfair 400000000 10000000 100000 1
# main xindex 200000000 10000000 100000 1
# main pgm 400000000 10000000 100000 1

# main xindex 400000000 10000000 100000 1
# main pgm 400000000 10000000 100000 1
# main fastfair 400000000 10000000 100000 1
# main alex 400000000 10000000 100000 1
# main lipp 400000000 10000000 100000 1

# main xindex 200000000 10000000 100000 1
# main lipp 150000000 10000000 100000 1
main nnbtree 400000000 10000000 100000 1
# main all 400000000 10000000 100000 1
# main alex 10000 1000 100000 1