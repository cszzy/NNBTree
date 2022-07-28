

function Runlocaltest() {
#local test
# 1-32 thread
echo 1 $mode $r_s 200m pmem /mnt/AEP0 0 > test.txt
./mlc --loaded_latency -o./test.txt > local-$r_s-$mode-1-thread.txt
sleep 10

for thread in {1..15}
do
    thread_num=$[$thread+1]
    echo 0-$thread $mode $r_s 200m pmem /mnt/AEP0 0 > test.txt
    thread+=1
    ./mlc --loaded_latency -o./test.txt > local-$r_s-$mode-$thread_num-thread.txt
    sleep 10
done

echo 0-15 $mode $r_s 200m pmem /mnt/AEP0 0 > test.txt
echo 32 $mode $r_s 200m pmem /mnt/AEP0 0 >> test.txt
./mlc --loaded_latency -o./test.txt > local-$r_s-$mode-17-thread.txt
sleep 10

for thread in {33..47}
do
    thread_num=$[$thread-15]
    echo 0-15 $mode $r_s 200m pmem /mnt/AEP0 0 > test.txt
    echo 32-$thread $mode $r_s 200m pmem /mnt/AEP0 0 >> test.txt
    ./mlc --loaded_latency -o./test.txt > local-$r_s-$mode-$thread_num-thread.txt
    sleep 10
done
}

function Runremotetest() {
#remote test
# 1-32 thread
echo 1 $mode $r_s 200m pmem /mnt/AEP1 0 > test.txt
./mlc --loaded_latency -o./test.txt > remote-$r_s-$mode-1-thread.txt
sleep 10

for thread in {1..15}
do
    thread_num=$[$thread+1]
    echo 0-$thread $mode $r_s 200m pmem /mnt/AEP1 0 > test.txt
    thread+=1
    ./mlc --loaded_latency -o./test.txt > remote-$r_s-$mode-$thread_num-thread.txt
    sleep 10
done

echo 0-15 $mode $r_s 200m pmem /mnt/AEP1 0 > test.txt
echo 32 $mode $r_s 200m pmem /mnt/AEP1 0 >> test.txt
./mlc --loaded_latency -o./test.txt > remote-$r_s-$mode-17-thread.txt
sleep 10

for thread in {33..47}
do
    thread_num=$[$thread-15]
    echo 0-15 $mode $r_s 200m pmem /mnt/AEP1 0 > test.txt
    echo 32-$thread $mode $r_s 200m pmem /mnt/AEP1 0 >> test.txt
    ./mlc --loaded_latency -o./test.txt > remote-$r_s-$mode-$thread_num-thread.txt
    sleep 10
done
}

mode='R'
r_s='rand'
Runlocaltest $mode $r_s
Runremotetest $mode $r_s

mode='R'
r_s='seq'
Runlocaltest $mode $r_s
Runremotetest $mode $r_s

mode='W6'
r_s='rand'
Runlocaltest $mode $r_s
Runremotetest $mode $r_s

mode='W6'
r_s='seq'
Runlocaltest $mode $r_s
Runremotetest $mode $r_s