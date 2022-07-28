#!/bin/bash

function get_latency()
{
    tail $1 -n 100 | grep "00000" | awk '{print $2}'
}

function get_bandwidth()
{
    tail $1 -n 100 | grep "00000" | awk '{print $3}'
}

type="remote"
mode='R'
r_s='seq'
res_file=""

for thread in {1..32}
do
    res_file="$type-$r_s-$mode-$thread-thread.txt"
    # get_latency $res_file
    get_bandwidth $res_file
done
