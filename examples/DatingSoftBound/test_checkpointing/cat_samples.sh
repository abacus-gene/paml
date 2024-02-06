#!/usr/bin/bash

cat test_ckpt/mcmc.txt | head -n 1 > complete_sample.txt
for i in `ls -d test_ckpt*`; do
    cat $i/mcmc.txt | tail -n +2 >> complete_sample.txt
done

