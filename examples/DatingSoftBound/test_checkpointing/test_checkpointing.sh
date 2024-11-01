#!/usr/bin/bash

cd test_ckpt
echo "ENTERING THE FIRST ROUND, WHERE WE WILL SAVE THE CHECKPOINT"
# run mcmctree
../../../../src/./mcmctree
cd ..
# copy to a new dir where we will resume
cp -r test_ckpt test_ckpt.2
cd test_ckpt.2
echo "ENTERING THE SECOND ROUND, WHERE WE WILL RESUME FROM THE PREVIOUS CHECKPOINT AND KEEP SAVING"
# replace the value of checkpoint
sed -i 's/checkpoint\ =\ 1/checkpoint\ =\ 2/g' mcmctree.ctl
# clean files from the previous round
rm FigTree.tre mcmc.txt out.txt SeedUsed
# run again mcmctree which should resume from mcmctree.ckpt
../../../../src/./mcmctree
cd ..
# copy to a new dir where we will resume
cp -r test_ckpt.2 test_ckpt.3
cd test_ckpt.3
echo "ENTERING THE THIRD ROUND, WHERE WE WILL RESUME FROM THE PREVIOUS CHECKPOINT AND KEEP SAVING"
# clean files from the previous round
rm FigTree.tre mcmc.txt out.txt SeedUsed
# run again mcmctree which should resume from mcmctree.ckpt
../../../../src/./mcmctree
cd ..
# copy to a new dir where we will resume
cp -r test_ckpt.3 test_ckpt.4
cd test_ckpt.4
echo "ENTERING THE FOURTH ROUND, WHERE WE WILL RESUME FROM THE PREVIOUS CHECKPOINT AND KEEP SAVING"
# clean files from the previous round
rm FigTree.tre mcmc.txt out.txt SeedUsed
# run again mcmctree which should resume from mcmctree.ckpt
../../../../src/./mcmctree
