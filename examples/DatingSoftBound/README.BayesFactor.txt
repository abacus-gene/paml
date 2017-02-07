Bayes factor calculation using MCMCTREE
Ziheng Yang

5 February 2017

Read the document BFDriverDOC.pdf, and replace bpp with MCMCtree. 

The folders bf1/ is for running the marginal likelihood calculation for clock=2
The folders bf2/ is for running the marginal likelihood calculation for clock=3


    cd bf1
    ../../../bin/BFdriver mcmctree.ctl 16

Check and edit the commands file, and specify the full path for mcmctree if necessary, e.g.,
     echo "../../../bin/mcmctree mcmctree.b$I.ctl > log.b$I.txt" > tmp.sh

Then run the commands to submit jobs onto the queue.

    source commands


Do the same thing in bf2/ for clock3, and collect results into the excel file: BF.Clock23.xlsx.
Note that the last line in the file commands is a grep command for retrieving the results after 
all jobs are finished:
grep BFbeta log.b*.txt
