#!/bin/sh

/home/currielab/tools_and_software/RAxML-8.1.24/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -n ml.raxml -s ml.phy -o SID5947 -f a -x 897543 -N 100 -p 345232 --silent -T 30
/home/currielab/tools_and_software/RAxML-8.1.24/raxmlHPC-PTHREADS-SSE3 -o SID5947 -f b -m GTRGAMMA -z RAxML_bootstrap.ml.raxml -t RAxML_bestTree.ml.raxml -n ml.BS_TREE -T 30
