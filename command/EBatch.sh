#!/bin/bash
cd ../deletion-maintenance/
./bicore -BiCore-Number-Batch $1 $2 $3 $4
./bicore -Multi-Edges-Delete $1 batch-rem.txt
cd ../insertion-maintenance/
./bicore -Multi-Edges-Insert $1 batch-ins.txt