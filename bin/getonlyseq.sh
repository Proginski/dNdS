#!/bin/bash
# I'm not doing this trick just for fun, but because there an insoluble bug with faOneRecord
echo "$2" > tmp
faSomeRecords $1 tmp tmp.fa
cat tmp.fa | tail -n +2
rm tmp.fa

