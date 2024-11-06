#!/bin/bash

faa=$(echo $1 | sed -E "s/(.*)\..*/\1/").faa
faTrans $1 ${1}_tmp
sed "/>/! s~Z~*~g" ${1}_tmp > $faa
rm ${1}_tmp
