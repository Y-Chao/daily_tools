#!/usr/bin/bash

read -p "Source files' suffix: " suffix

for sf in *.${suffix};do
    dir_name=${sf%.*}
    echo ${dir_name}
    mkdir ${dir_name}
    mv $sf ${dir_name}/init.${suffix}
done