#!/bin/bash

#批处理，也可以使用 parallel

input="./batch_date.txt"


while IFS= read -r line
do
  ./build/linux/x86_64/release/c $line
  #nohup ./build/linux/x86_64/release/c $line & 
done < "$input"

 
