#!/bin/bash

file="param"
indirect="IN"
outdirect="OUT"
savedirect="OUT/save"
log="OUT/log"
name="sym1"

if [ ! -d  $outdirect ]; then
	mkdir $outdirect
fi

if [ ! -d  $savedirect ]; then
	mkdir $savedirect
fi

maxsize=1024

if [ -e $log ]; then
	actualsize=$(du -k "$log" | cut -f 1)
	if [ $actualsize -ge $maxsize ]; then
		cp $log $log$(date +"%T")
		rm $log
		touch $log
	fi
else
	touch $log
fi

x=1
line_num=0
while read line
do
  l=$(sed -n "${x}p" "$indirect/$file")
  x=$(($x+1))
  line_num=$((line_num+1))
  jobname=$name.$line_num
  #echo $line
  start=`date +%s`
  ./progi $l $jobname
  end=`date +%s`
  runtime=$((end-start))
  echo runtime: $runtime secs >>$log
  sleep 1
done <$indirect/$file

#echo $line


