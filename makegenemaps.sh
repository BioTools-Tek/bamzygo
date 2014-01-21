#!/bin/bash

usage() {
	echo "Create a single map file containing the start and stop positions of exons genes"
	echo ""
	echo "`basename $0` (hg18|hg19) [extra args]"
	exit
}

ext=".genemap"
folder=""

if [ "$1" = "" ];then
  usage
fi

case "$1" in
hg18)
   folder="hg18"
   ;;
hg19)
   folder="hg19"
   ;;
*)
   usage
   ;;
esac

extrargs="${@:2}"

mkdir -p $folder

rm $folder/.*tmp.*

#for chrome in `seq 1 22` X Y; do
for chrome in `seq 1 22`; do
	chr="chr"$chrome
	
	file1=$folder/$folder$ext$chr"_tmp"
	file2=$folder/$folder$ext$chr"_tmp2"
	
	./bedtarget2 $chr --database $folder $extrargs >> $file1  &
	pid=$!;
	./bedtarget2 $chr --database $folder $extrargs >> $file2 2> /dev/null
	wait $pid;
	
	while [ "`diff -b $file1 $file2 | wc | awk '{print $2}'`" != "0" ]; do
		echo "\ntmp files not the same, retrying..."
		./bedtarget2 $chr --database $folder $extrargs  >> $file1 &
		pid2=$!
		./bedtarget2 $chr --database $folder $extrargs  >> $file2 2> /dev/null
		wait $pid2
	done;
#	echo -en "\ttmp files are equal!\n"
	
	cat $file2 >> $folder/$folder$ext".tmp"
	rm $file1
	rm $file2
	
	
done

#Filter out unneccasary lines
cat $folder/$folder$ext".tmp" | egrep "(browser)|(track)" -v > $folder/$folder$ext
rm $folder/$folder$ext".tmp"
dos2unix $folder/$folder$ext
mac2unix $folder/$folder$ext

echo "FIN" >&2
