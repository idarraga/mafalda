#!/bin/bash

if [ ${#@} -lt 1 ] ; then
    echo "use: $0 dir(string)"
    echo "John Idarraga <idarraga@cern.ch>"
    exit 1
fi

prefix="data"
dir=$1
newdir=${dir}_stripped

echo "[INFO] working at $PWD"
echo "[INFO] Making a copy of $dir into $newdir ..."
cp -Rf ./$dir ./$newdir
cd $newdir
nfiles=`ls ${prefix}*.txt | wc -l`
echo "[INFO] work on $nfiles files ${prefix}*.txt found inside $newdir"
echo "[INFO] removing trigger stamp suffix.  Filenames will be changed."
numberstrip=""
for a in ${prefix}*.txt
do
    # create the new filename
    newfn=${a:0:${#a}-14}
    newfn=${newfn}.txt
    mv $a $newfn
    
    # find the maximum number of frames
    numberstrip=${a:0:${#a}-16}
    numberstrip=${numberstrip#data}
done

# remove the zeroes at the left
maxn=`echo $numberstrip | sed 's/^[0]*//'`

echo "[INFO] max frame is : $maxn"

frameId=0
while [ $frameId -le $maxn ]
do
    # Build the name, recover the 6 digits structure
    framefn11=${prefix}`echo $frameId | awk '{ printf("%06d\n", $1) }'`

    # chip numbering goes
    # 1 2 3 4
    # 8 7 6 5
    # A python script will take care of this part
    python ../mergeOctopuce.py $framefn11
    
    # include dummy dsc file
    cp ../frame_dummyDSC.dsc "${framefn11}_assembly.txt.dsc"

    let frameId=$frameId+1
done

# move all files to its final location
cd ..
lastdir=${dir}_assembly
echo "[INFO] creating dir $lastdir"
mkdir $lastdir
mv ./${newdir}/*assembly* ./$lastdir
# erase temp dir
echo "[INFO] removing $newdir"
rm -Rf ./${newdir}
echo "[DONE] the assembly is located at $lastdir"

#for a in ${prefix}*.txt
#do
#    echo "[INFO] merging frame "
#done

#python ../mergeOctopuce