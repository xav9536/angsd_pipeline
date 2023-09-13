###this script will split your custom list of bamfiles in several lists by population

#variables
source 01_scripts/01_config.sh

# Make a list of bamfiles by pop, as given by pop.txt
POP_FILE1=02_info/pop.txt #pop list n° 1 

cat $POP_FILE1 | while read i
do
echo $i
grep $i 02_info/bam.filelist > 02_info/"$i"bam.filelist
done

##possible to re-do this for another list of pop, beware it is not the same names, or it will overwrite
#POP_FILE2=02_info/group.txt #pop list n° 2

#cat $POP_FILE2 | while read i
#do
#echo $i
#grep $i 02_info/bam.filelist > 02_info/"$i"bam.filelist
#done
