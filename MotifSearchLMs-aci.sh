if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-Analysis Name <Name for this Analysis>"
  echo "-LinkFolder <Link output folder location>"
  echo "-Metacluster1 <Number of metaclusters in SOM #1> "
  echo "-Metacluster2 <Number of metaclusters in SOM #2> "
  echo "-MotifDatabase <Motif IDs for deep bind> "
  echo "-ReferenceGenome <Reference genome in fa format> "
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi
while (( "$#" ));
do
  case "$1" in
    -Analysis) Analysis=$2;;
    -LinkFolder) LinkFolder=$2;;
    -Metacluster1) Metacluster1=$2;;
    -Metacluster2) Metacluster2=$2;;
    -MotifDatabase) MotifDatabase=$2;;
    -ReferenceGenome) ReferenceGenome=$2;;
    -Thresh) Thresh=$2;;
  esac

  shift
done

let rows=$Metacluster1-1
let cols=$Metacluster2-1
echo $rows
echo $cols
mkdir $Analysis
for i in `seq 0 $rows`; do
	for j in `seq 0 $cols`;
        do
                if [ ! -d $Analysis/Motifs_"$i"_"$j"_fimo ]; then
			echo "first if"
                	if [ -s $LinkFolder/Combo_"$i"_"$j" ]; then
				echo "second if"
				echo "#!/bin/bash" > run.sh
				echo "#PBS -l nodes=1:ppn=1" >> run.sh
				echo "#PBS -A open" >> run.sh
				echo "#PBS -l walltime=48:00:00" >> run.sh
				echo "#PBS -l mem=8gb">>run.sh

				echo 'cd ~/work/scripts/MotifAnalysis2' >> run.sh
				echo 'echo `pwd`' >> run.sh
				echo "echo $LinkFolder/Combo_"$i"_"$j >> run.sh
				echo "cut -f 1,2,3 $LinkFolder/Combo_"$i"_"$j" | sed 's/:/\t/g' | sed 's/-/\t/g' > $Analysis/Regions_"$i"_"$j".bed" >> run.sh
				echo "bedtools getfasta -fi $ReferenceGenome -bed $Analysis/Regions_"$i"_"$j".bed -fo $Analysis/Regions_"$i"_"$j".fa" >> run.sh
				echo "sed -i 's/\(.*\)/\U\1/g' $Analysis/Regions_"$i"_"$j".fa" >> run.sh
				echo "sed -i 's/CHR/chr/g' $Analysis/Regions_"$i"_"$j".fa" >> run.sh
				#echo "fimo --bgfile bgFile --qv-thresh --thresh $Thresh --parse-genomic-coord --max-stored-scores 20000000 --oc $Analysis/Motifs_"$i"_"$j"_fimo $MotifDatabase $Analysis/Regions_"$i"_"$j".fa" >> run.sh
				echo "fimo --bgfile hg38_bg --parse-genomic-coord --max-stored-scores 20000000 --oc $Analysis/Motifs_"$i"_"$j"_fimo $MotifDatabase $Analysis/Regions_"$i"_"$j".fa" >> run.sh
#				echo "../deepbind/deepbind $MotifDatabase < $Analysis/Regions_"$i"_"$j".fa > $Analysis/Motifs_"$i"_"$j >> run.sh
				numQueued=`qstat -u csj5166 | grep -oP "(Q|R)" | wc -l`
				echo $numQueued
				while [ $numQueued -ge 99 ]; do
					sleep 10m
					numQueued=`qstat -u csj5166 | grep -oP "(Q|R)" | wc -l`
					echo $numQueued
				done
				qsub run.sh
			fi
		fi
	done
done

