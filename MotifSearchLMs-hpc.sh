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
  echo "-BGFile <Location of background file>"
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
    -BGFile) BGFile=$2;;
  esac

  shift
done


let rows=$Metacluster1-1
let cols=$Metacluster2-1
mkdir $Analysis
#let i=37
for i in `seq 0 $rows`;
do
#let j=36
        for j in `seq 0 $cols`;
        do
                if [ ! -d $Analysis/Motifs_"$i"_"$j"_fimo ]; then
                	if [ -s $LinkFolder/Regions_"$i"_"$j" ]; then
                		echo "#$ -N MotifSearchFusion-$i-$j
				#$ -l mem_size=512
				#$ -pe openmp 1
				#$ -q sam128
				module load bedtools
				module load meme/4.12.0
				echo $LinkFolder/Regions_"$i"_"$j"
				bedtools getfasta -fi $ReferenceGenome -bed $LinkFolder/Regions_"$i"_"$j" -fo $Analysis/Motifs_"$i"_"$j".fa
				sed -i 's/\(.*\)/\U\1/g' $Analysis/Motifs_"$i"_"$j".fa
				sed -i 's/CHR/chr/g' $Analysis/Motifs_"$i"_"$j".fa
				fimo --bgfile $BGFile --qv-thresh --thresh $Thresh --parse-genomic-coord --max-stored-scores 100000 --oc $Analysis/Motifs_"$i"_"$j"_fimo $MotifDatabase $Analysis/Motifs_"$i"_"$j".fa;">$Analysis/run.sh
				let qstat=`qstat -u csjansen | wc -l`
				echo $qstat
				while [ $qstat -ge 900 ]; do
				        sleep 5m;
				        let qstat=`qstat -u csjansen | wc -l`;
				done;
	                	qsub $Analysis/run.sh
			        echo $i,$j
			fi;
		fi;
        done;
done;

