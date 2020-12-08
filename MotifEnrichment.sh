#$ -N Zscores
#$ -q sam,bio

if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-Analysis Name <Name for this Analysis>"
  echo "-LinkFolder <FusionBreakup output folder location>"
  echo "-Metacluster1 <Number of metaclusters in SOM #1> "
  echo "-Metacluster2 <Number of metaclusters in SOM #2> "
  echo "-Pval <Desired pval>"
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi

extraOptions=""

while (( "$#" ));
do
  case "$1" in
    -Analysis) Analysis=$2;;
    -LinkFolder) LinkFolder=$2;;
    -Metacluster1) Metacluster1=$2;;
    -Metacluster2) Metacluster2=$2;;
    -Pval) Pval=$2;;
    -Fimo) extraOptions=echo $extraOptions," -Fimo $2";;
  esac

  shift
done
echo $Pval
echo $extraOptions
MotifZscore/MotifZScore -Metacluster1 $Metacluster1 -Metacluster2 $Metacluster2 -FusionBreakup $LinkFolder -Analysis $Analysis/Motifs -ZScoreOutput $Analysis/ZScore_out -AllOutput $Analysis/AllOut -pval $Pval $extraOptions

