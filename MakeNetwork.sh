if [ $# -eq 0 ]    # Script invoked with no command-line args?
then
  echo "Usage: `basename $0` [required options] "
  echo "Required Options:"
  echo "-Analysis <Name for this Analysis>"
  echo "-LinkFolder <FusionBreakup output folder location>"
  echo "-Metacluster1 <Number of metaclusters in DNA SOM> "
  echo "-Metacluster2 <Number of metaclusters in RNA SOM> "
  echo "-Output <Output file location>"
  echo "-MotifFile <tsv of motif IDs to proper label> [optional]"
  exit          # Exit and explain usage.
                            # Usage: scriptname -options
                            # Note: dash (-) necessary
fi

Options=""


while (( "$#" ));
do
  case "$1" in
    -Analysis) Analysis=$2;;
    -LinkFolder) LinkFolder=$2;;
    -Metacluster1) Metacluster1=$2;;
    -Metacluster2) Metacluster2=$2;;
    -Output) Output=$2;;
    -MotifFile) Options="$Options -MotifFile $2";;
  esac

  shift
done
echo $Options
./MakeNetwork/MakeNetwork -Metacluster1 $Metacluster1 -Metacluster2 $Metacluster2 -FusionBreakup $LinkFolder -ZScoreFile $Analysis/ZScore_out -Output $Output $Options
