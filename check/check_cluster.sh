#!/usr/bin/env bash
#
# Call with "make testcluster"
#
# To get the result files call "./evalcheck_cluster.sh
# results/check.$TSTNAME.$BINNAME.$SETNAME.eval in directory check/
# This leads to result files
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.out
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.res
#  - results/check.$TSTNAME.$BINNMAE.$SETNAME.err

TSTNAME=$1
BINNAME=$2
SETNAME=$3
BINID=$4
TIMELIMIT=$5
NODELIMIT=$6
MEMLIMIT=$7
DISPFREQ=$8
VERSION=$9
LPS=${10}
QUEUE=${11}
CLIENTTMPDIR=${12}
OPT=${13}

ACCOUNT="dopt"
QUEUETYPE="srun"

# get current SCIP path
SCIPPATH=`pwd`

echo

# determine whether we want to run exclusively
if test "$OPT" == "opt"
then
    EXCLUSIVESTR="--exclusive"
else
    EXCLUSIVESTR=""
fi


if test ! -e $SCIPPATH/results
then
    mkdir $SCIPPATH/results
fi

# check if the settings file exists
if test ! -e $SCIPPATH/../settings/$SETNAME.set
then
    echo Skipping test since the settings file $SCIPPATH/../settings/$SETNAME.set does not exist.
    exit
fi
SETTINGS=$SCIPPATH/../settings/$SETNAME.set


# check if binary exists
if test ! -e $SCIPPATH/../$BINNAME
then
    echo Skipping test since the binary $BINNAME does not exist.
    exit
fi

# check if queue has been defined
if test "$QUEUE" == ""
then
    echo Skipping test since the queue name has not been defined.
    exit
fi

# we add 100% to the hard time limit and additional 600 seconds in case of small time limits
HARDTIMELIMIT=`expr \`expr $TIMELIMIT + 600\` + $TIMELIMIT`

# the srun queue requires a format duration HH:MM:SS (and optionally days) and megabytes.
#format is (d-)HH:MM:SS
TMP=`expr $HARDTIMELIMIT`
HARDTIMELIMIT=""
DIVISORS=(60 60 24)
for((i=0; i<=2; i++))
do
    printf -v HARDTIMELIMIT "%02d${HARDTIMELIMIT}" `expr ${TMP} % ${DIVISORS[i]}`
    # separate the numbers by colons except for the last (HH hours)
    if test $i -lt 2
    then
        HARDTIMELIMIT=":${HARDTIMELIMIT}"
    fi
    TMP=`expr ${TMP} / ${DIVISORS[i]}`
done
if test ${TMP} -gt 0
then
    HARDTIMELIMIT=${TMP}-${HARDTIMELIMIT}
fi

# we add 10% to the hard memory limit and additional 100mb to the hard memory limit
# values of MEMLIMIT are in MB
HARDMEMLIMIT=`expr \`expr $MEMLIMIT + 100\` + \`expr $MEMLIMIT / 10\``

EVALFILE=$SCIPPATH/results/check.$TSTNAME.$BINID.$QUEUE.$SETNAME.eval
echo > $EVALFILE

# counter to define file names for a test set uniquely
COUNT=1

echo "testset        " $TSTNAME
echo "queue          " $QUEUE
echo "memlimit       " $MEMLIMIT
echo "hard memlimit  " $HARDMEMLIMIT
echo "timelimit      " $TIMELIMIT
echo "hard timelimit " $HARDTIMELIMIT
echo

STORENAME=""
for i in `cat testset/$TSTNAME.test` DONE
do
  if test "$i" = "DONE"
  then
      break
  fi

  # determine network file
  NETFILE=$STORENAME".net"

  # check if we use a compressed version
  if test ! -f $SCIPPATH/$NETFILE
  then
      NETFILE=$STORENAME".net.gz"
  fi

  # try to use compressor file if it is present
  CSFILE=$STORENAME".cs"

  # store the name and read name of scenario in next loop
  if test "$STORENAME" = ""
  then
      STORENAME="$i"
      continue
  else
      STORENAME=""
  fi
  # set scenario file
  SCNFILE=$i

  # check if problem instance exists
  if test -f $SCIPPATH/$NETFILE
  then
      # compute display name for slurm
      TMPSCN=(`basename $SCNFILE .gz`)
      TMPSCN=(`basename $TMPSCN .scn`)
      TMPSCN=(${TMPSCN//\_/ })
      SHORTSCN=""
      # remove leading "nomination"
      for name in ${TMPSCN[@]}
      do
	  if [ $name == "nomination" ]
	  then
  	      continue
	  else
	      if [ -z $SHORTSCN ]
	      then
		  SHORTSCN=$name
	      else
  		  SHORTSCN=$SHORTSCN"_"$name
	      fi
	  fi
      done

      JOBNAME=$SHORTSCN
      BASENAME=$USER.$QUEUE.$TSTNAME.$COUNT"_"$SHORTSCN.$BINID.$SETNAME

      echo $SCIPPATH/results/$BASENAME >> $EVALFILE

      COUNT=`expr $COUNT + 1`

      export CLIENTTMPDIR
      export SETTINGS
      export BINNAME
      export TIMELIMIT
      export MEMLIMIT
      export NODELIMIT
      export DISPFREQ
      export BASENAME
      export SCIPPATH
      export NETFILE
      export SCNFILE
      export CSFILE

      sbatch --job-name=${JOBNAME} --mem=$HARDMEMLIMIT -p $QUEUE -A $ACCOUNT --time=${HARDTIMELIMIT} ${EXCLUSIVESTR} --output=/dev/null runcluster.sh
  else
      echo "input file "$SCIPPATH/$NETFILE" not found!"
  fi
done
