#!/bin/bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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

SETDIR=../settings

# get current SCIP path
SCIPPATH=`pwd`

if test ! -e results
then
    mkdir results
fi

OUTFILE=results/check.$TSTNAME.$BINID.$SETNAME.out
ERRFILE=results/check.$TSTNAME.$BINID.$SETNAME.err
RESFILE=results/check.$TSTNAME.$BINID.$SETNAME.res
TEXFILE=results/check.$TSTNAME.$BINID.$SETNAME.tex
SETFILE=results/check.$TSTNAME.$BINID.$SETNAME.set

MVORCP=mv

DATEINT=`date +"%s"`
if test -e $OUTFILE
then
    $MVORCP $OUTFILE $OUTFILE.old-$DATEINT
fi
if test -e $ERRFILE
then
    $MVORCP $ERRFILE $ERRFILE.old-$DATEINT
fi

uname -a >>$OUTFILE
uname -a >>$ERRFILE
date >>$OUTFILE
date >>$ERRFILE

# unset settings name if necessary
if test "$SETNAME" = "default"
then
    SETNAME=""
else
    if test ! -e $SETDIR/$SETNAME.set
    then
        echo "Skipping test since the settings file $SETDIR/${SETNAME}.set does not exist."
        exit 1
    fi

    tmp=$SETNAME
    SETNAME="-s "$SETDIR"/"$tmp".set"
fi

# process test instances
STORENAME=""
for i in `cat testset/$TSTNAME.test` DONE
do
    if test "$i" = "DONE"
    then
	break
    fi

    NETFILE=$STORENAME".net"
    CSFILE=$STORENAME".cs"
    if test "$STORENAME" = ""
    then
	STORENAME="$i"
	continue
    else
	STORENAME=""
    fi

    # determine filenames
    SCNFILE=$i

    # check if problem instance exists
    if test -f $SCIPPATH/$NETFILE
    then
	if test -f $SCIPPATH/$CSFILE
	then
	    FILES=$SCIPPATH/$NETFILE" "$SCIPPATH/$SCNFILE" -b "$SCIPPATH/$CSFILE
	else
	    FILES=$SCIPPATH/$NETFILE" "$SCIPPATH/$SCNFILE
	fi

	echo @01 $FILES ===========
	echo @01 $FILES ===========                >> $ERRFILE
	echo -----------------------------
	date
	date >>$ERRFILE
	echo -----------------------------
	date +"@03 %s"
	echo " ulimit -f 200000; ../$BINNAME $FILES $SETNAME -t $TIMELIMIT -mem $MEMLIMIT -n $NODELIMIT -d $DISPFREQ"
	bash -c " ulimit -f 200000; ../$BINNAME $FILES $SETNAME -t $TIMELIMIT -mem $MEMLIMIT -n $NODELIMIT -d $DISPFREQ" 2>>$ERRFILE
	date +"@04 %s"
	echo -----------------------------
	date
	date >>$ERRFILE
	echo -----------------------------
	echo
	echo =ready=
    else
	echo @02 FILE NOT FOUND: $NETFILE ===========
	echo @02 FILE NOT FOUND: $NETFILE =========== >>$ERRFILE
    fi
done | tee -a $OUTFILE

date >>$OUTFILE
date >>$ERRFILE

if test -e $DONEFILE
then
    ./evalcheck.sh $OUTFILE
fi
