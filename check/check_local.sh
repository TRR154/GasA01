#!/usr/bin/env bash
# filepath: check/check_local_parallel.sh
# Parallel version of check_local.sh with SCIP-style logging

TSTNAME=$1
BINNAME=$2
SETTINGS=$3
TIMELIMIT=$4
NODELIMIT=$5
MEMLIMIT=$6
DISPFREQ=$7
VERSION=$8
NCORES=$9   # Number of parallel jobs
LPS=${10}

if [ -z "$NCORES" ]; then
    NCORES=8
fi

mkdir -p results

# Big merged files
BIGOUT="results/check.$TSTNAME.$(basename $BINNAME).$SETTINGS.out"
BIGERR="results/check.$TSTNAME.$(basename $BINNAME).$SETTINGS.err"
rm -f "$BIGOUT" "$BIGERR"

TOTAL=$(wc -l < testset/$TSTNAME.test)

if [ ! -e "../settings/${SETTINGS}.set" ]; then
    echo "Skipping test since the settings file ../settings/${SETTINGS}.set does not exist."
    exit 1
fi

# keep both: full path (for SCIP) and short name (for filenames)
SETTINGS="../settings/${SETTINGS}.set"
SETNAME=$(basename "$SETTINGS" .set)

COUNT=0
while read -r BASENAME SCNFILE; do
    COUNT=$((COUNT+1))
    NETFILE="${BASENAME}.net"
    CSFILE="${BASENAME}.cs"

    (
        BASE_NET=$(basename "$BASENAME")
        BASE_SCN=$(basename "$SCNFILE" .scn)

        OUTFILE_INSTANCE="results/check.${TSTNAME}.$(basename $BINNAME).${SETNAME}.${BASE_NET}.${BASE_SCN}.out"
        ERRFILE_INSTANCE="results/check.${TSTNAME}.$(basename $BINNAME).${SETNAME}.${BASE_NET}.${BASE_SCN}.err"

        #echo "Creating logs: $OUTFILE_INSTANCE , $ERRFILE_INSTANCE"
        rm -f "$OUTFILE_INSTANCE" "$ERRFILE_INSTANCE"
        touch "$OUTFILE_INSTANCE" "$ERRFILE_INSTANCE"

        # Check if .net file exists
        if [ ! -f "$NETFILE" ]; then
            {
                echo "@02 FILE NOT FOUND: $NETFILE ==========="
            } >>"$ERRFILE_INSTANCE"
            exit 0
        fi

        # Build files argument
        if [ -f "$CSFILE" ]; then
            FILES="$NETFILE $SCNFILE -b $CSFILE"
        else
            FILES="$NETFILE $SCNFILE"
        fi
        echo "Starting: $FILES"

        # Log header to OUTFILE
        {
            echo "@01 $FILES ==========="
            echo "-----------------------------"
            date
            echo "-----------------------------"
            date +"@03 %s"
            echo " ulimit -f 200000; ../$BINNAME $FILES -s $SETTINGS -t $TIMELIMIT -mem $MEMLIMIT -d $DISPFREQ -n $NODELIMIT"
        } >>"$OUTFILE_INSTANCE"

        # Log header to ERRFILE
        {
            echo "@01 $FILES ==========="
            echo "-----------------------------"
            date
            echo "-----------------------------"
        } >>"$ERRFILE_INSTANCE"

        # Run SCIP
        bash -c "ulimit -f 200000; ../$BINNAME $FILES -s $SETTINGS -t $TIMELIMIT -n $NODELIMIT -mem $MEMLIMIT -d $DISPFREQ " \
            >>"$OUTFILE_INSTANCE" 2>>"$ERRFILE_INSTANCE"
        echo "Finished $COUNT out of $TOTAL: $FILES"

        # Log footer to OUTFILE
        {
            date +"@04 %s"
            echo "-----------------------------"
            date
            echo "-----------------------------"
            echo
            echo "=ready="
        } >>"$OUTFILE_INSTANCE"

        # Log footer to ERRFILE
        {
            date
            echo "-----------------------------"
        } >>"$ERRFILE_INSTANCE"

    ) &

    # Add PID to a job array
    pids+=($!)

# Wait if too many jobs
while [ "${#pids[@]}" -ge "$NCORES" ]; do
    # Remove finished PIDs from array
    for i in "${!pids[@]}"; do
        if ! kill -0 "${pids[i]}" 2>/dev/null; then
            unset 'pids[i]'
        fi
    done
    sleep 0.1
done

done < testset/$TSTNAME.test

wait
echo "All $COUNT jobs finished."

# Concatenate per-instance logs into one big file
while read -r BASENAME SCNFILE; do
    BASE_NET=$(basename "$BASENAME")
    BASE_SCN=$(basename "$SCNFILE" .scn)

    OUTFILE_INSTANCE="results/check.${TSTNAME}.$(basename $BINNAME).${SETNAME}.${BASE_NET}.${BASE_SCN}.out"
    ERRFILE_INSTANCE="results/check.${TSTNAME}.$(basename $BINNAME).${SETNAME}.${BASE_NET}.${BASE_SCN}.err"

    cat "$OUTFILE_INSTANCE" >>"$BIGOUT"
    cat "$ERRFILE_INSTANCE" >>"$BIGERR"
done < testset/$TSTNAME.test


# Automatically run evaluation
echo "Running evalcheck.sh on $BIGOUT ..."
./evalcheck.sh "$BIGOUT"
