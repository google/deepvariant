#!/bin/bash

# set -e -o pipefail

if [ $# -eq 0 ]; then
    echo "Usage: $0 <SOURCE DIRECTORY> <WORKING DIRECTORY>"
    exit
fi

DIRPATH=$1
WORKDIR=$2

wait_for_eval_script_to_wait() {
    while true; do
        look_at_tail_0=`tail -1 $WORKDIR/eval.log | grep "Terminating eval after 10000 seconds of no checkpoints" | wc -l`
        look_at_tail_1=`tail -1 $WORKDIR/eval.log | grep "Waiting for new checkpoint at" | wc -l`
        look_at_tail=$((look_at_tail_0 + look_at_tail_1))

        if [ $look_at_tail -eq 1 ]; then
            echo "Found eval script waiting for the next checkpoint"
            sleep 1
            break
        else
            echo "Eval script is running. Waiting ... "
            sleep 15
        fi
    done
}

add_checkpoint() {
    checkpoint_name=$1
    echo "all_model_checkpoint_paths: \"$checkpoint_name\"" >> $WORKDIR/checkpoint.so_far
    printf "model_checkpoint_path: \"$checkpoint_name\"\n$(cat $WORKDIR/checkpoint.so_far)" > $WORKDIR/checkpoint
}

steps=$(for model in `ls $DIRPATH/model.ckpt-[0-9]*.index`; do
    model_step=`echo $model | sed 's?.*model.ckpt-\([0-9]*\).index?\1?g'`
    echo $model_step
done | sort -n)

echo "Found the following checkpoint steps: $steps"

for model_step in $steps; do
    echo "Running evaluation for step $model_step"

    # Wait for any previous checkpoint evaluations to complete
    wait_for_eval_script_to_wait

    # Copy checkpoint files
    cp $DIRPATH/model.ckpt-${model_step}.data* $WORKDIR/.
    cp $DIRPATH/model.ckpt-${model_step}.meta $WORKDIR/.
    cp $DIRPATH/model.ckpt-${model_step}.index $WORKDIR/.

    # Add checkpoint
    add_checkpoint model.ckpt-${model_step}

    # Now simply wait for checkpoint evaluation to continue
    echo "Copied checkpoint. Waiting for evaluation on checkpoint to complete"
    sleep 10
done

echo "Evaluations completed. Results in $WORKDIR"
