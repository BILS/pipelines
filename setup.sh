#!/bin/bash

# Script to help with setting up the BILS Bpipe pipeline infrastructure

WORKING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

MODULE_DIR="$WORKING_DIR/modules"
BPIPE_FOLDER=( $WORKING_DIR/bin/bpipe* )
BPIPE_EXE=$( which bpipe 2>/dev/null )

if [[ -n "$BPIPE_LIB" ]]; then
    echo "\$BPIPE_LIB is already set to '$BPIPE_LIB' (good)"
else
    echo "\$BPIPE_LIB not set, add this to ~/.bash_profile:"
    echo "  export BPIPE_LIB=$MODULE_DIR"
fi

if [[ -n "$BPIPE_EXE" ]]; then
    echo "Bpipe executable already installed as '$BPIPE_EXE' (good)"
else
    echo "Bpipe executable not found in \$PATH, add this to ~/.bash_profile:"
    echo "  export PATH=\$PATH:$BPIPE_FOLDER/bin"
fi

if [[ -n "$BPIPE_BIN" ]]; then
	echo "\$BPIPE_BIN is already set to '$BPIPE_BIN' (good)"
else
	echo "\$BPIPE_BIN not set, add this to ~/.bash_profile:"
	echo " export BPIPE_BIN=$WORKING_DIR/bin"
fi
