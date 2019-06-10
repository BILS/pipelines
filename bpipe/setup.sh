#!/bin/bash

# Script to help with setting up the BILS Bpipe pipeline infrastructure

WORKING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

MODULE_DIR="$WORKING_DIR/modules"
BPIPE_FOLDER=( $WORKING_DIR/bin/bpipe* )
BPIPE_EXE=$( which bpipe 2>/dev/null )

if [[ -n "$BPIPE_LIB" ]]; then
    echo -e "(\033[32mgood\033[0m) \$BPIPE_LIB is already set to '$BPIPE_LIB'"
else
    echo -e "(\033[31mbad\033[0m) \$BPIPE_LIB not set, add this to ~/.bash_profile:"
    echo "  export BPIPE_LIB=$MODULE_DIR"
fi

if [[ -n "$BPIPE_EXE" ]]; then
    echo -e "(\033[32mgood\033[0m) Bpipe executable already installed as '$BPIPE_EXE'"
else
    echo -e "(\033[31mbad\033[0m) Bpipe executable not found in \$PATH, add this to ~/.bash_profile:"
    echo "export PATH=\$PATH:$BPIPE_FOLDER/bin"
fi

if [[ -n "$BPIPE_BIN" ]]; then
	echo -e "(\033[32mgood\033[0m) \$BPIPE_BIN is already set to '$BPIPE_BIN'"
else
	echo -e "(\033[31mbad\033[0m) \$BPIPE_BIN not set, add this to ~/.bash_profile:"
	echo " export BPIPE_BIN=$WORKING_DIR/bin"
fi
if [[ -n "$BPIPE_CONFIG" ]]; then
        echo -e "(\033[32mgood\033[0m) \$BPIPE_CONFIG is already set to '$BPIPE_CONFIG'"
else
        echo -e "(\033[31mbad\033[0m) \$BPIPE_CONFIG not set, add this to ~/.bash_profile:"
        echo " export BPIPE_CONFIG=$WORKING_DIR/config"
fi
if [[ -n "$BPIPE_PIPELINE" ]]; then
        echo -e "(\033[32mgood\033[0m) \$BPIPE_PIPELINE is already set to '$BPIPE_PIPELINE'"
else
        echo -e "(\033[31mbad\033[0m) \$BPIPE_PIPELINE not set, add this to ~/.bash_profile:"
        echo " export BPIPE_PIPELINE=$WORKING_DIR/pipelines"
fi
