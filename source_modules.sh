#!/bin/bash
MODULEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/src
#export PYTHONPATH=$MODULEDIR:$PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:$MODULEDIR"
echo "Included $MODULEDIR in PYTHONPATH"
