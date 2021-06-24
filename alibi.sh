#!/bin/bash

while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo "Usage: $0 [options]"
      echo " "
      echo "-h, --help                          show brief help"
      echo "-i                                  specify an input file to use"
      exit 0
      ;;
    -i)
      shift
      if test $# -gt 0; then
        export INPUT=$1
      else
        echo "no input file specified"
        exit 1
      fi
      shift
      ;;
    *)
      echo "$1 is not a recognized flag!"
      exit 1
      ;;
  esac
done

python ./bin/gfaSorter.py $INPUT
