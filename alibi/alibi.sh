#!/bin/bash
while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo "Usage: $0 [options]"
      echo " "
      echo "-h, --help                          show brief help"
      echo "-i                                  specify an input file to use"
      echo "-o                                  specify a file to store output in"
      exit 0
      ;;
    -i)
      shift
      if test $# -gt 0; then
        export INPUT=$1
        echo "$INPUT"
      else
        echo "no input file specified"
        exit 1
      fi
      shift
      ;;
    -o)
      shift
      if test $# -gt 0; then
        export OUTPUT=$1
        echo "$OUTPUT"
      else
        echo "no output file specified"
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

if [ -z $INPUT ]; then
  echo "no input file specified"
  exit 1
fi

if [ -z $OUTPUT ]; then
  echo "no output file specified"
  exit 1
fi

python gfaSorter.py $INPUT $OUTPUT
