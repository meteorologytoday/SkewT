#!/bin/bash
set -x
echo "Combine input file: $1"
echo "Combine output file: $2"

convert -gravity SouthWest \( "$1" -repage 0x0+0+0 \) \( logo.png -resize 30% \) -geometry +800+800 -composite "$2"
