#!/bin/sh

#echo "./clean_data.sh"
#./clean_data.sh

echo "nohup ./espic.x >> screenout 2>&1 &"
../bin/espic.x
