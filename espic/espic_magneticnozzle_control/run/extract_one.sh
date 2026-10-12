#!/bin/bash

LINES="170"
LINE2TI="100"

for line in ${LINES}
do
  ti=$(($line * $LINE2TI))
  echo "sed \"${line}q;d\" hist_phi.dat > one_phi_te${ti}.dat"
  sed "${line}q;d" hist_phi.dat > one_phi_te${ti}.dat
  echo "sed \"${line}q;d\" hist_chgDenst.dat > one_chgDenst_te${ti}.dat"
  sed "${line}q;d" hist_chgDenst.dat > one_chgDenst_te${ti}.dat
  echo "sed \"${line}q;d\" hist_velField.dat > one_velField_te${ti}.dat"
  sed "${line}q;d" hist_velField.dat > one_velField_te${ti}.dat
  echo "sed \"${line}q;d\" new_hist_velField.dat > new_one_velField_te${ti}.dat"
  sed "${line}q;d" new_hist_velField.dat > new_one_velField_te${ti}.dat
   echo "sed \"${line}q;d\" new_new_hist_velField.dat > new_new_one_velField_te${ti}.dat"
  sed "${line}q;d" new_new_hist_velField.dat > new_new_one_velField_te${ti}.dat
done
