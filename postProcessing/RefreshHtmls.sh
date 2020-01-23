#!/bin/bash

#https://unix.stackexchange.com/questions/187167/traverse-all-subdirectories-in-and-do-something-in-unix-shell-script


#If multiple notebook changed and you are on Unix
#run this script and it updates the notebooks using a tree walker
STR = /home/cneijssel/Documents/COMPASpop/popsynth/Papers/NeijsselEtAL/PostProcessing
echo $d
for d in $(find $STR -maxdepth 1 -type d)
do
  #Do something, the directory is accessible with $d:
  echo $d
  ipython nbconvert *.ipynb
done >output_file
