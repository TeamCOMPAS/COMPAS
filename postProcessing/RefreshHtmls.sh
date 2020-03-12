#!/bin/bash

#https://unix.stackexchange.com/questions/187167/traverse-all-subdirectories-in-and-do-something-in-unix-shell-script


#If multiple notebook changed and you are on Unix
#run this script and it updates the notebooks using a tree walker
STR = ./
echo $d
for d in $(find $STR -maxdepth 2 -type d)
do
  #Do something, the directory is accessible with $d:
  echo $d"/"
  #count the ipynb files
  #https://stackoverflow.com/questions/3856747/check-whether-a-certain-file-type-extension-exists-in-directory
  count=`ls -1 $d"/"*.ipynb 2>/dev/null | wc -l`
  echo $count
  if [ $count != "0" ];
    then
      ipython nbconvert  $d"/"*.ipynb 
    else
      # no files in directory so dont try to convert
      :
  fi
done


