#!/bin/bash

# created in October 2022 @christinadelta
# part of the Beads EEG analysis pipeline

# this simple script first creates directories (if not already created) for contrast images
# of Beads EEG data. Because the name of each image spm creates is the same for each subject
# the code will rename (add a subject number to each img file) and move the contrast images
# to the corresponding folder/dir

# --------------------------------

# first define paths and variables
basedir="/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/"
bashdir=$basedir$"shell_scripts"
spmfold="matlab_scripts/spm_analysis/"
spmdir=$basedir$spmfold
echo $spmdir # check that dirs were properly created
outputdir=$spmdir$"output/"
erpc=$spmdir$"erp_contrasts/"
tfrc=$spmdir$"tfr_contrasts/"
numsubs=40

# if erp and tfr contrasts dirs do not exist, create it
mkdir -p $erpc $tfrc

# create subdirs in the contrast directories (if they don't exist)
mkdir -p $erpc$"diffEasy" $erpc$"urnDraw" $erpc$"inter" $erpc$"urns" $erpc$"draws"
mkdir -p $tfrc$"diffEasy" $tfrc$"urnDraw" $tfrc$"inter" $tfrc$"urns" $tfrc$"draws"

# loop over subjects and start moving and renaming the contrast images
for i in {1..2}
do # go to this_sub dir grab the image files from the folder and store them in the corresponding constrast folder

  # 1. first go to the erp contrats and then to tfr contrasts
  # 1.a difficult vs easy contrast images
  # 1.b urn vs darw contrast images
  # 1.c interaction images
  # 1.d urns images
  # 1.e draws images

  # if subject num is smaller than 10 add a decimal
  if (( ${i} < 10 )); then
    subout=$outputdir$"sub-0${i}/"
    tmp="sub-0${i}"

  else
    subout=$outputdir$"sub-${i}/"
    tmp="sub-${i}"
  fi # end of if statement

  echo $tmp
done # end of for loop
