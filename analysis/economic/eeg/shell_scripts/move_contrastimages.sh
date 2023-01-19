#!/bin/bash

# created in January 2023 @christinadelta
# part of the OPTIMAL STOPPING EEG analysis pipeline

# this simple script first creates directories (if not already created) for contrast images
# of economic and face EEG data. Because the name of each image spm creates is the same for each subject
# the code will rename (add a subject number to each img file) and move the contrast images
# to the corresponding folder/dir

# --------------------------------

# first define paths and variables
basedir="/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/economic/eeg/"
bashdir=$basedir$"shell_scripts"
datadir="/Volumes/DeepSpaceStuff/optimal_stopping_data/economic/"
spmfold="spm_analysis/"
spmdir=$datadir$spmfold
echo $spmdir # check that dirs were properly created
outputdir=$spmdir$"output/"
erpc=$spmdir$"erp_contrasts/"
tfrc=$spmdir$"tfr_contrasts/"
numsubs=40

# if erp and tfr contrasts dirs do not exist, create it
mkdir -p $erpc $tfrc

# create subdirs in the contrast directories (if they don't exist)
mkdir -p $erpc$"acceptReject" $erpc$"reject" $erpc$"accept"
mkdir -p $tfrc$"acceptReject" $tfrc$"reject" $tfrc$"accept"

# loop over subjects and start moving and renaming the contrast images
for i in {1..40}
do # go to this_sub dir grab the image files from the folder and store them in the corresponding constrast folder
# 1. first go to the erp contrats and then to tfr contrasts
  # if subject num is smaller than 10 add a decimal
  if (( ${i} < 10 )); then
    subout=$outputdir$"sub-0${i}/"
    tmp="sub_0${i}_" # grab subject number

  else
    subout=$outputdir$"sub-${i}/"
    tmp="sub_${i}_"
  fi # end of if statement
  # test that you are in the correct dir and with the correct sub number
  echo $subout
  echo $tmp

  # first move and rename (add sub number) the erp contrast images
  tmpdir=$subout$"wra_maceerpfdfMspmeeg_${tmp}economic_block_01" # accept vs reject contrast dir
  tmpdir2=$subout$"wr_maceerpfdfMspmeeg_${tmp}economic_block_01" # only reject contrast dir
  tmpdir3=$subout$"wa_maceerpfdfMspmeeg_${tmp}economic_block_01" # only accept contrast dir

  # where to move?
  mv $tmpdir$"/condition_rejectVSaccept.nii" $erpc$"acceptReject/${tmp}condition_rejectVSaccept.nii" # accept vs reject contrast dir
  mv $tmpdir2$"/condition_reject.nii" $erpc$"reject/${tmp}condition_reject.nii" # only reject contrast dir
  mv $tmpdir3$"/condition_accept.nii" $erpc$"accept/${tmp}condition_accept.nii" # only accept contrast dir

  # second move and rename (add sub number) the tfr contrast images
  tmpdir=$subout$"wra_mPrtf_cetfrfdfMspmeeg_${tmp}economic_block_01" # accept vs reject contrast dir
  tmpdir2=$subout$"wr_mPrtf_cetfrfdfMspmeeg_${tmp}economic_block_01" # only reject contrast dir
  tmpdir3=$subout$"wa_mPrtf_cetfrfdfMspmeeg_${tmp}economic_block_01" # only accept contrast dir

  # where to move?
  mv $tmpdir$"/condition_rejectVSaccept.nii" $tfrc$"acceptReject/${tmp}condition_rejectVSaccept.nii" # accept vs reject contrast dir
  mv $tmpdir2$"/condition_reject.nii" $tfrc$"reject/${tmp}condition_reject.nii" # only reject contrast dir
  mv $tmpdir3$"/condition_accept.nii" $tfrc$"accept/${tmp}condition_accept.nii" # only accept contrast dir

done # end of for loop
