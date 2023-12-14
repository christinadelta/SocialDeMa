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
datadir="/Volumes/DeepSpaceStuff/optimal_stopping_data/beads/"
bashdir=$basedir$"shell_scripts"
spmfold="spm_analysis/"
spmdir=$datadir$spmfold
echo $spmdir # check that dirs were properly created
outputdir=$spmdir$"output/"
# erpc=$spmdir$"erp_contrasts/"
tfrc=$spmdir$"tfr_contrasts_gamma_sp/"
numsubs=40

# if erp and tfr contrasts dirs do not exist, create it
# mkdir -p $erpc $tfrc
mkdir -p $tfrc

# create subdirs in the contrast directories (if they don't exist)
# mkdir -p $erpc$"diffEasy" $erpc$"urnDraw" $erpc$"inter" $erpc$"urns" $erpc$"draws"
mkdir -p $tfrc$"diffEasy" $tfrc$"urnDraw" $tfrc$"inter" $tfrc$"urns" $tfrc$"draws"

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

  '''
  # first move and rename (add sub number) the erp contrast images
  tmpdir=$subout$"wde_maceerpfdfMspmeeg_${tmp}beads_block_01" # diff vs easy contrast dir
  tmpdir2=$subout$"wud_maceerpfdfMspmeeg_${tmp}beads_block_01" # urn vs draw contrast dir
  tmpdir3=$subout$"wi_maceerpfdfMspmeeg_${tmp}beads_block_01" # uinteraction contrast dir
  tmpdir4=$subout$"wu_maceerpfdfMspmeeg_${tmp}beads_block_01" # urns contrast dir
  tmpdir5=$subout$"wd_maceerpfdfMspmeeg_${tmp}beads_block_01" # draws contrast dir

  # where to move?
  mv $tmpdir$"/condition_DiffVsEasy.nii" $erpc$"diffEasy/${tmp}condition_diffeasy.nii" # diff vs easy
  mv $tmpdir2$"/condition_urnVSdraw.nii" $erpc$"urnDraw/${tmp}condition_urndraw.nii" # urns vs draws
  mv $tmpdir3$"/condition_interaction.nii" $erpc$"inter/${tmp}condition_interaction.nii" # interaction
  mv $tmpdir4$"/condition_onlyurn.nii" $erpc$"urns/${tmp}condition_onlyurn.nii" # only urns
  mv $tmpdir5$"/condition_onlydraw.nii" $erpc$"draws/${tmp}condition_onlydraw.nii" # only draws '''

  # second move and rename (add sub number) the tfr contrast images
  tmpdir=$subout$"wde_mPgamma_sp_rtf_cetfrfdfMspmeeg_${tmp}beads_block_01" # diff vs easy contrast dir
  tmpdir2=$subout$"wud_mPgamma_sp_rtf_cetfrfdfMspmeeg_${tmp}beads_block_01" # urn vs draw contrast dir
  tmpdir3=$subout$"wi_mPgamma_sp_rtf_cetfrfdfMspmeeg_${tmp}beads_block_01" # uinteraction contrast dir
  tmpdir4=$subout$"wu_mPgamma_sp_rtf_cetfrfdfMspmeeg_${tmp}beads_block_01" # urns contrast dir
  tmpdir5=$subout$"wd_mPgamma_sp_rtf_cetfrfdfMspmeeg_${tmp}beads_block_01" # draws contrast dir

  # where to move?
  mv $tmpdir$"/condition_DiffVsEasy.nii" $tfrc$"diffEasy/${tmp}condition_diffeasy.nii" # diff vs easy
  mv $tmpdir2$"/condition_urnVSdraw.nii" $tfrc$"urnDraw/${tmp}condition_urndraw.nii" # urns vs draws
  mv $tmpdir3$"/condition_interaction.nii" $tfrc$"inter/${tmp}condition_interaction.nii" # interaction
  mv $tmpdir4$"/condition_onlyurn.nii" $tfrc$"urns/${tmp}condition_onlyurn.nii" # only urns
  mv $tmpdir5$"/condition_onlydraw.nii" $tfrc$"draws/${tmp}condition_onlydraw.nii" # only draws

done # end of for loop
