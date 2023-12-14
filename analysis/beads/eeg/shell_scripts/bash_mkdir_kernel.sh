#!/bin/bash

# created in December 2023 by @christinadelta
# part of the Beads EEG analysis pipeline

# this simple .sh script creates directries and sub directories 
# to store smoothed TFR contrasts for each smoothing kernel that we want to try


# --------------------------------
# first define paths and variables
basedir="/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/"
datadir="/Volumes/DeepSpaceStuff/optimal_stopping_data/beads/"
bashdir=$basedir$"shell_scripts"
spmfold="spm_analysis/"
spmdir=$datadir$spmfold
echo $spmdir # check that dirs were properly created
tfra=$spmdir$"tfr_contr_smoothed_alpha/"
tfrb=$spmdir$"tfr_contr_smoothed_beta/"
tfrg=$spmdir$"tfr_contr_smoothed_gamma/"
tfrg_sp=$spmdir$"tfr_contr_smoothed_gamma_sp/"

# if if the tfr contrasts dirs do not exist, create them:
mkdir -p $tfra $tfrb $tfrg $tfrg_sp

###### START WITH ALPHA TF --- CREATE DIRS FOR EACH SMOOTHING KERNEL
# create subdirectories within each of the tfr directories to store the kernel-specific smoothed images
# KERNEL 1 -- default [8 8 8]
# create subdirs in the contrast directories (if they don't exist)
krnl_one=$tfra$"kernel_one/" 
mkdir -p $krnl_one
mkdir -p $krnl_one$"diffEasy" $krnl_one$"urnDraw" $krnl_one$"inter" $krnl_one$"urns" $krnl_one$"draws"

# KERNEL 2 [4 4 20]
# create subdirs in the contrast directories (if they don't exist)
krnl_two=$tfra$"kernel_two/" 
mkdir -p $krnl_two
mkdir -p $krnl_two$"diffEasy" $krnl_two$"urnDraw" $krnl_two$"inter" $krnl_two$"urns" $krnl_two$"draws"

# KERNEL 3 [4 8 30]
krnl_three=$tfra$"kernel_three/" 
mkdir -p $krnl_three
mkdir -p $krnl_three$"diffEasy" $krnl_three$"urnDraw" $krnl_three$"inter" $krnl_three$"urns" $krnl_three$"draws"

# KERNEL 4 [4 8 10]
krnl_four=$tfra$"kernel_four/" 
mkdir -p $krnl_four
mkdir -p $krnl_four$"diffEasy" $krnl_four$"urnDraw" $krnl_four$"inter" $krnl_four$"urns" $krnl_four$"draws"

# KERNEL 5 [8 10 15]
krnl_five=$tfra$"kernel_five/" 
mkdir -p $krnl_five
mkdir -p $krnl_five$"diffEasy" $krnl_five$"urnDraw" $krnl_five$"inter" $krnl_five$"urns" $krnl_five$"draws"


# ------------------------------------------ # 
###### BETA TF --- CREATE DIRS FOR EACH SMOOTHING KERNEL
# create subdirectories within each of the tfr directories to store the kernel-specific smoothed images
# KERNEL 1 -- default [8 8 8]
# create subdirs in the contrast directories (if they don't exist)
krnl_one=$tfrb$"kernel_one/" 
mkdir -p $krnl_one
mkdir -p $krnl_one$"diffEasy" $krnl_one$"urnDraw" $krnl_one$"inter" $krnl_one$"urns" $krnl_one$"draws"

# KERNEL 2 [4 4 20]
# create subdirs in the contrast directories (if they don't exist)
krnl_two=$tfrb$"kernel_two/" 
mkdir -p $krnl_two
mkdir -p $krnl_two$"diffEasy" $krnl_two$"urnDraw" $krnl_two$"inter" $krnl_two$"urns" $krnl_two$"draws"

# KERNEL 3 [4 8 30]
krnl_three=$tfrb$"kernel_three/" 
mkdir -p $krnl_three
mkdir -p $krnl_three$"diffEasy" $krnl_three$"urnDraw" $krnl_three$"inter" $krnl_three$"urns" $krnl_three$"draws"

# KERNEL 4 [4 8 10]
krnl_four=$tfrb$"kernel_four/" 
mkdir -p $krnl_four
mkdir -p $krnl_four$"diffEasy" $krnl_four$"urnDraw" $krnl_four$"inter" $krnl_four$"urns" $krnl_four$"draws"

# KERNEL 5 [8 10 15]
krnl_five=$tfrb$"kernel_five/" 
mkdir -p $krnl_five
mkdir -p $krnl_five$"diffEasy" $krnl_five$"urnDraw" $krnl_five$"inter" $krnl_five$"urns" $krnl_five$"draws"

# ------------------------------------------ # 
###### GAMMA ALL TF --- CREATE DIRS FOR EACH SMOOTHING KERNEL
# create subdirectories within each of the tfr directories to store the kernel-specific smoothed images
# KERNEL 1 -- default [8 8 8]
# create subdirs in the contrast directories (if they don't exist)
krnl_one=$tfrg$"kernel_one/" 
mkdir -p $krnl_one
mkdir -p $krnl_one$"diffEasy" $krnl_one$"urnDraw" $krnl_one$"inter" $krnl_one$"urns" $krnl_one$"draws"

# KERNEL 2 [4 4 20]
# create subdirs in the contrast directories (if they don't exist)
krnl_two=$tfrg$"kernel_two/" 
mkdir -p $krnl_two
mkdir -p $krnl_two$"diffEasy" $krnl_two$"urnDraw" $krnl_two$"inter" $krnl_two$"urns" $krnl_two$"draws"

# KERNEL 3 [4 8 30]
krnl_three=$tfrg$"kernel_three/" 
mkdir -p $krnl_three
mkdir -p $krnl_three$"diffEasy" $krnl_three$"urnDraw" $krnl_three$"inter" $krnl_three$"urns" $krnl_three$"draws"

# KERNEL 4 [4 8 10]
krnl_four=$tfrg$"kernel_four/" 
mkdir -p $krnl_four
mkdir -p $krnl_four$"diffEasy" $krnl_four$"urnDraw" $krnl_four$"inter" $krnl_four$"urns" $krnl_four$"draws"

# KERNEL 5 [8 10 15]
krnl_five=$tfrg$"kernel_five/" 
mkdir -p $krnl_five
mkdir -p $krnl_five$"diffEasy" $krnl_five$"urnDraw" $krnl_five$"inter" $krnl_five$"urns" $krnl_five$"draws"


# ------------------------------------------ # 
###### GAMMA range (38:42 Hz) TF --- CREATE DIRS FOR EACH SMOOTHING KERNEL
# create subdirectories within each of the tfr directories to store the kernel-specific smoothed images
# KERNEL 1 -- default [8 8 8]
# create subdirs in the contrast directories (if they don't exist)
krnl_one=$tfrg_sp$"kernel_one/" 
mkdir -p $krnl_one
mkdir -p $krnl_one$"diffEasy" $krnl_one$"urnDraw" $krnl_one$"inter" $krnl_one$"urns" $krnl_one$"draws"

# KERNEL 2 [4 4 20]
# create subdirs in the contrast directories (if they don't exist)
krnl_two=$tfrg_sp$"kernel_two/" 
mkdir -p $krnl_two
mkdir -p $krnl_two$"diffEasy" $krnl_two$"urnDraw" $krnl_two$"inter" $krnl_two$"urns" $krnl_two$"draws"

# KERNEL 3 [4 8 30]
krnl_three=$tfrg_sp$"kernel_three/" 
mkdir -p $krnl_three
mkdir -p $krnl_three$"diffEasy" $krnl_three$"urnDraw" $krnl_three$"inter" $krnl_three$"urns" $krnl_three$"draws"

# KERNEL 4 [4 8 10]
krnl_four=$tfrg_sp$"kernel_four/" 
mkdir -p $krnl_four
mkdir -p $krnl_four$"diffEasy" $krnl_four$"urnDraw" $krnl_four$"inter" $krnl_four$"urns" $krnl_four$"draws"

# KERNEL 5 [8 10 15]
krnl_five=$tfrg_sp$"kernel_five/" 
mkdir -p $krnl_five
mkdir -p $krnl_five$"diffEasy" $krnl_five$"urnDraw" $krnl_five$"inter" $krnl_five$"urns" $krnl_five$"draws"