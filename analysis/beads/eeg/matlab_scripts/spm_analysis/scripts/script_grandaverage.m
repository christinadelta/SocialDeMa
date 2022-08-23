spm('defaults', 'eeg');

S = [];
%%
S.D = [
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_01_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_02_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_03_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_04_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_05_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_06_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_07_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_08_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_09_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_11_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_12_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_13_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_14_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_15_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_16_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_17_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_18_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_19_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_20_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_21_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_22_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_23_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_24_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_25_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_26_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_27_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_28_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_29_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_30_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_31_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_32_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_33_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_34_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_35_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_36_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_37_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_38_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_39_beads_block_01.mat'
       '/Users/christinadelta/githubstuff/rhul_stuff/SocialDeMa/analysis/beads/eeg/matlab_scripts/spm_analysis/averages/maceerpfdfMspmeeg_sub_40_beads_block_01.mat'
       ];
%%
S.outfile = 'grand_average';
S.weighted = 1;
D = spm_eeg_grandmean(S);


