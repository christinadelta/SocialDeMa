# Beads-task EEG preprocessing/analysis pipeline

Formal preprocessing and analyses of the EEG data (Beads task). 

* Conversion to spm object
* Define channels (EEG + EOG)
* Look for bad channels – this is an important step and needs to be done for each subject before creating the montage. If a subject has a bad channel during recording, we need to specify it at the next step and exclude it from re-referencing.  
* Create montage (EEG + vertical, horizontal EOG) & re-reference (average for now) – use modified montage.m function
* High pass filtering: cutoff 0.5 (or 0.1)
* Downsampling: 256Hz
* Low pass filtering: cutoff 30
* Epoching (with baseline correction) – use modified trialdef.m function
	* Probability: easy/diff
	* Choices: draw/urn
* Merging (concatenating) blocks 
* Preparing – Channel coordinates for fiducials: will use the default spm Biosemi64 montage
* Artefact detection/removal
* **Evoked analysis (ERPs)**:
	* Average epochs/trials:
	* Compute Contrast conditions (use the average – contrast method):
		* Urn vs Draw [-1 1 -1 1]
		* Difficult vs Easy [1 1 -1 -1]
		* Interaction [-1 1 1 -1]
		* Only urns [0 1 0 1]
		* Only draws [1 0 1 0]
	* Convert to 3D images 
* 2nd level Mass Univariate analysis using the contrasts [urn vs draw, diff vs easy, interaction] and the entire peristimulus time [-500 800]
* Individual differences analysis:
	* 2nd-level Massive Univariate analysis – Use urn choices contrasts and array of number of draws as covariate 
* Association of evoked responses and action values 
	* Compute grand average (over subjects). Using the grand average file, select a specific time-range and sensors that show larger response to urn choices and draw choices 
	* For each participant, crop the urn choice trials around the time-range and for sensors specified at the above step and extract the data from the meeg object as a new dataset. 
	* Compute the difference between action values (dAQ) for drawing again and choosing urns
	* Linear Regression (draw-by draw/ 1st level/ within subjects) – using the diff in (dAQ) as regressor
	* Use the 1st level beta weights to model a 2nd level (one sample) t-test 
	* Optional: draw-by-draw (epochs) correlation with dAQ → get one r value per participant → convert r values to z (Fisher’s r-to-z transformation) and run (one sample t-test  at group level (to test if overall r values differ from 0). 
* **Time-Frequency Representation Analysis (TFRs)**: 
	* Wavelet estimation (factor 7) for power 1-55Hz – output is 2 files (phase and power)
	* Baseline rescaling (-500ms to -50ms) – only power file
	* Average epochs/trials:
		* Over beta range [13 to 30 Hz] 
		* Over time using the entire peristimulus time
	* Traditional robust average of conditions 
	* Compute Contrast conditions (use the average – contrast method) 
		* Urn vs Draw [-1 1 -1 1]
		* Difficult vs Easy [1 1 -1 -1]
		* Interaction [-1 1 1 -1]
		* Only urns [0 1 0 1]
		* Only draws [1 0 1 0]
	* Convert to 3D images 
* 2nd level Mass Univariate analysis using the contrasts [urn vs draw, diff vs easy, interaction] and the entire peristimulus time [-500 800] and all frequencies
* Individual differences analysis:
	* 2nd-level Massive Univariate analysis – Use urn choices contrasts and number of draws as covariate 
* Association of beta power and action values: 
	* Average TFR-rescaled object over specific time and frequency range:
		* Compute grand average (over subjects). Using the grand average file, select a specific time-range, frequency (beta range) and sensors that show larger response to urn choices and draw choices
		* Average subject meeg object over frequency (beta range – slow [13-20Hz] and fast [21-30Hz])
		* Average subject meeg object over time selected using the GrandMean object 
	* Compute the difference between action values (dAQ) for drawing again and choosing urns
	* Linear Regression (draw-by draw/ 1st level/ within subjects) with beta range data as dependent variable – using the diff in (dAQ) as regressor
	* Use the 1st level beta weights to model a 2nd level (one sample) t-test 
	* Optional: draw-by-draw (epochs in beta range) correlation with dAQ → get one r value per participant → convert r values to z (Fisher’s r-to-z transformation) and run (one sample t-test  at group level (to test if overall r values differ from 0). 











	

