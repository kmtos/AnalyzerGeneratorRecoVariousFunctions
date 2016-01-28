To RECO:

	cd RECO/

Then, you must edit line 32 of STEP_1_cfg.py, line 31 of STEP_2_cfg.py, and line 27 of STEP_3_cfg.py to the location of your input files. 
On line 8 of STTEP_3_cfg.py, change this to anything you want. The typical string is "PAT", but I needed to change it so I could use "PAT" for another producer plot. If you plan on using a producer after RECO step3, then leave it as PAT1. 
Also, on lines 60-65, put any collections you would like to keep there from the RECO step2.

Next, edit the paths in the SUBMIT_.sh files.

You shouldn't have to edit RECO.sh. This is the main script you will use. Sometimes, if a lot of the jobs fail, what you can do is get a list of the numbers that failed and put them in the variable "files" in RECO_fix.sh. 
I don't usually type them in by hand but use the bash shell command "diff" and some nice vi tricks to then copy into that variable.

Then you should be ready to go. Execute this command:
	./RECO.sh STEP_1_cfg SUBMIT_1 <Name of your choice> <queue of your choice> <Starting number for files>  <number of files generated in the gen sample>

The RECO_fix.sh is run the same way, but the number of files generated does absolutely nothing.

