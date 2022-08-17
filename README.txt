These are the MSRESOLVE Workshop materials, put together for a workshop in August 2022.
At the end of this workshop, participants will know how to use MSRESOLVE in a basic way to get concentrations from mass spectrometry data,
how to use several of the features, as well as where to find the settings of the many features.
Different directories exist and will also be made by the user for different examples of tasks that can be performed with MSRESOLVE.


PARTICIPANTS ARE EXPECTED TO HAVE ANACONDA INSTALLED AND TO HAVE TESTED THAT MSRESOLVE RUNS PRIOR TO THIS WORKSHOP.
IF THEY HAVE NOT ALREADY DONE SO, PARTICIPANTS SHOULD FOLLOW THE INSTRUCTIONS IN 000-ParticpantWorkshopPreparationInstructions.docx

If you have not already done so, install the required dependencies using the following command from spyder or an anaconda prompt: pip install -r requirements.txt
If pip does not work, first install pip and spyder from the anaconda interface or typing in "conda install pip" and "conda install spyder"
Python files can be run from an anaconda prompt by typing "python runfile.py" , or can be run by pressing the play button after opening in spyder.
From spyder, you should restart the kernel between runs.

00 An introductory slide.

   Before working on the exercises, we now need to get an anaconda window or spyder window open.
   

01  "Run" a premade typical analysis: 
    Navigate into that directory and run the runfile.py.  Various outputs will come out in that directory.
    ScaledConcentrations.tsv  and the graphs directory are the main outputs.  Section 2.3 of the Quickstart explains what the various files are.
    Right click on ScaledConcentrations.tsv, click on properties, and click "change" to change the default program for opening the tsv to be excel. Note: you may need to find out where excel is on your computer.
    TSV files can also be dragged into excel when excel is already open.
    As noted in Section 3.1 of the Quickstart document, for this analysis there is no 'detectable' Ethene or Actylene (those are 'within error'). There is only Ethane. The abilty to have error bars are an important feature of MSRESOLVE

02 Creating a reference pattern file.
   This section is OPTIONAL for MSRESOLVE, but is useful, as it describes how to retrieve reference patterns from online by use of JDXConverter (another program).
   In this directory and your anaconda prompt, use 'pip install -r requirements.txt'  (this is now for the JDXConverter.py requirements)
   Open the MoleculesInfo.tsv 
   To add more molecueles, add rows with additional molecule names. 
   When adding new molecules, it is possible to leave the number of electrons, molecular weight, etc. as blank: JDXConverter will retrieve the missing information from online
   as long as the molecule exists at NIST Webbook.
   
   Add the name "methanol" at the bottom and save
   For acetone, remove the number of electrons and molecular weight, and save.
   
   To run the program, use an anaconda prompt because JDXConverter does not work from spyder.
   Use 'python JDXConverter.py' 
   During the run, you can just press enter twice since we are not using any special option.
   
   Open the directory "OutputFiles" and open the .tsv file for ConvertedSpectra with the highest number (should be 1).
   
   Now you know how to retrieve spectra and meta data from online.
   For use in any MSRESOLVE run, you can delete any undesired columns.
   NOTE: it is always best to have a reference pattern measured for *each* molecule. 
   However, using such retrieved reference patterns (which come from NIST) is reasonable for any molecules that does not have a directly measured pattern.
   MSRESOLVE does have features for extracting reference patterns from measured data.
    
03 How to do a typical analysis.
   Directory 01 has been copied to make directory 3.  However, all of the features normally used have been turned off.
   
   Open the UserInput.py file. Normally, a user must put their settings into this file before each analysis.
   It will be useful to use "ctrl+f" in the UserInput file.
   Change the data input file:
   UserChoices['inputFiles']['dataToAnalyzeFileName'] = 'Collected_GasPulse.csv'
   
   Now try to run runfile.py
   
   Look at the graphs. Note that there is some kind of early time behavior that is not the ethane pulse. Let's get rid of that.
   Turn on the "Time Range" feature, and put in a start and finish time that will remove that early region while keeping the main pulse.
   
   Now try to run runfile.py
   
   Let's add a baseline correction because that's a good practice (even though this data set does not seem to need it).

   Turn on the baseline feature, by changing:
   UserChoices['linearBaselineCorrectionSemiAutomatic']['on'] = ___
   UserChoices['linearBaselineCorrectionSemiAutomatic']['earlyBaselineTimes'] = ___
   UserChoices['linearBaselineCorrectionSemiAutomatic']['lateBaselineTimes'] =  ___
   Note: This is similar to section 3.2 of the quick start guide.
   
04 How to do a typical analysis continued:
  [DELETE 04 DIRECTORY IF IT EXISTS] 
   Copy Directory 3 and paste it. Then rename the new directory to start with "04"
   
   We are now going to use the sequential linear subtraction (SLS) feature of MSRESOLVE. This is a 'special' feature of MSRESOLVE, and in various cases, it produces more accureate results than the 'inverse' way (inverse is what most other programs would use). It will be useful to see some of the intermediate calculations that MSRESOLVE can export to file.
   
   First, change the user input to:
   UserChoices['dataAnalysisMethods']['solverChoice'] = 'sls'
   UserChoices['ExportAtEachStep']['on'] = 'yes'
   
   Now try to run runfile.py
   
   You will see an error message. Scroll up, and read what the messages printed out were. There is a warning that tells the user that MSRESOLVE was unable to solve the problem with the current settings. The warning notes that:
    "You may want to try using referenceMassFragmentFilterThreshold within the feature Reference Mass Fragmentation Threshold,"
    
   Why is SLS not working, and why is that recommendation made? To understand this, let's consider how SLS works:
   The way that SLS works is by approximating small fragments as 'zero' in order to solve what remains. Then, MSRESOLVE back corrects for those approximations before giving final outputs.
   
   For the example we're working on , let's open an exported reference pattern during the analysis, since that will be standardized to 100. This time, there is only one:
   Exported0ReferencePatternOriginalForCorrectionValues.csv
   
   What size of small fragmentation pattern peak do you think is necessary for SLS to work?
    
   Let's turn on the relevant feature and put that number in:
   
   UserChoices['applyReferenceMassFragmentsThresholds']
   UserChoices['applyReferenceMassFragmentsThresholds']['referenceMassFragmentFilterThreshold'] = _____
   UserChoices['applyReferenceMassFragmentsThresholds']['referenceSignificantFragmentThresholds'] = ___
   
   After running MSRESOLVE, compare to the inverse outputs in 03.
   We see that the chemicals are different in final concentrations, but that the conclusions aren't changed in this case.
   In general, the solutions can be affected by the solver, and 'sls' with careful choices should be used when possible.
     
05 How to do a typical analysis continued:   
    [DELETE 05 DIRECTORY IF IT EXISTS] 
   Make a copy of directory 4. Rename it to 05.
   ctrl + f  for sls, then change the solving for answer choice to 'inverse'.
   Now run MSRESOLVE.
   If comparing to 04, you will see that for **this** problem, and these settings, the output should be the same or similar.

   
06 In this example, we will use the concentration finder feature 
    [DELETE 06 DIRECTORY IF IT EXISTS] 
   Copy directory 05. (solver is inverse)
   UserChoices['concentrationFinder']['on'] = 'no'  #Change this to 'yes'.
   Consider the case that a signal of 1E-7 on m27 is known to be associated with 1 Torr of Ethane.  In this case, we can calibrate the concentrations accordingly.
   Turn the feature on:
   UserChoices['concentrationFinder']['on'] = 'yes'
   ['massNumberTSC_List'] = [27]
   ['moleculeSignalTSC_List'] = [1.0e-7]
   ['moleculeConcentrationTSC_List'] = 1
   Now run MSRESOLVE.
   
   The graphs that come out now include resolvedConcentrations. These are the same graphs as before, it's just the units that are different.
   The data is in: ResolvedConcentrations.tsv
       
   We can also use a different molecule, or separate factors for separate molecules.
   
   
07 Now let's make an example where want to add a molecule to our reference file. Let's say Ethanol
   First make a copy of an existing directory (such as 06).
   For the molecule we want to add, we can copy the data out of the JDXConverter output and place it in the reference file.
    --> first check the UserInput.py file to find out what the name of the reference file is.
    --> then open that file so that you can add the missing reference pattern.
    --> You can put in a full spectrum, but for this workshop to be quick, let's just paste in the intensities for masses 27, 29, and 31.
    --> after pasting, you'll see the molecule fragmentation patterns are not even on the same scale. That is fine. MSRESOLVEE will standardize them.
    --> Because we are using reference uncertainties from file, we need to update that file, too. The uncertainties file is *absolute** uncertainties, so make the uncertainties of the extent desired. 5% of the pattern? Or a fixed amount for all mass fragments? You can choose what you would like.
    We used 'inverse' for this example, because using sls for this example would be less easy. Using sls is usually best, but not always.
   
08 Tuning correction:
   Make a copy of directory 06 and rename it to 08.   
   For directory 08, we will implement the tuning corrector feature of MSRESOLVE.
   This feature allows for more accurate solved concentrations when using external reference patterns, as it corrects external fragment patterns for differences between different spectrometers.
   MSRESOVLE does this by comparing the spectra of some molecules that have been measured on both spectrometers and performing a tuning correction for the unmeasured molecule(s).
   
   The requirements for such analysis are as follows:
     UserChoices['inputFiles']['dataToAnalyzeFileName'] -- your measured data to convert to concentrations (as usual)
     UserChoices['inputFiles']['referencePatternsFileNamesList'] -- your reference file, with molecules measured on your spectrometer (as usual)
     UserChoices['tuningCorrection']['on'] -- We will need to change this to yes.
     UserChoices['tuningCorrection']['tuningCorrectPatternInternalVsExternal'] = 'External' #we are going to correct apattern from an external pattern, so we will not change this.
     UserChoices['tuningCorrection']['referenceFileStandardTuningAndForm'] -- #We will provide NIST patterns to MSRESOLVE as 'standard' tuning. If MSRESOLVE has a standard tuning file provided, that allows for even better tuning correction. 
     UserChoices['tuningCorrection']['referenceFileExistingTuningAndForm'] =[] #This is the pattern that will be pattern tuning corrected. We are not going to change this line because the standard tuning file will automatically be used for the existing pattern if none is provided.
     UserChoices['tuningCorrection']['referenceFileDesiredTuningAndForm'] =[] #This is the pattern that will be matched: leave this blank to use your same reference pattern for desired tuning.
   
   For this example we are going to include ethanol, again, but this time we are going to include the tuning corrected version.
   Let's take the following steps to get prepared.
     (a) Copy ConvertedSpectra1.tsv from the JDXConverter directory and place it in this one.
     (b) MSRESOLVE capitalization must be exact. Let's check to make sure the capitalization of the molecule names in our referencePatternsFile matches this external pattern, for the molecules that are in both files and will be used for tuning calibration:
             ConvertedSpectra1                ExtractedReferencePattern
             ethyne                           Acetylene
             ethene                           Ethylene
             ethane                           Ethane
        We see that the names used and capitaliation are different in ConvertedSpectra1.tsv. Let's change ExtractedReferencePattern to be the best practices naming and capitalization.
        
    Now we have everything ready and just need to change our input file according to the requirements noted above!     
     (1) populate our dataToAnalyzeFileName -- This is already done.
     (2) populate our referencePatternsFileNamesList -- This is already done
     (3) UserChoices['tuningCorrection']['on'] -- We will need to change this to yes.
     (4) UserChoices['tuningCorrection']['tuningCorrectPatternInternalVsExternal'] = 'External' #we are going to correct apattern from an external pattern, so we will not change this.
     (5) UserChoices['tuningCorrection']['referenceFileStandardTuningAndForm'] = ['ConvertedSpectra1.tsv', 'xyyy'] #Giving the external standard tuning pattenrn.
     (6) UserChoices['tuningCorrection']['referenceFileExistingTuningAndForm'] =[] #we won't change this because the standard tuning pattern will be used automatically if we leave this blank.
     (7) UserChoices['tuningCorrection']['referenceFileDesiredTuningAndForm'] =[] #we won't change this because the desired tuning is our internal refernce pattern, and that will be used by default if we leave this blank.
     (8) We need to do one more thing: ConvertedSpectra1 has many molecules, and we don't want to search for signs of all of them in our data -- that would be unsolvable. So after the tuning correction is done, we need to look for only the four molecules we want to solve for (not all of the molecules from the external patterns).
        UserChoices['specificMolecules']['on'] = 'yes'
        UserChoices['specificMolecules']['chosenMoleculesNames'] = ['ethane', 'ethyne', 'ethene', 'ethanol']

    Now run MSRESOLVE by the runfile. 
    The excel file in the SOLUTIONS directory (EthanolTuningCorrected08VersusUncorrected07.xlsx) shows that the change in the pattern from the tuning correction is small. Accordingly, if we check the scaledConcentrations graph for Ethanol in directory 07 and also in 08, we see that there is not much change. Whether a tuning correction makes a big difference in the concentrations or not is system specific. It should be performed when possible.
   
After the workshop, attendees may want to go through the "MSRESOLVE_QUICKSTART_AND_EXAMPLE_ANALYSIS" document, which is included with the workshop files.