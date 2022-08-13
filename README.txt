These are the MSRESOLVE Workshop materials, put together for a workshop in August 2022.
Different directories exist and will also be made by the user for different examples of tasks that can be performed with MSRESOLVE.

PARTICIPANTS ARE EXPECTED TO HAVE ANACONDA INSTALLED AND TO HAVE TESTED THAT MSRESOLVE RUNS PRIOR TO THIS WORKSHOP.
IF THEY HAVE NOT ALREADY DONE SO, PARTICIPANTS SHOULD FOLLOW THE INSTRUCTIONS IN 000-ParticpantWorkshopPreparationInstructions.docx

If you have not already done so, install the required dependencies using the following command from spyder or an anaconda prompt: pip install -r requirements.txt
If pip does not work, first install pip and spyder from the anaconda interface or typing in "conda install pip" and "conda install spyder"
Python files can be run from an anaconda prompt by typing "python runfile.py" , or can be run by pressing the play button after opening in spyder.


00 An introductory slide.

   Before working on the exercises, we now need to get an anaconda window or spyder window open.
   

01  "Run" a premade typical analysis: 
    Navigate into that directory and run MSRESOLVE.py.  Various outputs will come out. Section 2.3 of the Quickstart explains what the various files are.
    ScaledConcentrations.tsv  and the graphs directory are the main outputs.
    Right click on ScaledConcentrations.tsv, click on properties, and click "change" to change the default program for opening the tsv to be excel. Note: you may need to find out where excel is on your computer.
    TSV files can also be dragged into excel when excel is already open.
    As noted in Section 3.1 of the Quickstart document, for this analysis there is no 'detectable' Ethene or Actylene (those are 'within error'). There is only Ethane. The abilty to have error bars are an important feature of MSRESOLVE

02 Creating a reference pattern file.
   Open the MoleculesInfo.tsv 
   To add more molecueles, add rows with additional molecule names. 
   When adding new molecules, it is possible to leave the number of electrons, molecular weight, etc. as blank: JDXConverter will retrieve the missing information from online
   as long as the molecule exists at NIST Webbook.
   
   Choose a molecule with a common name (such as acetone) and try removing the number of electrons and molecular weight.
   Then run JDXConverter.py. During the run, you can just press enter twice since we are not using any special options.
   If it does not run, use 'pip install -r requirements.txt'  (this is now for the JDXConverter.py requirements)
   
   Open the output files and open the tsv file.
   
   For use in any MSRESOLVE run, you can delete any undesired columns.
   NOTE: it is always best to have a reference pattern measured for *each* molecule. However, using individual NIST patterns is reasonable for any molecules that do not have a directly measured pattern.
   MSRESOLVE does have features for extracting reference patterns from data, which may occur in a later workshop exercise.
    
03 How to do a typical analysis.
   Directory 01 has been copied to make directory 3.  However, all of the features normally used have been turned off.
   
   Open the UserInput.py file. Normally, a user must put their settings into this file before each analysis.
   It will be useful to use "ctrl+f" in the UserInput file.
   Change the data input file:
   UserChoices['inputFiles']['dataToAnalyzeFileName'] = 'Collected_GasPulse.csv'
   
   Now try to run MSRESOLVE.py
   
   The uncertainties are not able to be calculated, but let's come back to that point a bit later.
   
   Look at the graphs. Note that there is some kind of early time behavior that is not the ethane pulse. Let's get rid of that.
   Turn on the "Time Range" feature, and put in a start and finish time that will remove that early region while keeping the main pulse.
   
   Now try to run MSRESOLVE.py
   
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
   
   Now try to run MSRESOLVE.py
   
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
   
   Follow the Section 3.2 of Quickstart -- plotting the data and changing user input settings.
   #TODO: change to a better data set. Probably use example 1 with some of the settings removed (like baseline correction etc.) so the person needs to do it again?
   
   referenceMassFragmentFilterThreshold from 0 to 1.01
   linearBaselineCorrectionSemiAutomatic to on and put in times
  
05 How to do a typical analysis continued:   
    [DELETE 05 DIRECTORY IF IT EXISTS] 
   Make a copy of directory 4. Rename it to 05.
   ctrl + f  for sls, then change the solving for answer choice to 'inverse'.
   Now run MSRESOLVE.
   If comparing to 04, you will see that for **this** problem, and these settings, the output is now the same.
     
   
06 In this example, we will use the concentration finder feature 
    [DELETE 06 DIRECTORY IF IT EXISTS] 
   Copy directory 05.
   Change back to sls.
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
    --> after pasting, you'll see the molecule fragmentation patterns are not even on the same scale. That is fine. MSRESOLVEE will standardize them.
    --> Because we are using reference uncertainties from file, we need to update that file, too. The uncertainties file is *absolute** uncertainties, so make the uncertainties of the extent desired. 5% of the pattern? Or a fixed amount for all mass fragments?
   We use 'inverse' for this example, because sls doesn't work well on this example.
   
08 Tuning correction:
   In general, different mass spectrometers behave differently. So let's correct the fact that the ethanol was collected in a different spectrometer.
   ConvertedSpectra1.tsv is from the JDXConverter.
   This time we intentionally start with **not** having ethanol in our regular file.
   But we need some common molecules between files, so we open ConvertedSpectra1.tsv and make sure to rename Ethylene, Ethane, and Acetylene to have the same capitalization.
   MSRESOLVE will *not* recognize the molecules unless the names are exactly the same.
   Next, in the UserInput, we will change tuningCorrection to yes and 
   UserChoices['tuningCorrection']['referenceFileStandardTuningAndForm'] = ['ConvertedSpectra1.tsv', 'xyyy']
   MSRESOLVE will then also use this for the existing tuning.

C:\Users\fvs\Documents\GitHub\MSRESOLVESG\UnitTests\TuningCorrector shows test_4.py:
MSRESOLVE.G.referencePatternsFileNamesList = ['AcetaldehydeMeasured.tsv']
MSRESOLVE.G.tuningCorrectPatternInternalVsExternal = 'External'
MSRESOLVE.G.tuningCorrection = 'yes'
MSRESOLVE.G.createMixedTuningPattern = True
MSRESOLVE.G.referenceFileExistingTuningAndForm = ['ReferenceLiterature.tsv','xyyy']
MSRESOLVE.G.referenceFileDesiredTuningAndForm =['ReferenceCollected.tsv','xyyy']

creates a mixed reference pattern where butanal comes from the    ReferenceLiterature.tsv
   --> in 08-Typical-Analysis, tried changing the field name to "Source:" from whatever it was.
   
   Currently, now see: ExportedReferencePatternExternalTuningCorrected.tsv 
   and ExportedReferencePatternStandardForCorrectionValuesMixedStandardTuningRearranged.tsv <-- this is strange, because it only has the original 3 molecules and Source:	NIST Webbook	NIST Webbook	NIST Webbook

   
After the workshop, attendees may want to go through the "MSRESOLVE_QUICKSTART_AND_EXAMPLE_ANALYSIS" document, which is included with the workshop files.