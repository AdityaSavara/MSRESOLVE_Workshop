These are the MSRESOLVE Workshop materials, put together for a workshop in August 2022.
It is intended for people to use these directories alongside the "Quickstart" document (which is included).
Different directories have different examples of tasks that can be performed with MSRESOLVE.

Python files can be run from an anaconda prompt with typing "python MSRESOLVE.py" , or can be run by pressing the play button after opening in spyder.


01 Has a typical analysis: Navigate into that directory and run MSRESOLVE.py.  Various outputs will come out. Section 2.3 of the Quickstart explains what the various files are.
    ScaledConcentrations.csv  and the graphs directory are the main outputs.
    As noted in section, Section 3.1 , there is no 'detectable' Ethene or Ethyne (Actylene). There is only Ethane. The error bars are an important feature of MSRESOLVE
    #TODO: Collected_Reference.csv  should be changed to Collected_Pulse.csv
    
02 is a learning opportunity for how to do a typical analysis.  This is following section 3.2 of the quick start guide.
   A dataset for an example analysis has already been made.
   Open the UserInput file. Normally, a user must edit this file.
   Follow the Section 3.2 of Quickstart -- plotting the data and changing user input settings.
   #TODO: change to a better data set. Probably use example 1 with some of the settings removed (like baseline correction etc.) so the person needs to do it again?
   
03 Way to make a reference file.
   Open the MoleculesInfo.tsv . Note that there are molecule names. Additional molecule names can be added. 
   It is possible to leave the number of electrons and molecular weight as blank: JDXConverter will retrieve the missing information from online
   as long as the molecule exists at NIST Webbook.
   Choose a molecule with a common name (such as acetone) and try removing the number of electrons and molecular weight.
   Then run JDXConverter.py  -- press enter twice since we are not using any special options.
   Open the output files and open the tsv file.
   Delete any undesired columns.
   
04 In this example, we will use the concentration finder feature 
    UserChoices['concentrationFinder']['on'] = 'no'  #Change this to 'yes'.
    Consider the case that a signal of 1E-7 on m27 is known to be associated with 1 Torr of Ethane.  In this case, we can calibrate the concentration.
    
    
    We can also use a different molecule: like m27