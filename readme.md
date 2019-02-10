# hbtepLib

I developed this code to download and analyze data from [HBTEP](http://sites.apam.columbia.edu/HBT-EP/) at [Columbia University](http://www.columbia.edu/). 

The code was originially developed for personal use, but I working on making it available to rest of our research group.  

You can contact me, John Brooks, at `jwb2159 at columbia dot edu` with any questions.  

## 1. Getting started with git

There is a somewhat steep learning curve with git.  I find it easier to start with a GUI interface like [SmartGit](https://www.syntevo.com/smartgit/) and eventually working down to the command line.  

Here are a few tutorials to get your started:
1.  [SmartGit beginners guide](https://www.youtube.com/watch?v=gB8OmhRJ0D8)
2.  [Git Tutorial for Beginners: Command-Line Fundamentals](https://www.youtube.com/watch?v=HVsySz-h9r4)
3.  [Atlassian version control tutorial](https://www.atlassian.com/git/tutorials/what-is-version-control)

## 2. Installing the library
The library is already installed on spizter, and you can skip this step if you only plan to use it there.  If you wish to install the library on your personal computer, then this section is for you.  

#### 2A. Installation on Ubuntu from the command line

This assumes you have the mdsplus python library already installed.  If you do not, see the later section below.  

1. Make sure git is installed on your computer.  You can check by running `$ git --version`
2. Clone the repository:  Navigate to a directory in which you want to place the library.  Then type `$ git clone https://github.com/jwbrooks0/hbtepLib`  . An hbtepLib folder should now exist there.
3. Add this new directory to PYTHON_PATH.  You can check your PYTHON_PATH env. variable with the command `$ echo $PYTHON_PATH`.  If the directory in which you placed hbtepLib is not there, you'll need to add it.  Open ~/.profile and add the address to the hbtep folder.  Then run `$ source ~/.profile` to enable it for the current session.  Check the results by again typing `$ echo $PYTHON_PATH` .
4. Make sure you can import the new library.  Enter the python environment by typing `$ python`.  Then type `>>> import hbtepLib as hbt`.  You should get an error about a missing _hbtPreferences.py file.  This means everything is doing ok.
5. Make a copy _hbtPreferences.template and rename it as _hbtPreferences.py.  Then fill out the variables with your details.  
6. If you did step 5 correctly, re-typing `>>> import hbtepLib as hbt` in the python environment should work without error.   

You now have access to the library.  Refer the tutorials above for how to manage the repo.  

#### 2B. Installation on Windows

This assumes you do not have mdsplus already installed on your computer.  

1.  Install python and whatever IDE (front-end gui for python) you like.  I recommend [https://www.spyder-ide.org/](spyder).  Others recommend jupiter notebooks.  The terminal interface is also *ok*, but it lacks much of the useful functionality that spyder will provide.  
2.  Download the latest stable release of [http://www.mdsplus.org/index.php/Latest_Windows_Distributions](MDSplus) for Windows
3.  At your python terminal, type `import MDSplus as mds`.  If successful, your MDSplus library is correctly installed.  
4.  Install a git client.  The internet seems to recommend [https://www.syntevo.com/smartgit/](SmartGit). 
5.  Clone this repository: https://github.com/jwbrooks0/hbtepLib
6.  Add the directory containing hbtepLib to PYTHONPATH.  Spyder has a PYTHONPATH manager under tools (and restart spyder).  
7.  At the terminal, type `import hbtepLib as hbt`.  It should give an error with something like `_hbtPreferences.py file not found` if correctly installed.  


## 3.  Setting up the library

After you have installed the library on your personal computer, you need to configure the preferences file.

1.   Start by making a copy of _hbtPreferences.template and renaming it as _hbtPreferences.py.  Then fill out `_HBT_SERVER_ADDRESS = ""` and `_HBT_SERVER_NAME = ""` with the correct values.  If you don't know what to put in these variables, please ask someone at the lab.  
2.  Return to python command line and type  `import hbtepLib as hbt`.  There should be no error.  Then type, `hbt.get.tpData(95996,plot=True)` and you should get a plot of the triple probe data.

## 4. Quick tutorial on using the library

1.  At the python terminal, type `import hbtepLib as hbt`.  If the library is installed correctly, you should receive no errors.
2.  Next, play around with autocomplete (in spyder) to explore the various commands within this library.  Read the *Tips for Spyder* sections below.  
3.  hbtepLib is a container for several sublibraries.  At the python terminal, type `hbt.` and then press tab to explore the libraries.  `get` is the primary library associated with downloading and analyzing HBT-EP data.  `process` has a number of data analysis tools (curve fitting, filtering, etc).  
4.  Next, type `hbt.get.` and press tab for autocomplete.  Most every diagnostic on HBT-EP has a function here to download the data.  In addition, many analysis codes (mode analysis, edge q, etc) also have functions.  
5.  Now let's start downloading some of the data.  At the python terminal, type `data=hbt.get.ipData(96530)`.  Now, type `data.` and hit tab.  You'll see that all of the data associated with the plasma current is stored here, including a few functions (such as one to plot itself).  `data.ip` is the I_p data, and `data.time` is the time (in seconds) associated with it.  
6.  Let's try a different example.  Now, type `data=hbt.get.qStarData(96496,plot=True)`.   As before, `data` is an object that contains all of the data as well as some subfunctions.  This time, however, a plot should automatically pop up showing q_star.  
7.  Another example.  This time, type `hbt.get.nModeData(96530,tStart=0.003,tStop=0.005,plot=True)`.  This will download and plot the data associated with the n=1 mode and also narrow the time window between 3 and 5 ms.  Most functins should allow a tStart and tStop.  
8.  Keep using auto complete and the Help menu in spyder to explore the various functionality of each sub-library.  Feel free to open up the various code files to see how they work.  

## Notes

You need to create a preferences.py file in the hbtep directory and populate it correctly if you want to use this code locally on your computer.  Check out the Step 3 above.  

## Tips for Spyder

I like Spyder for my pthon IDE for two reasons.  

1. It has tab-auto-complete of the library and sublibraries which make finding the function you want very easy.  

2. It has the object inspector (now called Help) which will instantly show any documentation (instructions) associated with any command you are writing.  In my code base, I've documented a lot (but not all) of the functions in order to make them easier to use.  To make use of this second feature, I recommend going to Preferences -> Help and selecting "Automatic connections" for both "editor" and IPython-Console"

## Tips for installing MDSplus

Installing MDSplus is typically the hardest part of this entire process.  Start [here](http://www.mdsplus.org/index.php?title=Downloads&open=82992333933183500305&page=Software%2FDownloads) and carefully follow their instructions.

##### Alternate method for installing MDSplus

I've only tried this once, but it worked on the first go for me.  It was also MUCH easier than doing the above.  

  1. Download the [MDSplus sources](http://www.mdsplus.org/index.php?title=Downloads&open=82992333933183500305&page=Software%2FDownloads#MDSplus_Sources) stable file.  [Direct link](https://github.com/MDSplus/mdsplus/archive/stable.zip).  
  2.  Extract the compressed files and navigate to ```/mdsplus-stable/mdsobjects``` 
  3.  Rename the ```plasma``` directory to ```MDSplus```.  This is case sensitive so be careful.
  4.  Move this directory to whereever you keep your python libraries.  
  5.  Add the directory where the MDSplus folder sits to your PYTHONPATH environment.  Spyder has an option under ```Tools``` to do this.  For Ubuntu, modfiy your .bashrc file as described [here](https://stackoverflow.com/questions/3402168/permanently-add-a-directory-to-pythonpath).

