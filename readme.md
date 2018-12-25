# HBTEP

This is my personal code for data analysis on [HBTEP](http://sites.apam.columbia.edu/HBT-EP/) at [Columbia University](http://www.columbia.edu/).  

For now, the intent of this project is for my personal use, but longterm, I wish to be able to pass on a functional code base to future students working on HBTEP. 

You can contact me, John Brooks, at <jwb2159@columbia.edu> with any questions.  

## Getting started with git

I recommend the following youtube video:
[Git Tutorial for Beginners: Command-Line Fundamentals](https://www.youtube.com/watch?v=HVsySz-h9r4)

Also, this tutorial provides a much more thorough introduction to version control
[Atlassian version control tutorial](https://www.atlassian.com/git/tutorials/what-is-version-control)

## Install library from command line on Ubuntu

1. Make sure git is installed on your computer.  You can check by running `$ git --version`
2. Clone the repository.  Navigate to a directory in which you want to place the library.  Then type `$ git clone https://github.com/jwbrooks0/hbtepLib`  . An hbtepLib folder should now exist there.
3. Add this new directory to PYTHON_PATH.  You can check your PYTHON_PATH env. variable with the command `$ echo $PYTHON_PATH`.  If the directory in which you placed hbtepLib is not there, you'll need to add it.  Open ~/.profile and add the address to the hbtep folder.  Then run `$ source ~/.profile` to enable it for the current session.  Check the results by again typing `$ echo $PYTHON_PATH` .
4. Make sure you can import the new library.  Enter the python environment by typing `$ python`.  Then type `>>> import hbtepLib as hbt`.  You should get an error about a missing _hbtPreferences.py file.  This means everything is doing ok.
5. Make a copy _hbtPreferences.template and rename it as _hbtPreferences.py.  Then fill out the variables with your details.  
6. If you did step 5 correctly, re-typing `>>> import hbtepLib as hbt` in the python environment should work without error.   

You now have access to the library.  Make sure you watch the [Git Tutorial for Beginners: Command-Line Fundamentals](https://www.youtube.com/watch?v=HVsySz-h9r4) for instructions on how to work with the repository.  

## Installation on Windows

1.  Install python and whatever IDE you like.  I recommend [https://www.spyder-ide.org/](spyder).
2.  Download the latest stable release of [http://www.mdsplus.org/index.php/Latest_Windows_Distributions](MDSplus) for Windows
3.  At your python terminal, type `import MDSplus as mds`.  If successful, your MDSplus library is correctly installed.  
4.  Install a git client.  The internet seems to recommend [https://www.syntevo.com/smartgit/](SmartGit). 
5.  Clone this repository: https://github.com/jwbrooks0/hbtepLib
6.  Add the directory containing hbtepLib to PYTHONPATH.  Spyder has a PYTHONPATH manager under tools (and restart spyder).  
7.  At the terminal, type `import hbtepLib as hbt`.  It should give an error with something like `_hbtPreferences.py file not found` if correctly installed.  
8.  Make a copy _hbtPreferences.template and rename it as _hbtPreferences.py.  Then fill out the variables with your details.  If you don't know what to put in these variables, please ask someone at the lab.  
9.  At the terminal, again type  `import hbtepLib as hbt`.  There should be no error.  Then type, `hbt.get.tpData(95996,plot=True)` and you should get a plot.
10.  Play around with autocomplete to explore the various commands within this library.  I've done a decent job commenting the various functions, so please play around with tab-auto-complete and the object inspector (renamed to 'help' recently) in Spyder.  Note that if you run the library from command line, it will work just fine, but you won't have access to tab-auto-complete and the object inspector.

## Notes

You need to create a preferences.py file in the hbtep directory and populate it correctly if you want to use this code locally on your computer.  Check out the instructions above.  

## Tips for Spyder

I like Spyder for my pthon IDE for two reasons.  1) It has tab-auto-complete of the library and sublibraries which make finding the function you want very easy.  2) It has the object inspector (now called Help) which will instantly show any documentation (instructions) associated with any command you are writing.  In my code base, I've documented a lot (but not all) of the functions in order to make them easier to use.  To make use of this second feature, I recommend going to Preferences -> Help and selecting "Automatic connections" for both "editor" and IPython-Console"

