# HBTEP

This is my personal code for data analysis on [HBTEP](http://sites.apam.columbia.edu/HBT-EP/) at [Columbia University](http://www.columbia.edu/).  

For now, the intent of this project is for my personal use, but longterm, I wish to be able to pass on a functional code base to future students working on HBTEP. 

You can contact me, John Brooks, at <jwb2159@columbia.edu> with any questions.  

## Getting started with git

I recommend the following youtube video:
[Git Tutorial for Beginners: Command-Line Fundamentals](https://www.youtube.com/watch?v=HVsySz-h9r4)

## Installation from command line

1. Make sure git is installed on your computer.  You can check by running `$ git --version`
2. Clone the repository.  Navigate to a directory in which you want to place the library.  Then type `$ git clone https://github.com/jwbrooks0/hbtepLib`  . An hbtepLib folder should now exist there.
3. Add this new directory to PYTHON_PATH.  You can check your PYTHON_PATH env. variable with the command `$ echo $PYTHON_PATH`.  If the directory in which you placed hbtepLib is not there, you'll need to add it.  Open ~/.profile and add the address to the hbtep folder.  Then run `$ source ~/.profile` to enable it for the current session.  Check the results by again typing `$ echo $PYTHON_PATH` .
4. Make sure you can import the new library.  Enter the python environment by typing `$ python`.  Then type `>>> import hbtepLib as hbt`.  You should get an error about a missing _hbtPreferences.py file.  This means everything is doing ok.
5. Make a copy _hbtPreferences.template and rename it as _hbtPreferences.py.  Then fill out the variables with your details.  
6. If you did step 5 correctly, re-typing `>>> import hbtepLib as hbt` in the python environment should work without error.   

You now have access to the library.  Make sure you watch the [Git Tutorial for Beginners: Command-Line Fundamentals](https://www.youtube.com/watch?v=HVsySz-h9r4) for instructions on how to work with the repository.  


## Notes

You need to create a preferences.py file in the hbtep directory if you want to use this code.  This file should contain the variable:

_HBT_SERVER_ADDRESS = 'str'

which points at our data server.  Because I don't wish for this address to be public knowledge, you will need to get this address from me.  This is also why preferences.py is an ignored file with the project.  

When you clone the library from github, DO NOT INCLUDE www in the address.  I don't know why, but I'm not able to push if the project was setup with www in the address.  
