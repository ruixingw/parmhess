# rxcclib


## What is it?

*rxcclib* is a code package that defines molecules, MM functions, and methods for file parsing. It is the fundamental library for our other programs such as Parmhess, Tsubasa, and Katachi.

## Dependencies

*rxcclib* requires following packages.
- [Python3](https://www.python.org/) ([Anaconda](https://www.continuum.io/downloads) is strongly recommended)
- [numpy](http://www.numpy.org/) (Included in Anaconda)

## Install

1. Clone this repository by ``` git clone https://github.com/ruixingw/rxcclib.git```.
2. 
    a. Move *rxcclib* directory to your Python install's **site-packages** directory. If you use Anaconda, it is like ```anaconda/lib/python3.5/site-packages```.

    b. Or, You can move *rxcclib* to your own directory, and set $PYTHONPATH environment. For example, if *rxcclib* is put at ```/home/user/mypythonpath/rxcclib```, then ```export PYTHONPATH=/home/user/mypythonpath```(bash) or ```setenv PYTHONPATH "/home/user/mypythonpath"```(C shell).
