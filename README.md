# KazanViewer
Visualization and data processing GUI for EPR and NMR
  
Boris Epel* and Alexey Silakov**
  
* University of Chicago
** Pennsylvania State University
Many people in the world are using Matlab™ software package by Mathworks™ for the data processing and analysis. And many of them take the advantage of its console interface. However, sometimes the point-and-click type of operation is more convenient, especially for a big amount of datasets and routine processing. Matlab has possibilities to create the Graphical User Interface. In Kazan Viewer we have combined the calculation power of Matlab with the elaborated GUI .

Here we present the Matlab™ GUI, which help us to process data that are coming out of pulse and CW EPR, ENDOR, ELDOR and NMR experiments. Kazan Viewer behaves as a stand-alone application though it can easily exchange data with the Matlab shell. It is a freeware and has an open code so feel free to use and modify it. We just ask to leave our credits and comments in the text and inform us about useful additions. The latest source code (Matlab scripts) is available here or below.

KAZAN Viewer has main window(viewer) and data-processing plugins. Viewer has file load/save possibilities, 1D/2D plots and set of assisting tools. The file loader understands many data formats related to EPR and NMR experiments. Most of of the standard plot options are available through the menu. The common for EPR and NMR data representations like ppm or g-factor scale are implemented. Viewer can also produce standard figures, which can be used f.e. for printing.

The main advantage of KAZAN Viewer is its flexibility. It is as flexible as Matlab is. User can load data from the shell or utilize easy-to-write custom GUIs - plugins. Plugins take data from Viewer and generate output sets of data, which are visualized then by Viewer. They can have the routines of any complexity, execute external programs etc. Many of plugins generate the Matlab script, which could be used outside the viewer. Some day the rules for plugin writing will be summarized in the help. Meanwhile one can use plugins that are written by us as example. Data comparison, FFT, simple math, EPR/ENDOR and HYSCORE simulations and some other plugins are included.
   
Description
   
More or less standard file browser with possibilities to memorize favorite directories (manually); 

Support for most common file standards including:
 

... Bruker XEPR / XWINNMR;
... ASCII(most common);
... TecMag (MacNMR);
... RMN;
... Weizmann .d00 (Jaap software);
... SpecMan .d01
... Bruker FTIR Opus format (partially)
 

Real/Imaginary/Magnitude plot for data visualization;

1D/stacked 2D/contour/density, real/imaginary/magnitude plots for data visualization;

Custom problem-dependent plugins (additional GUIs) that can handle specific processing task;

In the package we include several scripts for simulating various pulse EPR measurements such as ESEEM (and HYSCORE), EPR (for various approximations) and ENDOR.

File history;
 
...and many other things - just try it.
  
  
Download


Below under "Attachments" you will find some older versions of Kazan Viewer; use arrows on the right to download.

Install
1) Unpack the archive to a conveniet directory using an archive managers such as WinZip or WinRar (there are various programs in Linux for that as well).
 
2) In MatLab (you have it, don't you?) use "Set Path" menu or "addpath" command to add the directory to the search list.
 
2.5) If your system is not "win32", you will need to recompile "*.c" routines from the KV for proper work of some simulation programs. Assuming that "mex" had been configured correctly {check Mathworks.com for problems}, just type mex <filename>.c. 
3) Type "kazan" to start the program
4) Have fun

NOTE. Most of the functions of Kazan Viewer are independent of any side programs / scripts. However, some plugins are meant to be an interface for simulation programs such as from Easyspin. Moreover, some of our simulation routines utilize basic functions from Easyspin package (e.g. "sop", "planck", "lshape" etc).
NOTE 2. Some of the oldest KV simulation programs might not work anymore due to changes in definition/interface of some functions from Easyspin. I try to trace that stuff but maninly for those KV routines that I frequently use myself. Just let me know if you have any problems ...
