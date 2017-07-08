# Data-driven Modeling of Structured Populations: A Practical Guide to the Integral Projection Model

This is the companion git repository to "Data-driven Modeling of Structured Populations: A Practical Guide to the Integral Projection Model" (2016) by Ellner, Childs, and Rees. It contains R code to implement the many different examples in the book (plus a few other bits and pieces that never made it in). There are two ways to get hold of all of this. If you understand how git works you can just clone the repository and away you go. If you are not familiar with git, then download everything as a ZIP file and unzip it in a sensible location on your machine. 

The scripts have been written to use relative paths, assuming R's working directory has been set to the folder that the script is in. So for example if you have created a directory structure that includes "C:ipmbook/Rcode/c2", "C:/ipmbook/Rcode/c3", and so on, then before running scripts from the c3, inside R do setwd("C:/ipmbook/Rcode/c3"). On a Mac, it's even easier. If you have associated R with .R files, and start R by double-clicking on a script, R's working directory will be set to the folder that the script is in, and the scripts in that folder should "just work". 

A word of warning: we are not professional programmers. Like most scientists we are largely self-taught, which means: 1) we do not always write the most beautiful or efficient R code; 2) we each have a slightly different approach to getting things done in R. We have tried to be consistent, but occasionally our styles diverge. We are at least reasonably sure that everything works as advertised, but do let us know if you spot a problem. 

Happy IPMing!

Steve, Dylan and Mark

## Organisation

The organisation is fairly self-explanatory. Each chapter that contains one or more examples has its own subdirectory inside the *Rcode* directory (called *c2*, *c3*, *c4*, ...). Each file inside these chapter folders is either a self contained R script, or a script that depends on other scripts (included via `source()`). The *figures* folder contains a parrallel set of empty subdirectories (*c2*, *c3*, *c4*, etc). These are just here so that you can run the examples to produce the figures in the book "as is", i.e. the figures will get dumped into the appropriate directory.

## Corrections

There will inevitably be typos and omissions in a book of this size, and possibly some errors in the code examples. Please let one of us know if you think you've found one. We will use this repository to keep a living list of any problems that are discovered, as well as to note any changes we make to the R code examples post publication. There is also a file Errata.tex that lists all known errata in the book itself. [Direct link to the pdf](https://github.com/ipmbook/first-edition/raw/master/Errata.pdf)

