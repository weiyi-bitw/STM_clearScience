STM_clearScience
================

Source code for creating figures in STM paper

###Re-create predictions in the DREAM7 BCC challenge

The script for training the model using METABRIC and making predictions on OSLOVAL datasets is at:

https://raw.github.com/weiyi-bitw/BCCModels/master/runModel.r

But since the BCC package and predictiveModeling package have not been updated, if you are using R 3.0 and up, you will need to install them yourself from the source:


```
untar(download.packages(
"predictiveModeling", 
destdir="~", 
repos="http://depot.sagebase.org/CRAN/prod/2.15", 
type="source")[,2])

untar(download.packages(
"BCC", 
destdir="~", 
repos="http://depot.sagebase.org/CRAN/prod/2.15", 
type="source")[,2])
```

And do:

```
install.packages( c("~/predictiveModeling","~/BCC"), 
                  repos=NULL)
```

If you have any question regarding the model and the DreamBox7 package, please contact me.
If you have any question regarding installing BCC and predictiveModeling, please contact the Synapse support team (http://support.sagebase.org/sagebase) .
