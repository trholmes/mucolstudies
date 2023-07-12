# mucolstudies
Collection of scripts for performing mucol studies. 

For a very simple version of reading an slcio file, see `makeMuonPlots_simple.py`. 

For a slightly more advanced version that handles making multiple histograms more elegantly, see `makeMuonPlots.py`. It may be helpful to look at these side-by-side if you're trying to learn what's going on here.

Don't look at `makeMuonPlots_k4hep.py` -- it doesn't work yet. I currently can only use RDataFrame to read these files; anything else results in a segfault. Looking into it!

For the SLCIO files, to get a list of the colletions you can access, run `anajob <filename>` on one of your files.
Use the `COLLECTION NAME` to access them. 
To understand what kinds of functions you can use on the particles in these collections, look at the `COLLECTION TYPE` and look up its functions [here](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/namespaceEVENT.html).
