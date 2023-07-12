# mucolstudies
Collection of scripts for performing mucol studies. 

For a very simple version of reading an slcio file, see `makeMuonPlots_simple.py`. 

For a slightly more advanced version that handles making multiple histograms more elegantly, see `makeMuonPlots.py`. It may be helpful to look at these side-by-side if you're trying to learn what's going on here.

Don't look at `makeMuonPlots_k4hep.py` -- it doesn't work yet. I currently can only use RDataFrame to read these files; anything else results in a segfault. Looking into it! 
