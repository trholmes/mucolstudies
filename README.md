# mucolstudies
Collection of scripts for performing Muon Collider studies. These scripts assume you're inside the singularity image used in the [Fermilab Muon Collider tutorial](https://mcdwiki.docs.cern.ch/tutorials/fermilab2022/computing_setup/) which does some path mapping for you on the Snowmass cluster. In particular, you should run:

`singularity run -B /collab/project/snowmass21/data/muonc:/data -B /work/$USER /cvmfs/unpacked.cern.ch/registry.hub.docker.com/infnpd/mucoll-ilc-framework\:1.6-centos8`

and then from inside the image, run:

`source /opt/ilcsoft/muonc/init_ilcsoft.sh`

(FYI, putting these two in one shell script and sourcing it will not work; it won't execute the second command until you exit the singularity. You can also update to a newer release, which you can find [here](https://confluence.infn.it/display/muoncollider/Software).)

For a very simple version of reading an slcio file, see `makeMuonPlots_simple.py`. 

For a slightly more advanced version that handles making multiple histograms more elegantly, see `makeMuonPlots.py`. It may be helpful to look at these side-by-side if you're trying to learn what's going on here.

Don't look at `makeMuonPlots_k4hep.py` -- it doesn't work yet. I currently can only use RDataFrame to read these files; anything else results in a segfault. Looking into it!

For the SLCIO files, to get a list of the colletions you can access, run `anajob <filename>` on one of your files.
Use the `COLLECTION NAME` to access them. 
To understand what kinds of functions you can use on the particles in these collections, look at the `COLLECTION TYPE` and look up its functions [here](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/namespaceEVENT.html).

To understand in more depth how reconstruction works, you can look at the code [here](https://github.com/MuonColliderSoft/DDMarlinPandora/tree/master/src). To understand what options were passed to this code, you'll need to look at the steering files that were used to run it. Fede's are [here](https://github.com/madbaron/SteeringMacros/tree/master/Reco), but you'd have to know which one he used. For more generic ones, have a look at the ones used in the
[tutorial](https://github.com/MuonColliderSoft/MuC-Tutorial/tree/master/reconstruction).
