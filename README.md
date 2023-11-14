# mucolstudies
Collection of scripts for performing Muon Collider studies. These scripts assume you're inside the singularity image used in the [Fermilab Muon Collider tutorial](https://mcdwiki.docs.cern.ch/tutorials/fermilab2022/computing_setup/) which does some path mapping for you on the Snowmass cluster. In particular, you should run:

```
apptainer build k4toroid.sif docker://madbaron/k4test-ubuntu:latest
apptainer run --no-home -B /tank/data/snowmass21/muonc:/data -B /home/$USER -B /work/$USER /home/$USER/mucol/k4toroid.sif
```

and then from inside the image, run:

```
source  /opt/spack/share/spack/setup-env.sh
spack env activate -V k4prod
source /opt/spack/opt/spack/linux-ubuntu22.04-x86_64/gcc-11.3.0/mucoll-stack-2023-07-30-ysejlaccel4azxh3bxzsdb7asuzxbfof/setup.sh
export MARLIN_DLL=$MARLIN_DLL:/opt/MyBIBUtils/lib/libMyBIBUtils.so
```

(FYI, putting these two in one shell script and sourcing it will not work; it won't execute the second command until you exit the singularity. You can also update to a newer release, which you can find [here](https://confluence.infn.it/display/muoncollider/Software), but I ran into trouble running with 1.7 with my current setup.)

For a very simple version of reading an slcio file, see `makeMuonPlots_simple.py`. 

For a slightly more advanced version that handles making multiple histograms more elegantly, see `makeMuonPlots.py`. It may be helpful to look at these side-by-side if you're trying to learn what's going on here.

Both of these scripts will put plots in a `plots/` directory - make one if you don't have it already or it will crash.

To work with key4hep files, look at `makeMuonPlots_k4hep.py`. This is just a skeleton file; it doesn't show you as much as even the `makeMuonPlots_simple.py` file when it comes to e.g. plotting, but it does give a simple recipe for how to open and loop over these files with the key4hep syntax. Be sure to read the comments at the top of this file for specific instructions on how to run the docker image for this.

For the SLCIO files, to get a list of the colletions you can access, run `anajob <filename>` on one of your files.
Use the `COLLECTION NAME` to access them. 
To understand what kinds of functions you can use on the particles in these collections, look at the `COLLECTION TYPE` and look up its functions [here](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/namespaceEVENT.html).

To understand in more depth how reconstruction works, you can look at the code [here](https://github.com/MuonColliderSoft/DDMarlinPandora/tree/master/src). To understand what options were passed to this code, you'll need to look at the steering files that were used to run it. Fede's are [here](https://github.com/madbaron/SteeringMacros/tree/master/Reco), but you'd have to know which one he used. For more generic ones, have a look at the ones used in the
[tutorial](https://github.com/MuonColliderSoft/MuC-Tutorial/tree/master/reconstruction).
