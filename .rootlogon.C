{
    printf("Running rootlogon.C\n");
    //gROOT->LoadMacro("/home/trholmes/utils/tdrstyle.C"); 
    //setTDRStyle();
    if (gROOT->LoadMacro("/home/trholmes/utils/AtlasStyle.C") == 0) {
      SetAtlasStyle();
    } else {
      printf("Canceling rootlogon.C\n");
    }
}
