{
    printf("Running rootlogon.C\n");
    //gROOT->LoadMacro("/scratch/trholmes/utils/tdrstyle.C"); 
    //setTDRStyle();
    gROOT->LoadMacro("/scratch/trholmes/utils/AtlasStyle.C");
    SetAtlasStyle();
}
