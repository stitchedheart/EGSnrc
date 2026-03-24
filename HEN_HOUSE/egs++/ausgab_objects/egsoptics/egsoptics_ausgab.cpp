/*
##############################################################################
# 
# This is the cpp file for my PhD project, EGSoptics
# It is a implementation of optical and scintillation physics in EGS
# 
# Currently, the project is in its infancy with a lot of room for improvement
# 
#
#
##############################################################################
*/

#include <fstream>
#include <string>
#include <cstdlib>

#include "egsoptics_ausgab.h"
#include "egs_input.h"
#include "egs_functions.h"
#include "egs_ausgab_object.h"

EGS_Optics::EGS_Optics(const std::string &Name,EGS_ObjectFactory *f) 
    : EGS_AusgabObject(Name,f), scint_region(-1), min_edep(0.0), 
    birksCst(10.0), scintEff(10.0)
    
{

}

void EGS_Optics::setParameters(EGS_Input *input){

    if(!input) return;

    input->getInput("scint_region", scint_region);
    input->getInput("min_edep", min_edep);
    input->getInput("birksCst", birksCst);
    input->getInput("scintEff", scintEff);


    egsInformation("Ausgab initialized\n");
    egsInformation("scintillator region: %d\n", scint_region);
    egsInformation("Minimum energy: %4f\n", min_edep);
    egsInformation("Birk's constant: %4f\n", birksCst);
    egsInformation("Scintillator efficiency: %4f\n", scintEff);

}

EGS_Optics::~EGS_Optics(){

}

void EGS_Optics::setApplication(EGS_Application *App){
    EGS_AusgabObject::setApplication(App);

    if(!app) return;
}

extern "C" {

    EGS_AusgabObject *createAusgabObject(EGS_Input *input, 
                                         EGS_ObjectFactory *f) {
        
        EGS_Optics *result = new EGS_Optics("", f);    

        if (!result) return 0;

        result -> setParameters(input);

        return result;

}

}