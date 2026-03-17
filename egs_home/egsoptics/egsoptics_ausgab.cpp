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

EGS_Optics::EGS_Optics(const string &Name,EGS_ObjectFactory *f) 
    : EGS_AusgabObject(Name,f)
    
{
    //Initialize variables here later
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
        
        result -> setName("OpticsParams");

        return result;

}

}