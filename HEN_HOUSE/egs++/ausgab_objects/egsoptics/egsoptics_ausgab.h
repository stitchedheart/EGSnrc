/*
##############################################################################
# 
# This is the header file for my PhD project, EGSoptics
# It is a implementation of optical and scintillation physics in EGS
# 
# Currently, the project is in its infancy with a lot of room for improvement
#
#
#
##############################################################################
*/

#ifndef EGS_OPTICS_
#define EGS_OPTICS_

#include "egs_ausgab_object.h"
#include "egs_base_geometry.h"
#include "egs_application.h"
#include <string>

class EGS_Optics : public EGS_AusgabObject{

public:

    EGS_Optics(const string &Name="", EGS_ObjectFactory *f = 0);

    ~EGS_Optics();

    int processEvent(EGS_Application::AusgabCall iarg){
        
        if (!app) return 0;

        int ir = app->top_p.ir;  //Index region for the particle
        EGS_Float edep = app->getEdep();  //Energy deposited

        static bool firstPrint = true;
        if(firstPrint && iarg == 0) {
            std::cout << "EGS_Optics ausgab object successfully called!" << std::endl;
            firstPrint = false;  // prevent further prints
        }

        if (ir >= 0 && edep > 0){
            //This is where code things need to go.
            //Update ir to be the scintillator region 
            //and edep to be the user inputted value
            //for minimum energy for optical photons
            //to be produced.

            egsInformation("Energy deposited in step is: %.4f MeV", edep);
        }
        return 0;
    }

    bool needsCall(EGS_Application::AusgabCall iarg) const{
        return true;
        //return iarg == EGS_Application::AusgabCall::AFTER_STEP;
    }

    void setApplication(EGS_Application *App);
};

#endif

