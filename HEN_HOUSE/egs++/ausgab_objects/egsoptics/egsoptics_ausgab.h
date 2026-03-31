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
#include "egs_advanced_application.h"
#include "egs_base_geometry.h"
#include "egs_application.h"
#include <string>
#include <random>

class EGS_Optics : public EGS_AusgabObject{

public:

    EGS_Optics(const std::string &Name="", EGS_ObjectFactory *f = 0);

    ~EGS_Optics();

    void setParameters(EGS_Input *input);

    int processEvent(EGS_Application::AusgabCall iarg){
        
        if (!app) return 0;

        int ir = app->top_p.ir;  //Index region for the particle
        EGS_Float edep = app->getEdep();  //Energy deposited in MeV
        double stepLength = app->getTVSTEP(); //in cm

        static bool firstPrint = true;
        if(firstPrint && iarg == 0) {
            std::cout << "EGS_Optics ausgab object successfully called!" << std::endl;
            firstPrint = false;  // prevent further prints
        }

        if (ir == scint_region && edep >= min_edep){
            //This is where code things need to go.
            //Update ir to be the scintillator region 
            //and edep to be the user inputted value
            //for minimum energy for optical photons
            //to be produced.

            double dEdx = edep/stepLength;

            double photonNum = (scintEff*dEdx*stepLength)/(1+(birksCst*dEdx));

            egsInformation("Scintillation hit in region %d, Edep = %.8f MeV\n", ir, edep);
            egsInformation("Mean number of photons created is: %.4f\n", photonNum);

            static std::mt19937 rng(std::random_device{}());

            std::poisson_distribution<int> poisson(photonNum);
            int actualPhotons = poisson(rng);

            egsInformation("Actual number of photons created is: %d\n", actualPhotons);
            
            if(edep>1e-6){
                double photonsMeV = actualPhotons/edep;
                egsInformation("Photons per MeV: %.2f\n", photonsMeV);
            }

            static double totalPhoton = 0.0;
            static double totalEdep = 0.0;

            totalPhoton += actualPhotons;
            totalEdep += edep;

            if (totalEdep > 0){
                double runningAvg = totalPhoton/totalEdep;
                egsInformation("Running average photon/MeV: %.0f\n", runningAvg);
            }

        }
        return 0;
    }

    bool needsCall(EGS_Application::AusgabCall iarg) const{
        return true;
        //return iarg == EGS_Application::AusgabCall::AFTER_STEP;
    }

    void setApplication(EGS_Application *App);

private:
        int scint_region;
        double min_edep;
        double birksCst;
        double scintEff;
};

#endif

