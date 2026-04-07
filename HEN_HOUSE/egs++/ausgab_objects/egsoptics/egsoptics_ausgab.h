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

    void readSpectrum(const std::string &filename);

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

            int actualPhotons = computePhotons(edep, stepLength);

            egsInformation("Scintillation hit in region %d, Edep = %.8f MeV\n", ir, edep);

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

            for(int i=0; i<actualPhotons; ++i){
                double photonEnergy = samplePhotonEnergy();
                egsInformation("Photon %d energy = %.10f MeV\n", i+1, photonEnergy);
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
        int computePhotons(double edep, double stepLength);
        double wavelength_nm_to_energy_MeV(double wavelength_nm);
        double samplePhotonEnergy();

        int scint_region;
        double min_edep;
        double birksCst;
        double scintEff;

        //For spectrum sampling
        std::vector<double> wavelengths; //in nm
        std::vector<double> cumulative;  //CDF 0 to 1
        std::mt19937 rng{std::random_device{}()};  // random number generator
        std::string spectrumFile;
        
        
};

#endif

