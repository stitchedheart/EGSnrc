/*
###############################################################################
#
#  EGSnrc egs++ tutor4pp application
#  Copyright (C) 2024 National Research Council Canada
#
#  This file is part of EGSnrc.
#
#  EGSnrc is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Affero General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  EGSnrc is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
#  more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with EGSnrc. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
#  Authors:         Frederic Tessier, 2024
#
#  Contributors:
#
###############################################################################
#
#  This application prints detailed history, stack and particle information
#  for tutorial purposes. It was inspired by the "monitor a section in detail"
#  part of an EGSnrc course laboratory that used to rely on tutor4 with the
#  IWATCH option turned on.
#
###############################################################################
*/


#include "egs_advanced_application.h"
#include "egs_interface2.h"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "egs_vector.h"
#include <fstream>
#include <string>

int getRegionIdForLabel(const std::string& inputFile, const std::string& targetLabel) {
    std::ifstream file(inputFile);  // Open the input file
    if (!file) {
        egsInformation("Error: Unable to open input file: %s\n", inputFile.c_str());
        return -1;  // Return error code
    }

    std::string line;
    while (std::getline(file, line)) {
        // Look for lines that define labels (set label = scint 2 in this case)
        if (line.find("set label =") != std::string::npos) {
            // Extract the label and region ID
            size_t labelPos = line.find("set label =") + 12;  // Position after "set label ="
            size_t regionIdPos = line.find_last_of(" ") + 1;   // Position of region ID after label

            std::string label = line.substr(labelPos, regionIdPos - labelPos - 1);  // Extract the label name
            int regionId = std::stoi(line.substr(regionIdPos));  // Extract the region ID

            // If this line corresponds to the target label, return the region ID
            if (label == targetLabel) {
                egsInformation("Success; Found label '%s' with region ID %d in file: %s\n", targetLabel.c_str(), regionId, inputFile.c_str());
                return regionId;
            }
        }
    }

    egsInformation("Error: Label '%s' not found in file %s\n", targetLabel.c_str(), inputFile.c_str());
    return -1;  // Return error code if label is not found
}


//Generates random floating point between 0 and 1.0
double randUniform(){
    //rand() returns a value between 0 and RAND_MAX and dividing by RAND_MAX normalizes
    return (double)rand() / RAND_MAX;
}

// Calculates the refractive index of a material as a function of its wavelength
// The model used is a simple empirical dispersion relation (two-term Cauchy equation):
// Describes how n varies with the wavelength 
double refractiveIndex(double wavelength_nm){
    double A = 1.32;
    double B = 3.0e4;
    double n = A + B / (wavelength_nm * wavelength_nm);
    return n;
}

double getRI(int region) {
    switch(region) {
        case 0: return 1.0;   // Vacuum
        case 1: return 2.0;   // Ta (tantalum, opaque, but placeholder as 2)
        case 2: return 1.58;  // Scintillator
        case 3: return 1.47;  // Glass
        default: 
            egsWarning("Unknown region %d, returning default of 1.0 (vacuum)\n", region);
            return 1.0;  // Default refractive index (vacuum)
    }
}

// Computes the reflection (R) and transmission (T) coefficients
// The calculation uses the Fresnel equations for unpolarized light
void handleOpticalBoundary(double n1, double n2, double cos_theta_i, double &R, double &T){
   
    // Clamp input to avoid domain errors in sqrt()
    if (cos_theta_i < -1.0) cos_theta_i = -1.0;
    if (cos_theta_i > 1.0) cos_theta_i = 1.0;

    // Sine of the incident angle
    double sin_theta_i = sqrt(1.0 - cos_theta_i*cos_theta_i);
    
    // Apply Snell’s law
    double sin_theta_t = (n1/n2) * sin_theta_i;

    // Check for total internal reflection (no transmission possible)
    if (sin_theta_t >= 1.0){
        R = 1.0; //Full reflection
        T = 0.0; //No transmission
        return;
    }

    // Cosine of the transmitted angle
    double cos_theta_t = sqrt(1.0 - sin_theta_t*sin_theta_t);

    // Fresnel reflection coefficients for s- and p-polarizations
    double Rs = pow((n1*cos_theta_i - n2*cos_theta_t) / (n1*cos_theta_i + n2*cos_theta_t), 2);
    double Rp = pow((n1*cos_theta_t - n2*cos_theta_i) / (n1*cos_theta_t + n2*cos_theta_i), 2);
    
    // Average for unpolarized light
    R = 0.5 * (Rs+Rp);

    //R is absolute value, but since squared, ignore that
    //R = (n1-n2)/(n1+n2);
    //R *= R;
    R = pow((n1-n2)/(n1+n2),2); 

    // Transmission is simply 1 - R (energy conservation assumption)
    T = 1.0 - R;
    egsInformation("Reflection Coefficient: R = %.3f, Transmission Coefficient: T = %.3f\n", R, T);
}

class APP_EXPORT egsoptics_Application : public EGS_AdvancedApplication {
    public:
        egsoptics_Application(int argc, char **argv) :
            EGS_AdvancedApplication(argc,argv), save_case(-1), pmt_hits(0) {srand(time(NULL)); scintillatorRegion = getRegionIdForLabel("slabModNew.egsinp", "scint");}
        int ausgab(int iarg);
        int initScoring();
        int save_case;

        int pmt_hits = 0;

        int scintillatorRegion;

        struct OpticalPhoton {
            double x, y, z;    // position in cm
            double u, v, w;    // cosines direction
            double wavelength; // wavelength in nm
            double weight;     // tracking weight
            bool alive;        // is the photon still active
            int region;        // used to track the region of the new particle
        };

        void testOpticalPhotonTransport(double x0, double y0, double z0, double Edep);
};

// init scoring
int egsoptics_Application::initScoring () {
    for (int call=AfterTransport; call<UnknownCall; ++call) setAusgabCall((AusgabCall)call, true);
    return 0;
}

// ausgab
int egsoptics_Application::ausgab (int iarg) {

    // index of the current particle and region
    int     np = the_stack->np - 1;
    int     ir = the_stack->ir[np] - 2; 
    int     iq = the_stack->iq[np];
    double  E  = the_stack->E[np];
    double  x  = the_stack->x[np];
    double  y  = the_stack->y[np];
    double  z  = the_stack->z[np];
    double  u  = the_stack->u[np];
    double  v  = the_stack->v[np];
    double  w  = the_stack->w[np];
    double  wt = the_stack->wt[np];
    int     lt = the_stack->latch[np];
    int     npold = the_stack->npold - 1;


    if (current_case != save_case) {
        egsInformation("\n\nHISTORY #%d\n", current_case);
        egsInformation("============================================================================================================================================\n");
        egsInformation("%5s  %-40s %3s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n", "iarg", "event", "NP", "charge", "energy", "region", "x", "y", "z", "u", "v", "w", "weight", "latch");
        egsInformation("============================================================================================================================================\n");
        egsInformation("\n%5s  %-40s ", "", "Incident particle");
        egsInformation("%3d%8d%8.3f%8d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8x\n\n", np+1, iq, E, ir, x, y, z, u, v, w, wt, lt);
        save_case = current_case;
    }

    // Get the region of the incident particle
    int incidentRegion = the_stack->ir[the_stack->np - 1] - 2;

    //egsInformation("The scintillator region is: %d", scintillatorRegion);

    // Checking to see if the particle is in a scintillator region, if not leave loop prematurely 
    if (incidentRegion != scintillatorRegion){
        return 0;
    } else {
        double Edep = getEdep();
        double step = getTVSTEP();
        //Calls testing code I wrote if the energy deposited is > 0
        if (Edep > 0.0){
            //egsInformation(">>> Scintillation event: Edep = %.4f MeV at (%.2f, %.2f, %.2f)/n", Edep, x, y, z);
            testOpticalPhotonTransport(x, y, z, Edep);
        }
    }

    // Echo selected interactions to terminal
    int echo = true;
    string str;
    switch (iarg) {
        case    BeforeTransport:    echo = false;                                   break;
        case    EgsCut:             str = "Energy below Ecut or Pcut";              break;
        case    PegsCut:            str = "Energy below AE or AP";                  break;
        case    UserDiscard:        str = "User discard";                           break;
        case    ExtraEnergy:        str = "Extra Energy deposited";                 break;
        case    AfterTransport:     echo = false;                                   break;
        case    BeforeBrems:        str = "Bremsstrahlung about to occur";          break;
        case    AfterBrems:         echo = false;                                   break;
        case    BeforeMoller:       str = "Moller scattering about to occur";       break;
        case    AfterMoller:        echo = false;                                   break;
        case    BeforeBhabha:       str = "Bhabha scattering about to occur";       break;
        case    AfterBhabha:        echo = false;                                   break;
        case    BeforeAnnihFlight:  str = "Annihilation in fligth about to occur";  break;
        case    AfterAnnihFlight:   echo = false;                                   break;
        case    BeforeAnnihRest:    str = "Annihilation at rest about to occur";    break;
        case    AfterAnnihRest:     echo = false;                                   break;
        case    BeforePair:         str = "Pair production about to occur";         break;
        case    AfterPair:          echo = false;                                   break;
        case    BeforeCompton:      str = "Compton scattering about to occur";      break;
        case    AfterCompton:       echo = false;                                   break;
        case    BeforePhoto:        str = "Photoelectric effect about to occur";    break;
        case    AfterPhoto:         echo = false;                                   break;
        case    BeforeRayleigh:     str = "Rayleigh scattering about to occur";     break;
        case    AfterRayleigh:      echo = false;                                   break;
        case    FluorescentEvent:   str = "Fluorescence event just occured";        break;
        case    CosterKronigEvent:  str = "Coster-Kronig event just occured";       break;
        case    AugerEvent:         str = "Auger event just occured";               break;
        case    BeforePhotoNuc:     str = "Photonuclear event about to occur";      break;
        case    AfterPhotoNuc:      echo = false;                                   break;
        default:                    echo = false;
    }
    if (echo) {
        egsInformation("--------------------------------------------------------------------------------------------------------------------------------------------\n");
        egsInformation("%5d  %-40s ", iarg, str.c_str());
        egsInformation("%3d%8d%8.3f%8d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8x\n", np+1, iq, E, ir, x, y, z, u, v, w, wt, lt);
    }

    // print particle stack after selected interactions
    switch (iarg) {
        case    AfterBrems:
        case    AfterMoller:
        case    AfterBhabha:
        case    AfterAnnihFlight:
        case    AfterAnnihRest:
        case    AfterPair:
        case    AfterCompton:
        case    AfterPhoto:
        case    AfterRayleigh:
        case    FluorescentEvent:
        case    CosterKronigEvent:
        case    AugerEvent:
        case    AfterPhotoNuc:
                    egsInformation("%5d  %-40s \n", iarg, "STACK after interaction:");
                    for (int i=0; i<=np; i++) {
                        ir = the_stack->ir[i] - 2; 
                        iq = the_stack->iq[i];
                        E  = the_stack->E[i];
                        x  = the_stack->x[i];
                        y  = the_stack->y[i];
                        z  = the_stack->z[i];
                        u  = the_stack->u[i];
                        v  = the_stack->v[i];
                        w  = the_stack->w[i];
                        wt = the_stack->wt[i];
                        lt = the_stack->latch[i];
                        egsInformation("%5s  %13d", ".", i+1);
                        string str = "";
                        if (i == npold) str += " <- NPold";
                        if (i == np)    str += " <- NP";
                        egsInformation("%-28s", str.c_str());
                        egsInformation("%3s%8d%8.3f%8d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8x\n", "", iq, E, ir, x, y, z, u, v, w, wt, lt);
                    }
                    break;
        default:    break;
    }

    return 0;
}

//--------------------------------------------------------------------------
//   Tests stochastic transport of optical photons in a simple scintillator
//   geometry. Simulates photon creation, scattering,
//   absorption, and boundary reflection/transmission.
//--------------------------------------------------------------------------

void egsoptics_Application::testOpticalPhotonTransport(double x0, double y0, double z0, double Edep){
    //egsInformation("\n--- Starting stochastic optical photon transport test ---\n");

    // Get the region of the incident particle
    int incidentRegion = the_stack->ir[the_stack->np - 1]-2;

    // Scintillation yield: photons produced per MeV of deposited energy
    double scintYield = 10000.0; //eventually would be read from input, for now using average value of photons per interaction
    //int n_photons = (int)(Edep*scintYield);
    int n_photons = 1;

    // If no photons are produced, exit early
    if (n_photons == 0) {
        //egsInformation("No photons generated for Edep = %.4f MeV\n", Edep);
        return;   // exit this function early
    }

    egsInformation("\n--- Starting stochastic optical photon transport test ---\n");
    egsInformation("NDGeometry region ID at photon start = %d\n", the_stack->ir[the_stack->np - 1]-2);


    //----------------------------------------------------------------------
    // Step 1: Generate initial photons
    //----------------------------------------------------------------------
    
    bool hitPMT = false;
    std::vector<OpticalPhoton> photons;
    photons.reserve(n_photons);

    for (int i = 0; i < n_photons; ++i){
        OpticalPhoton p;
        pmt_hits = 0;

        // Initial position — at the site of energy deposition
        const double epsilon = 1e-4;
        p.x = x0, p.y = y0, p.z = z0;

        // Random isotropic direction (theta, phi) in 3D
        double theta = acos(1-2*randUniform()); // cos(theta) uniformly distributed in [-1,1]
        double phi = 2 * M_PI * randUniform();
        p.u = sin(theta) * cos(phi), p.v = sin(theta) * sin(phi), p.w = cos(theta);

        // Random emission wavelength between 415 and 440 nm (most common wavelengths for scintillator light), fixed to 425 nm for now
        p.wavelength = 425;
        p.weight = 1.0;
        p.alive=true;

        // Set the region of the photon to be the same as the incident particle's region
        p.region = incidentRegion;

        photons.push_back(p);
    }

    egsInformation("Generated %d optical photons for Edep = %.4f MeV in region %d at (%.2f, %.2f, %.2f)\n", n_photons, Edep, incidentRegion, x0, y0 ,z0);

    //----------------------------------------------------------------------
    // Step 2: Transport parameters
    //----------------------------------------------------------------------

    double mfp_scatter = 0.1; // scattering mean free path (cm)
    double mfp_absorb = 3.0;  // absorption mean free path (cm)

    int photon_id = 0;
    int total_escaped = 0;
    int total_absorbed = 0;

    //----------------------------------------------------------------------
    // Step 3: Photon transport loop
    //----------------------------------------------------------------------

    for (auto &ph : photons) {
        photon_id++;
        int steps = 0;
        bool escaped = false;

    //  Loop over photon transport steps until photon dies (absorption or escape),
    //  max number of steps reached (safety limit)
        int currentreg = ph.region;
        while (ph.alive && steps < 100) {
            steps++;

            // Sample free path to next scattering event
            double step = -mfp_scatter * log(randUniform());
            double rawStep = step;

            // Prepare EGS_Vector for position and direction
            EGS_Vector pos = { ph.x, ph.y, ph.z };
            EGS_Vector dir = { ph.u, ph.v, ph.w };

            // howfar() returns distance to next boundary and region beyond it
            EGS_Float t_boundary = 1e30;
            int med;
            int prevreg = currentreg;
            int newreg = howfar(currentreg, pos, dir, t_boundary, &med);
            currentreg = newreg;

            egsInformation("current region = %d, new region = %d, intended step = %.3f, t_boundary = %.3f\n", prevreg, currentreg, rawStep, t_boundary);

            if (currentreg < 0) {
                // geometry error or particle exits geometry
                ph.alive = false;
                break;
            }

            //Limit the physical step to the geometric boundary distance.
            bool hitBoundary = (rawStep > t_boundary);

            if (hitBoundary) {
                step = t_boundary;
                //currentreg = newreg;
            }
            else{step = rawStep;}

            // Move photon
            ph.x += ph.u * step;
            ph.y += ph.v * step;
            ph.z += ph.w * step;

            egsInformation("Photon intended step: %.3f\n", step);

            //Absorption check with simple exponential attenuation
            if (randUniform() > exp(-step / mfp_absorb)) {
                ph.alive = false;
                break;
            }

            //If intended step reached the boundary exactly, apply Fresnel reflection / transmission.

            if (hitBoundary) {

                // Refractive indices
                double n1 = getRI(prevreg);
                double n2 = getRI(newreg);
                double cos_i = fabs(ph.w);
                //double cos_i = fabs(ph.u * nx + ph.v * ny + ph.w * nz);

                double R = 0.0, T = 0.0;
                handleOpticalBoundary(n1, n2, cos_i, R, T);

                // ----- Reflection -----
                if (randUniform() < R) {
                    // Reverse direction (simple reflection)
                    ph.u = -ph.u;
                    ph.v = -ph.v;
                    ph.w = -ph.w;

                    ph.x += ph.u * 1e-6;
                    ph.y += ph.v * 1e-6;
                    ph.z += ph.w * 1e-6;

                    egsInformation("Photon reflected at boundary\n");

                    ph.region = prevreg;
                } 
                
                // ----- Transmission -----
                else {
        
                    ph.x += ph.u * 1e-6;
                    ph.y += ph.v * 1e-6;
                    ph.z += ph.w * 1e-6;

                    ph.region = newreg;

                    if (newreg == 3) {
                        pmt_hits++;
                        hitPMT=true;
                        ph.alive = false;
                        escaped = false;
                        egsInformation("Photon reached PMT\n");
                    }
                    else if (newreg == 0) {
                        ph.alive = false;
                        escaped = true;
                        egsInformation("Photon escaped geometry\n");
                    }
                }
            }

            // Only scatter isotropically if photon is still in scintillator
            //bool forceDirect = false;
            //if (ph.alive && ph.region == 2 && !forceDirect) {
                //double theta = acos(1 - 2 * randUniform());
                //double phi   = 2 * M_PI * randUniform();
                //ph.u = sin(theta) * cos(phi);
                //ph.v = sin(theta) * sin(phi);
                //ph.w = cos(theta);
            //}

        }   

        egsInformation("Photon %d: wavelength=%.1f nm, steps=%d, final=(%.2f, %.2f, %.2f)\n",
                   photon_id, ph.wavelength, steps, ph.x, ph.y, ph.z);
        egsInformation("Photon direction: u = %.3f, v = %.3f, w = %.3f\n", ph.u, ph.v, ph.w);

        if (escaped) {total_escaped++;}
        else if (ph.region == 3) {}
        else         {total_absorbed++;}
    }

    //----------------------------------------------------------------------
    // Step 4: Simulation summary
    //----------------------------------------------------------------------
    egsInformation("Summary: %d photons, %d escaped, %d absorbed, %d hit PMT", photons.size(), total_escaped, total_absorbed, pmt_hits);
    egsInformation("--- Optical photon test finished ---\n");
}

APP_MAIN (egsoptics_Application);
