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

    // Transmission is simply 1 - R (energy conservation assumption)
    T = 1.0 - R;
}

class APP_EXPORT egsoptics_Application : public EGS_AdvancedApplication {
    public:
        egsoptics_Application(int argc, char **argv) :
            EGS_AdvancedApplication(argc,argv), save_case(-1) {srand(time(NULL));}
        int ausgab(int iarg);
        int initScoring();
        int save_case;

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

    // Checking to see if the particle is in a scintillator region, if not leave loop prematurely 
    if (incidentRegion != 2){
        return 0;
    } else {
        double Edep = getEdep();
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
    int incidentRegion = the_stack->ir[the_stack->np - 1] - 2;

    // Scintillation yield: photons produced per MeV of deposited energy
    double scintYield = 10000.0; //eventually could be read from input, for now using average value of photons per interaction
    //int n_photons = (int)(Edep*scintYield);
    int n_photons = 1;

    // If no photons are produced, exit early
    if (n_photons == 0) {
        //egsInformation("No photons generated for Edep = %.4f MeV\n", Edep);
        return;   // exit this function early
    }

    egsInformation("\n--- Starting stochastic optical photon transport test ---\n");

    //----------------------------------------------------------------------
    // Step 1: Generate initial photons
    //----------------------------------------------------------------------

    std::vector<OpticalPhoton> photons;
    photons.reserve(n_photons);

    for (int i = 0; i < n_photons; ++i){
        OpticalPhoton p;

        // Initial position — at the site of energy deposition
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

    double mfp_scatter = 0.05; // scattering mean free path (cm)
    int mfp_absorb = 50.0;   // absorption mean free path (cm)

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
        while (ph.alive && steps < 100) {
            steps++;

            // Sample free path to next scattering event
            double step = -mfp_scatter * log(randUniform());

            // Prepare EGS_Vector for position and direction
            EGS_Vector pos = { ph.x, ph.y, ph.z };
            EGS_Vector dir = { ph.u, ph.v, ph.w };

            // howfar() returns distance to next boundary and region beyond it
            EGS_Float t_boundary = 1e30;
            int newRegion = ph.region;
            int status = howfar(ph.region, pos, dir, t_boundary, &newRegion);

            egsInformation(" Step %d: pos=(%.3f,%.3f,%.3f), dir=(%.3f,%.3f,%.3f), intended_step=%.3f, t_boundary=%.3f, region=%d\n",
               steps, ph.x, ph.y, ph.z, ph.u, ph.v, ph.w, step, t_boundary, ph.region);

            if (status < 0) {
                // geometry error or particle exits geometry
                ph.alive = false;
                break;
            }

            //Limit the physical step to the geometric boundary distance.
            if (step > (double)t_boundary) {
                step = (double)t_boundary;
            }

            // Move photon
            ph.x += ph.u * step;
            ph.y += ph.v * step;
            ph.z += ph.w * step;

            
                //Absorption check with simple exponential attenuation
            if (randUniform() > exp(-step / mfp_absorb)) {
                ph.alive = false;
                break;
            }

            //If intended step reached the boundary exactly, apply Fresnel reflection / transmission.

            if ((double)t_boundary <= step) {
                // Refractive indices for current/next region (placeholder)
                double n1 = 1.34;  // incident medium (you may fetch from region)
                double n2 = 1.0;   // next medium (you may fetch from newRegion)
                // For now assume normal incidence cosθ = |w|
                double cos_theta_i = fabs(ph.w);

                double R = 0.0, T = 0.0;
                handleOpticalBoundary(n1, n2, cos_theta_i, R, T);

                if (randUniform() < R) {
                    // reflection
                    ph.w = - ph.w;
                    // step back slightly to avoid immediate re‑crossing
                    ph.x -= ph.u * 1e-6;
                    ph.y -= ph.v * 1e-6;
                    ph.z -= ph.w * 1e-6;
                } else {
                    // transmission (escape)
                    ph.region = newRegion;
                    ph.alive = false;
                    escaped = true;
                    break;
                }
            }

            //Isotropic scattering
            double theta = acos(1 - 2 * randUniform());
            double phi   = 2 * M_PI * randUniform();
            ph.u = sin(theta) * cos(phi);
            ph.v = sin(theta) * sin(phi);
            ph.w = cos(theta);
        }

    egsInformation("Photon %d: wavelength=%.1f nm, steps=%d, final=(%.2f, %.2f, %.2f), %s\n",
                   photon_id, ph.wavelength, steps, ph.x, ph.y, ph.z,
                   escaped ? "escaped" : "absorbed");

    if (escaped) total_escaped++;
    else         total_absorbed++;
}

    //----------------------------------------------------------------------
    // Step 4: Simulation summary
    //----------------------------------------------------------------------
    egsInformation("Summary: %d photons, %d escaped, %d absorbed\n", photons.size(), total_escaped, total_absorbed);
    egsInformation("--- Optical photon test finished ---\n");
}

APP_MAIN (egsoptics_Application);
