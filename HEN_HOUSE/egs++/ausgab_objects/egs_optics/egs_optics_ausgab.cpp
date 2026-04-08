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

#include "egs_optics_ausgab.h"
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

    input->getInput("scintillator regions", scint_region);
    input->getInput("mininum energy deposited", min_edep);
    input->getInput("birks constant", birksCst);
    input->getInput("scintillator efficiency", scintEff);

    input->getInput("scintillator spectrum file", spectrumFile);

    egsInformation("Ausgab initialized\n");
    egsInformation("scintillator region: %d\n", scint_region);
    egsInformation("Minimum energy: %4f\n", min_edep);
    egsInformation("Birk's constant: %4f\n", birksCst);
    egsInformation("Scintillator efficiency: %4f\n", scintEff);
    egsInformation("Spectrum file: %s/n", spectrumFile.c_str());

    if(!spectrumFile.empty()) {
        readSpectrum(spectrumFile);
    } else {
        egsInformation("No spectrum file provided, default spectrum will be used.\n");
    }   
}

EGS_Optics::~EGS_Optics(){

}

void EGS_Optics::setApplication(EGS_Application *App){
    EGS_AusgabObject::setApplication(App);

    if(!app) return;
}

int EGS_Optics::computePhotons(double edep, double stepLength) {

    if (stepLength <= 0) return 0;

    double dEdx = edep / stepLength;

    double meanPhotons = (scintEff * dEdx * stepLength) / (1 + (birksCst * dEdx));

    egsInformation("Mean number of photons created is: %.4f\n", meanPhotons);

    std::poisson_distribution<int> poisson(meanPhotons);

    return poisson(rng);
}

void EGS_Optics::readSpectrum(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        egsInformation("Could not open spectrum file: %s\n", filename.c_str());
        return;
    }

    wavelengths.clear();
    cumulative.clear();

    std::vector<double> intensities;
    double wl, inten;

    // Read wavelength and intensity pairs from the file
    while (file >> wl >> inten) {
        wavelengths.push_back(wl);
        intensities.push_back(inten);
    }

    if (wavelengths.empty()) {
        egsInformation("Spectrum file %s is empty or invalid.\n", filename.c_str());
        return;
    }

    // Build the cumulative distribution function (CDF)
    cumulative.resize(intensities.size());
    double total = std::accumulate(intensities.begin(), intensities.end(), 0.0);
    if (total == 0.0) total = 1.0; // avoid division by zero

    cumulative[0] = intensities[0] / total;
    for (size_t i = 1; i < intensities.size(); ++i) {
        cumulative[i] = cumulative[i-1] + intensities[i] / total;
    }
    cumulative.back() = 1.0; // ensure last value is exactly 1

    egsInformation("Spectrum loaded with %zu points.\n", wavelengths.size());
}

double EGS_Optics::wavelength_nm_to_energy_MeV(double wavelength_nm) {
    const double h = 4.135667696e-15;  // eV·s
    const double c = 2.99792458e8;     // m/s
    double energy_eV = h * c / (wavelength_nm * 1e-9);  // nm -> m
    return energy_eV * 1e-6;  // eV -> MeV
}

double EGS_Optics::samplePhotonEnergy() {
    // Make sure of valid spectrum
    if (wavelengths.empty() || cumulative.empty()) {
        return 0.0;
    }

    // Pick a random number between 0 and 1
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double r = dist(rng);

    // Find the wavelength corresponding to this random probability
    auto it = std::lower_bound(cumulative.begin(), cumulative.end(), r);
    size_t index = std::distance(cumulative.begin(), it);

    double wl = wavelengths[index];  // wavelength in nm

    // Convert wavelength to photon energy in MeV
    return wavelength_nm_to_energy_MeV(wl);
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