/**
 * @file OMSimIBD.hh
 * @brief This file contains the class in charge of generating positrons from IBD of supernova antineutrino interactions
 *
 * @author Cristian Jesus Lozano Mariscal
 * @version 10.1
 * @ingroup SN
 */


#ifndef OMSimIBD
#define OMSimIBD_1
 
#include "G4VUserPrimaryGeneratorAction.hh"

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "OMSimSNTools.hh"

class G4ParticleGun;
class G4Event;

/**
 * @class OMSimIBD
 * @brief Class in charge of generating positrons from IBD of supernova antineutrino interactions
 *
 * Supported by OMSimSNTools, this class generates the corresponding positrons after the inverse beta decay
 * and weights the interactions corresponding to the cross section and the generated volume
 *
 * @ingroup sngroup
 */
class OMSimIBD : public G4VUserPrimaryGeneratorAction
{
public:

	OMSimIBD(G4ParticleGun*);    
	~OMSimIBD(){};

	void GeneratePrimaries(G4Event* anEvent);

    G4int                  nPoints0;     //tabulated function
    G4int                  nPoints_lum;     //tabulated function
    std::vector<G4double>  x_lum;
    std::vector<G4double>  f_lum;           //f(x)
    std::vector<G4double>  a_lum;           //slopes
    std::vector<G4double>  Fc_lum;          //cumulative of f
    G4int                  nPoints1;     //tabulated function
    std::vector<G4double>  fixFe_X;
    std::vector<G4double>  fixFe_Y;           //f(x)
    std::vector<G4double>  fixFe_a;           //slopes
    std::vector<G4double>  fixFe_Fc;          //cumulative of f
    G4int fixE_nPoints;
    G4int                  angdist_nPoints;    
    std::vector<G4double>  angdist_x;
    std::vector<G4double>  angdist_y;           
    
    G4double ThresholdEnergy;
    G4double Delta;
    G4double y2;
    G4double NTargets;
    G4double me;
    G4double mp;
    G4double mn;
    G4double Gf;
    G4double consg;
    
    //here model data will be stored
	std::vector<G4double> mNuBar_time;
	std::vector<G4double> mNuBar_luminosity;
	std::vector<G4double> mNuBar_meanenergy;
	std::vector<G4double> mNuBar_meanenergysquare;

private:     
    G4ParticleGun* ParticleGun;

	void AngularDistribution(G4double Enubar);
    G4double PositronEnergy(G4double Enubar, G4double costheta);
    G4double TotalCrossSection(G4double energy);

    OMSimSNTools mSNToolBox;
    G4double mFixedenergy;
    G4double mFixedenergy2;
    G4double mAlpha;
};

#endif
