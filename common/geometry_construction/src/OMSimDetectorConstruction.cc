#include "OMSimDetectorConstruction.hh"
#include "OMSimPDOM.hh"
#include "OMSimLOM16.hh"
#include "OMSimLOM18.hh"
#include "OMSimDEGG.hh"
#include "OMSimMDOM.hh"
#include "OMSimCommandArgsTable.hh"
#include "OMSimHitManager.hh"
#include "G4SDManager.hh"

OMSimDetectorConstruction::OMSimDetectorConstruction()
    : mWorldSolid(0), mWorldLogical(0), mWorldPhysical(0)
{
}

OMSimDetectorConstruction::~OMSimDetectorConstruction()
{
    delete mData;
}

/**
 * @brief Constructs the selected detector from the command line argument and returns the physical world volume.
 * @return Pointer to the physical world volume
 */
G4VPhysicalVolume *OMSimDetectorConstruction::Construct()
{
    mData = new InputDataManager();
    mData->searchFolders();

    constructWorld();
    constructDetector();

    return mWorldPhysical;
}

void OMSimDetectorConstruction::setSensitiveDetector(G4LogicalVolume *pLogVol, G4VSensitiveDetector *pSD)
{
    auto lSDManager = G4SDManager::GetSDMpointer();
    lSDManager->AddNewDetector(pSD);
    SetSensitiveDetector(pLogVol, pSD);
}