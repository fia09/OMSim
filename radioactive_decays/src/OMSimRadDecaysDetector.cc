#include "OMSimRadDecaysDetector.hh"
#include "OMSimMDOM.hh"
#include "OMSimPDOM.hh"
#include "OMSimLOM16.hh"
#include "OMSimLOM18.hh"
#include "OMSimDEGG.hh"
#include "OMSimCommandArgsTable.hh"
#include "OMSimHitManager.hh"
#include "G4SDManager.hh"
#include "OMSimSensitiveDetector.hh"
#include "G4Box.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4NistManager.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4AssemblyVolume.hh"


double hc_eVnm = 1239.84193;
/**
 * @brief Constructs the world volume (sphere).
 */
void OMSimRadDecaysDetector::constructWorld()
{
    mWorldSolid = new G4Box("World", 370*cm, (370)*cm, 370*cm);
    mWorldLogical = new G4LogicalVolume(mWorldSolid, mData->getMaterial("Ri_Air"), "World_log", 0, 0, 0);
    mWorldPhysical = new G4PVPlacement (0, G4ThreeVector(0.,0.,0.), mWorldLogical, "World_phys", 0, false, 0);
    mWorldLogical->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 0.0, 0.0, 0.0)));

    G4bool checkOverlaps = true;

    G4Box *lBoxPool = new G4Box("Pool", (346.8/2) * cm, (168.6/2+1)* cm, 144/2 * cm);  // Pool
    G4Box *lWaterSolid = new G4Box("Water", 345.8/2 * cm, 167.6/2 * cm, 143/2 * cm);  // Water
    G4LogicalVolume *lPoolLogical = new G4LogicalVolume(lBoxPool, mData->getMaterial("NoOptic_Absorber"), "Pool", 0, 0, 0);
    G4PVPlacement *lPoolPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0 * cm), lPoolLogical, "Pool_phys", mWorldLogical, false, 0, checkOverlaps);
    lPoolLogical->SetVisAttributes(G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.2)));

    mWaterLogical = new G4LogicalVolume(lWaterSolid, mData->getMaterial("Ri_Water"), "Water_log", 0, 0, 0);
    G4VPhysicalVolume *lWaterPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0* cm), mWaterLogical, "Water_phys", lPoolLogical, false, 0, checkOverlaps);
    G4VisAttributes *Water_vis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.2));
    mWaterLogical->SetVisAttributes(Water_vis);

    G4LogicalBorderSurface *watersurface = new G4LogicalBorderSurface("Water_skin", lWaterPhysical, lPoolPhysical, mData->getOpticalSurface("Refl_tank_plastic"));

    // balls
    G4NistManager *MatDatBase = G4NistManager::Instance();
    double ball_rad = 2.5; // radius of ball in cm
    double dip_depth = 0;  // dip depth of ball into water in cm (0 == half in water)
    G4VisAttributes *Alu_vis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.9, 1.0));
    double ambient_temperature = (20 + 273.15) * kelvin;
    double ambient_pressure = 1.013 * bar;
    G4Material *Mat_Absorber = new G4Material("Absorber Black Paint", 1.0 * g / cm3, 1, kStateSolid, ambient_temperature);
    Mat_Absorber->AddMaterial(MatDatBase->FindOrBuildMaterial("G4_AIR"), 100.0 * perCent);
    G4bool gBallGrid = true;

    // make a dense ball grid
    if (gBallGrid)
    {
        G4Ellipsoid *lBall_solid = new G4Ellipsoid("ball", ball_rad * cm, ball_rad * cm, ball_rad * cm, (-ball_rad - 5) * cm, dip_depth * cm);
        G4LogicalVolume *lBall_logical = new G4LogicalVolume(lBall_solid, Mat_Absorber, "Ball logical");

        G4OpticalSurface *ball_optical = new G4OpticalSurface("ball optical");
        G4LogicalSkinSurface *ballsurface = new G4LogicalSkinSurface("ball_skin", lBall_logical, ball_optical);
        G4double BallPhotonEnergy[2] = {hc_eVnm / 800 * eV, hc_eVnm / 248 * eV};
        G4double BallReflectivity[2] = {0.04, 0.04};
        G4MaterialPropertiesTable *ballOpticalProperties = new G4MaterialPropertiesTable();
        ballOpticalProperties->AddProperty("REFLECTIVITY", BallPhotonEnergy, BallReflectivity, 2);
        ball_optical->SetMaterialPropertiesTable(ballOpticalProperties);
        lBall_logical->SetVisAttributes(Alu_vis);

        // Define one layer as one assembly volume
        G4AssemblyVolume *assemblyDetector = new G4AssemblyVolume();

        // rotation not really needed for the pool but maybe for formalities oh well
        //  Rotation and translation of a plate(ball) inside the assembly
        G4RotationMatrix Ra;
        G4ThreeVector Ta;
        G4Transform3D Tr;
        // Rotation of the assembly inside the world
        G4RotationMatrix Rm;
        // make one ball to multiply later (layer could also consist of multiple objects)
        Ta.setX((-345.8/2 + ball_rad) * cm);
        Ta.setY((-167.6/2 + ball_rad) * cm);
        Ta.setZ(143/2 * cm - dip_depth * cm);
        Tr = G4Transform3D(Ra, Ta);
        assemblyDetector->AddPlacedVolume(lBall_logical, Tr);

        // number of balls next to each other (y)
        int yseries_len = std::floor((167.6)/(2*ball_rad));
        int xseries_len = std::floor((345.8-2*ball_rad)/(ball_rad*std::sqrt(3))+1);

         for (unsigned int i = 0; i < xseries_len; i++)
        {
            if (i % 2 == 0)
            {
                for (unsigned int j = 0; j < yseries_len; j++)
                {
                    G4ThreeVector Tm((i * 2 * ball_rad * std::sqrt(3) / 2) * cm, j * 2 * ball_rad * cm, 0 * cm);
                    Tr = G4Transform3D(Rm, Tm);
                    assemblyDetector->MakeImprint(mWaterLogical, Tr);
                }
            }
            else
            {
                for (unsigned int j = 0; j < yseries_len - 1; j++)
                {
                    G4ThreeVector Tm((i * 2 * ball_rad * std::sqrt(3) / 2) * cm, (j * 2 * ball_rad + ball_rad) * cm, 0 * cm);
                    Tr = G4Transform3D(Rm, Tm);
                    assemblyDetector->MakeImprint(mWaterLogical, Tr);
                }
            }
        }
    }
}



/**
 * @brief Constructs the selected detector from the command line argument.
 */
void OMSimRadDecaysDetector::constructDetector()
{
    OMSimHitManager &lHitManager = OMSimHitManager::getInstance();
    bool lPlaceHarness = OMSimCommandArgsTable::getInstance().get<bool>("place_harness");

    OMSimOpticalModule *lOpticalModule = nullptr;

    switch (OMSimCommandArgsTable::getInstance().get<G4int>("detector_type"))
    {

    case 0:
    {
        log_critical("No custom detector implemented!");
        break;
    }
    case 1:
    {
        log_info("Constructing single PMT");
        OMSimPMTConstruction *lPMTManager = new OMSimPMTConstruction(mData);
        lPMTManager->selectPMT("argPMT");
        lPMTManager->construction();
        lPMTManager->placeIt(G4ThreeVector(0, 0, 0), G4RotationMatrix(), mWorldLogical, "_0");
        lHitManager.setNumberOfPMTs(1, 0);
        lPMTManager->configureSensitiveVolume(this, "/PMT/0");
        break;
    }
    case 2:
    {
        lOpticalModule = new mDOM(mData, lPlaceHarness);
        break;
    }
    case 3:
    {

        lOpticalModule = new pDOM(mData, lPlaceHarness);
        break;
    }
    case 4:
    {

        lOpticalModule = new LOM16(mData, lPlaceHarness);
        break;
    }
    case 5:
    {

        lOpticalModule = new LOM18(mData, lPlaceHarness);
        break;
    }
    case 6:
    {
        lOpticalModule = new DEGG(mData, lPlaceHarness);
        break;
    }
    }

    if (lOpticalModule)
    {
        lOpticalModule->placeIt(G4ThreeVector(140*cm, 0*cm, 0*cm), G4RotationMatrix(), mWaterLogical, "");
        lOpticalModule->configureSensitiveVolume(this);
        mOpticalModule = lOpticalModule;
    }
}
