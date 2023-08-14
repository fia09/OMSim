/**
 * @file OMSim.cpp
 * @brief Implementation of the OMSim class.
 * @ingroup common
 */

#include "OMSim.hh"
#include <boost/program_options.hpp>

namespace po = boost::program_options;
/**
 * @brief Construct a new OMSim object.
 *
 * This constructor initializes the Geant4 run manager, visualization manager,
 * and navigator. It also sets up the command line arguments for the simulation.
 */
OMSim::OMSim() : mGeneralArgs("General options")
{
    mStartingTime = clock() / CLOCKS_PER_SEC;
    mRunManager = new G4RunManager;
    mVisManager = new G4VisExecutive;
    mNavigator = new G4Navigator();
 
    mGeneralArgs.add_options()
    ("help", "produce help message")
    ("output_file,o", po::value<std::string>()->default_value("output"), "filename for output")
    ("numevents,n", po::value<G4int>()->default_value(0), "number of photons emitted per angle")
    ("visual,v", po::bool_switch()->default_value(false), "shows visualization of module after run")
	("save_args", po::bool_switch()->default_value(true), "if true a json file with the args and seed is saved")
	("seed", po::value<long>(), "seed for random engine. If none is given a seed from CPU time is used")
	("environment", po::value<G4int>()->default_value(0), "medium in which the setup is emmersed [AIR = 0, ice = 1, spice = 2]")
    ("detail_pmt",  po::bool_switch(), "if given, simulate PMT with internal reflections and thin photocathode")
    ("glass", po::value<G4int>()->default_value(0), "DEPRECATED. Index to select glass type [VITROVEX = 0, Chiba = 1, Kopp = 2, myVitroVex = 3, myChiba = 4, WOMQuartz = 5, fusedSilica = 6]")
	("gel", po::value<G4int>()->default_value(1), "DEPRECATED. Index to select gel type [Wacker = 0, Chiba = 1, IceCube = 2, Wacker_company = 3]")
	("reflective_surface", po::value<G4int>()->default_value(0), "DEPRECATED. Index to select reflective surface type [Refl_V95Gel = 0, Refl_V98Gel = 1, Refl_Aluminium = 2, Refl_Total98 = 3]")
	("pmt_model", po::value<G4int>()->default_value(0), "DEPRECATED. R15458 (mDOM) = 0,  R7081 (DOM) = 1, 4inch (LOM) = 2, R5912_20_100 (D-Egg)= 3");
}



OMSim::~OMSim()
{   
    delete mNavigator;
    delete mVisManager;
    delete mRunManager;
    double lFinishtime = clock() / CLOCKS_PER_SEC;
	G4cout << "Computation time: " << lFinishtime - mStartingTime << " seconds." << G4endl;
}

/**
 * @brief UIEx session is started for visualisation.
 */
void OMSim::startVisualisation()
{
	OMSimUIinterface &lUIinterface = OMSimUIinterface::getInstance();
    char *argumv[] = {"all", NULL};
    G4UIExecutive *UIEx = new G4UIExecutive(1, argumv);
    lUIinterface.applyCommand("/control/execute ../aux/init_vis.mac");
    UIEx->SessionStart();
    delete UIEx;

}


/**
 * @brief Ensure that the output directory for the simulation results exists.
 *
 * If the output directory does not exist, this function will create it.
 *
 * @param pFilePath The path to the output directory.
 */
void OMSim::ensureOutputDirectoryExists(const std::string &pFilePath)
{
    std::filesystem::path full_path(pFilePath);
    std::filesystem::path dir = full_path.parent_path();
    if (!dir.empty() && !std::filesystem::exists(dir))
    {
        std::filesystem::create_directories(dir);
    }
}


/**
 * @brief Initialize the simulation.
 *
 * This function sets up the necessary Geant4 components.
 */
void OMSim::initialiseSimulation()
{
    OMSimCommandArgsTable &lArgs = OMSimCommandArgsTable::getInstance();
    ensureOutputDirectoryExists(lArgs.get<std::string>("output_file"));

    std::string lFileName = lArgs.get<std::string>("output_file") + "_args.json";
    if (lArgs.get<bool>("save_args"))
        lArgs.writeToJson(lFileName);

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine(lArgs.get<long>("seed"), 3));



    mDetector = new OMSimDetectorConstruction();

    mRunManager->SetUserInitialization(mDetector);

    mPhysics = new OMSimPhysicsList;
    mRunManager->SetUserInitialization(mPhysics);

    mVisManager->Initialize();

    mGenAction = new OMSimPrimaryGeneratorAction();
    mRunManager->SetUserAction(mGenAction);

    mRunAction = new OMSimRunAction();
    mRunManager->SetUserAction(mRunAction);

    mEventAction = new OMSimEventAction();
    mRunManager->SetUserAction(mEventAction);

    mTracking = new OMSimTrackingAction();
    mRunManager->SetUserAction(mTracking);

    mStepping = new OMSimSteppingAction();
    mRunManager->SetUserAction(mStepping);

    mRunManager->Initialize();

    OMSimUIinterface &lUIinterface = OMSimUIinterface::getInstance();
    lUIinterface.setUI(G4UImanager::GetUIpointer());


    mNavigator->SetWorldVolume(mDetector->mWorldPhysical);
    mNavigator->LocateGlobalPointAndSetup(G4ThreeVector(0., 0., 0.));

    mHistory = mNavigator->CreateTouchableHistory();

    lUIinterface.applyCommand("/control/execute ", lArgs.get<bool>("visual"));

}