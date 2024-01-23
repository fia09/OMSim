/**
 * @file OMSimHitManager.hh
 * @brief Defines structures and classes related to optical module photon hit management.
 *
 * This file describes the `HitStats` structure which stores detailed photon hit parameters
 * and the `OMSimHitManager` class which centralizes and manages photon hit data across
 * different optical modules. The manager ensures unified access to photon hit information,
 * with functionalities to append, retrieve, and manipulate hit data.
 *
 * @ingroup common
 */

#ifndef OMSimHitManager_h
#define OMSimHitManager_h 1

#include "OMSimPMTResponse.hh"

#include <G4ThreeVector.hh>
#include <fstream>

/**
 * @struct HitStats
 * @brief A structure of vectors to store information about detected photons.
 * @ingroup common
 */
struct HitStats
{
    std::vector<G4long> event_id;                         ///< ID of the event
    std::vector<G4double> hit_time;                       ///< Time of detection.
    std::vector<G4double> photon_flight_time;             ///< Photon flight time.
    std::vector<G4double> photon_track_length;            ///< Length of the photon's path before hitting.
    std::vector<G4double> photon_energy;                  ///< Energy of the detected photon.
    std::vector<G4int> PMT_hit;                           ///< ID of the PMT that detected the photon.
    std::vector<G4ThreeVector> photon_direction;          ///< Momentum direction of the photon at the time of detection.
    std::vector<G4ThreeVector> photon_local_position;     ///< Local position of the detected photon within the PMT.
    std::vector<G4ThreeVector> photon_global_position;    ///< Global position of the detected photon.
    std::vector<G4double> event_distance;                 ///< Distance between generation and detection of photon.
    std::vector<OMSimPMTResponse::PMTPulse> PMT_response; ///< PMT's response to the detected photon, encapsulated as a `PMTPulse`.
};

/**
 * @class OMSimHitManager
 * @brief Manages detected photon information.
 *
 * This class is responsible for storing, managing, and providing access to
 * the hit information related to detected photons across multiple optical modules.
 * The manager can be accessed via a singleton pattern, ensuring a unified access
 * point for photon hit data.
 *
 * The hits are stored using 'OMSimHitManager::appendHitInfo' in 'OMSimTrackingAction::PostUserTrackingAction'.
 * The analysis manager of each study is in charge of writing the stored information into a file (see for example 'OMSimEffectiveAreaAnalyisis::writeScan' or 'OMSimDecaysAnalysis::writeHitInformation').
 * Once the data is written, do not forget calling 'OMSimHitManager::reset'!.
 *
 * @ingroup common
 */
class OMSimHitManager
{

public:
    /**
     * @brief Returns the singleton instance of the OMSimHitManager.
     * @return A reference to the singleton instance.
     */
    static OMSimHitManager &getInstance()
    {
        static OMSimHitManager instance;
        return instance;
    }

    void appendHitInfo(
        G4double globalTime,
        G4double localTime,
        G4double trackLength,
        G4double energy,
        G4int PMTHitNumber,
        G4ThreeVector momentumDirection,
        G4ThreeVector globalPos,
        G4ThreeVector localPos,
        G4double distance,
        OMSimPMTResponse::PMTPulse response,
        G4int moduleIndex = 0);
    void reset();
    std::vector<double> countHits(int moduleIndex = 0);
    void setNumberOfPMTs(int pNumberOfPMTs, int moduleIndex = 0);
    HitStats getHitsOfModule(int moduleIndex = 0);
    void sortHitStatsByTime(HitStats &lHits);
    std::vector<int> calculateMultiplicity(const G4double pTimeWindow, int moduleNumber = 0);
    std::map<G4int, HitStats> mModuleHits; ///< Map of a HitStats containing hit information for each simulated optical module

    G4int getNextDetectorIndex() {return ++mCurrentIndex; }
private:
    std::map<G4int, G4int> mNumPMTs; ///< Map of number of PMTs in the used optical modules
    G4int mCurrentIndex = -1;
    OMSimHitManager() = default;
    ~OMSimHitManager() = default;
    OMSimHitManager(const OMSimHitManager &) = delete;
    OMSimHitManager &operator=(const OMSimHitManager &) = delete;
};

#endif
