#ifndef OMSimDecaysAnalysis_h
#define OMSimDecaysAnalysis_h 1

#include <G4ThreeVector.hh>
#include <fstream>


/**
 * @struct DecayStats
 * @brief A structure to store information about decays.
 * @ingroup radioactive
 */
struct DecayStats
{
    std::vector<G4long> event_id;                 ///< ID of the event
    std::vector<G4String> isotope_name;           ///< Isotope name and energy level
    std::vector<G4double> decay_time;             ///< Time of the decay after (possibly) randomising inside simulation time window
    std::vector<G4ThreeVector> decay_position;    ///< Global position of the decay.
};



class OMSimDecaysAnalysis
{
public:
    static OMSimDecaysAnalysis &getInstance()
    {
        static OMSimDecaysAnalysis instance;
        return instance;
    }
    void appendDecay(G4String pParticleName, G4double pDecayTime, G4ThreeVector pDecayPosition);
    void setOutputFileName(G4String pName);
    void writeMultiplicity();
    void writeDecayInformation();
    void writeHitInformation();

    void reset();
    

private:
    DecayStats mDecaysStats;
    
    std::fstream mDatafile;

    G4String mHitsFileName;
    G4String mDecaysFileName;
    G4String mMultiplicityFileName;

    OMSimDecaysAnalysis() = default;
    ~OMSimDecaysAnalysis() = default;
    OMSimDecaysAnalysis(const OMSimDecaysAnalysis &) = delete;
    OMSimDecaysAnalysis &operator=(const OMSimDecaysAnalysis &) = delete;
};

#endif
