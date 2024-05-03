#include "OMSimEffectiveAreaAnalyisis.hh"
#include "OMSimCommandArgsTable.hh"
#include "OMSimHitManager.hh"

/**
 * @brief Writes the header line to the output file.
 */
void OMSimEffectiveAreaAnalyisis::writeHeader()
{
	mDatafile.open(mOutputFileName.c_str(), std::ios::out | std::ios::app);
	mDatafile << "# Phi(deg)"
			  << "\t"
			  << "Theta(deg)"
			  << "\t"
			  << "hits[1perPMT]"
			  << "\t"
			  << "total_hits"
			//  << "\t"
			//  << "EA_Total(cm^2)"
			//  << "\t"
			//  << "EA_Total_error(cm^2)"
			  << "\t" << G4endl;
	mDatafile.close();
}

/**
 * @brief Calculates the effective area based on the number of hits and beam properties.
 * @param pHits The number of hits.
 * @return Returns a structure with the effective area and its uncertainty.
 */
effectiveAreaResult OMSimEffectiveAreaAnalyisis::calculateEffectiveArea(double pHits)
{
	OMSimCommandArgsTable &lArgs = OMSimCommandArgsTable::getInstance();
	G4double lBeamRadius = lArgs.get<G4double>("radius") / 10.; // in cm
	G4double lNumberPhotons = lArgs.get<int>("numevents");
	G4double lBeamArea = CLHEP::pi * lBeamRadius * lBeamRadius;
	G4double lEA = pHits * lBeamArea / lNumberPhotons;
	G4double lEAError = sqrt(pHits) * lBeamArea / lNumberPhotons;
	return {lEA, lEAError};
}
/**
 * @brief Writes a scan result to the output file.
 * @param pPhi The phi angle used in the scan to be written to the output file.
 * @param pTheta The theta angle used in the scan to be written to the output file.
 * @param pTimeWindow The time window for multiplicity calculation.
 */
void OMSimEffectiveAreaAnalyisis::writeScan(G4double pPhi, G4double pTheta)
{
	std::vector<double> lHits = OMSimHitManager::getInstance().countHits();
	mDatafile.open(mOutputFileName.c_str(), std::ios::out | std::ios::app);
	mDatafile << pPhi << "\t" << pTheta << "\t";
	
	G4double lTotalHits = 0;
	for (const auto &hit : lHits)
	{
		mDatafile << hit << "\t";
		lTotalHits = hit; // last element is total nr of hits
	}
	//effectiveAreaResult lEffectiveArea = calculateEffectiveArea(lTotalHits);

	//mDatafile << lEffectiveArea.EA << "\t" << lEffectiveArea.EAError << "\t";
	mDatafile << G4endl;
	mDatafile.close();
}

/**
 * @brief Calls calculateMultiplicity and writes the results to the output file.
 */
void OMSimEffectiveAreaAnalyisis::writeMultiplicity(G4double pTimeWindow)
{
	HitStats lHits = OMSimHitManager::getInstance().getHitsOfModule();

	if (lHits.event_id.size() > 0)
	{
		std::vector<int> lMultiplicity = OMSimHitManager::getInstance().calculateMultiplicity(pTimeWindow);
		mDatafile.open(mOutputFileNameMultiplicity.c_str(), std::ios::out | std::ios::app);
		for (const auto &value : lMultiplicity)
		{
			mDatafile << value << "\t";
		}
		mDatafile << G4endl;
		mDatafile.close();
	}
}

/**
 * @brief Write data of the hits to the output file.
 */
void OMSimEffectiveAreaAnalyisis::writeHitInformation()
{
	HitStats lHits = OMSimHitManager::getInstance().getHitsOfModule();

	if (lHits.event_id.size() > 0)
	{
		OMSimHitManager::getInstance().sortHitStatsByTime(lHits);
		mDatafile.open(mOutputFileNameInfo.c_str(), std::ios::out | std::ios::app);
		for (int i = 0; i < (int)lHits.event_id.size(); i++)
		{
			// mDatafile << lHits.event_id.at(i) << "\t";
			// mDatafile << std::setprecision(13);
			mDatafile << lHits.hit_time.at(i) / ns << "\t";
			// mDatafile << std::setprecision(4);
			mDatafile << lHits.PMT_hit.at(i) << "\t";
			/*	mDatafile << lHits.photon_energy.at(i) << "\t";
				mDatafile << lHits.photon_global_position.at(i).x() << "\t";
				mDatafile << lHits.photon_global_position.at(i).y() << "\t";
				mDatafile << lHits.photon_global_position.at(i).z() << "\t";
				mDatafile << lHits.PMT_response.at(i).PE << "\t";
				mDatafile << lHits.PMT_response.at(i).TransitTime << "\t";
				mDatafile << lHits.PMT_response.at(i).DetectionProbability << "\t"; */
			mDatafile << G4endl;
		}

		mDatafile << G4endl;
		mDatafile.close();
	}
}