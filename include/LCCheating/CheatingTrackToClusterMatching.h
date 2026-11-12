/**
 *  @file   LCContent/include/LCCheating/CheatingTrackToClusterMatching.h
 *
 *  @brief  Header file for the cheating track to cluster matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_CHEATING_TRACK_TO_CLUSTER_MATCHING_ALGORITHM_H
#define LC_CHEATING_TRACK_TO_CLUSTER_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lc_content {

/**
 *  @brief CheatingTrackToClusterMatching class
 */
class CheatingTrackToClusterMatching : public pandora::Algorithm {
public:
  /**
   *  @brief  Default constructor
   */
  CheatingTrackToClusterMatching();

private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  // Opt-in. A track that can form neither a charged pfo with a cluster (CanFormPfo) nor a
  // clusterless one (CanFormClusterlessPfo) is dropped by TrackPreparation, but an association made here still blocks
  // its cluster from becoming a neutral pfo (PfoCreation only uses clusters without associated tracks), so the cluster
  // is lost. Example: photon conversions in the tracker (first hit beyond MaxBarrelTrackerInnerRDistance).
  bool m_shouldSkipTracksThatCannotFormPfo; ///< Whether to leave such tracks unassociated (default false: old behavior)
};

} // namespace lc_content

#endif // #ifndef LC_CHEATING_TRACK_TO_CLUSTER_MATCHING_ALGORITHM_H
