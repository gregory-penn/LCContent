/**
 *  @file   LCContent/src/LCCheating/CheatingTrackToClusterMatching.cc
 *
 *  @brief  Implementation of the cheating track to cluster matching algorithm class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCCheating/CheatingTrackToClusterMatching.h"

using namespace pandora;

namespace lc_content {

CheatingTrackToClusterMatching::CheatingTrackToClusterMatching() : m_shouldSkipTracksThatCannotFormPfo(false) {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingTrackToClusterMatching::Run() {
  // Read current lists
  const TrackList* pCurrentTrackList = NULL;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCurrentTrackList));

  const ClusterList* pClusterList = NULL;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

  // Clear any existing track - cluster associations
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveCurrentTrackClusterAssociations(*this));

  // Construct a map from mc particle to tracks
  typedef std::map<const MCParticle*, TrackList> TracksPerMCParticle;
  TracksPerMCParticle tracksPerMCParticle;
  unsigned int nSkippedTracks(0); // tracks left unassociated by ShouldSkipTracksThatCannotFormPfo

  for (TrackList::const_iterator iter = pCurrentTrackList->begin(), iterEnd = pCurrentTrackList->end(); iter != iterEnd;
       ++iter) {
    try {
      const Track* const pTrack = *iter;

      // A track that TrackPreparation will not use for a pfo must not claim a cluster; its cluster then stays
      // unassociated and PfoCreation can make a neutral pfo from it
      if (m_shouldSkipTracksThatCannotFormPfo && !pTrack->CanFormPfo() && !pTrack->CanFormClusterlessPfo()) {
        ++nSkippedTracks;

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
          std::cout << "CheatingTrackToClusterMatching: track skipped (CanFormPfo 0, CanFormClusterlessPfo 0), p at dca "
                    << pTrack->GetMomentumAtDca().GetMagnitude() << " GeV" << std::endl;

        continue;
      }

      const MCParticle* const pMCParticle(MCParticleHelper::GetMainMCParticle(pTrack));

      TracksPerMCParticle::iterator itTracksPerMCParticle(tracksPerMCParticle.find(pMCParticle));

      if (tracksPerMCParticle.end() == itTracksPerMCParticle) {
        if (!tracksPerMCParticle.insert(TracksPerMCParticle::value_type(pMCParticle, TrackList(1, pTrack))).second)
          throw StatusCodeException(STATUS_CODE_FAILURE);
      } else {
        itTracksPerMCParticle->second.push_back(pTrack);
      }
    } catch (StatusCodeException&) {
    }
  }

  // Construct a map from mc particle to clusters
  typedef std::map<const MCParticle*, ClusterList> ClustersPerMCParticle;
  ClustersPerMCParticle clustersPerMCParticle;

  for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd;
       ++iter) {
    try {
      const Cluster* const pCluster = *iter;
      const MCParticle* const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

      ClustersPerMCParticle::iterator itClustersPerMCParticle(clustersPerMCParticle.find(pMCParticle));

      if (clustersPerMCParticle.end() == itClustersPerMCParticle) {
        if (!clustersPerMCParticle.insert(ClustersPerMCParticle::value_type(pMCParticle, ClusterList(1, pCluster)))
                 .second)
          throw StatusCodeException(STATUS_CODE_FAILURE);
      } else {
        itClustersPerMCParticle->second.push_back(pCluster);
      }
    } catch (StatusCodeException&) {
    }
  }

  // Make the track to cluster associations
  MCParticleList mcParticleList;
  for (const auto& mapEntry : tracksPerMCParticle)
    mcParticleList.push_back(mapEntry.first);
  mcParticleList.sort(PointerLessThan<MCParticle>());

  for (const MCParticle* const pMCParticle : mcParticleList) {
    const TrackList& trackList(tracksPerMCParticle.at(pMCParticle));

    ClustersPerMCParticle::const_iterator itClustersPerMCParticle(clustersPerMCParticle.find(pMCParticle));

    if (clustersPerMCParticle.end() == itClustersPerMCParticle)
      continue;

    const ClusterList& clusterList = itClustersPerMCParticle->second;

    if (trackList.empty() || clusterList.empty())
      continue;

    for (TrackList::const_iterator itTrack = trackList.begin(), itTrackEnd = trackList.end(); itTrack != itTrackEnd;
         ++itTrack) {
      // If the mc particle is associated with multiple clusters, can only associate to highest energy cluster (clusters
      // should be merged)
      const Cluster* pHighestEnergyCluster = NULL;
      float highestEnergy(-std::numeric_limits<float>::max());

      for (ClusterList::const_iterator itCluster = clusterList.begin(), itClusterEnd = clusterList.end();
           itCluster != itClusterEnd; ++itCluster) {
        const float clusterEnergy((*itCluster)->GetHadronicEnergy());

        if (clusterEnergy > highestEnergy) {
          highestEnergy = clusterEnergy;
          pHighestEnergyCluster = *itCluster;
        }
      }

      if (NULL == pHighestEnergyCluster)
        throw StatusCodeException(STATUS_CODE_FAILURE);

      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                               PandoraContentApi::AddTrackClusterAssociation(*this, *itTrack, pHighestEnergyCluster));
    }
  }

  if (m_shouldSkipTracksThatCannotFormPfo && PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
    std::cout << "CheatingTrackToClusterMatching: " << nSkippedTracks << " of " << pCurrentTrackList->size()
              << " tracks left unassociated (ShouldSkipTracksThatCannotFormPfo)" << std::endl;

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingTrackToClusterMatching::ReadSettings(const TiXmlHandle xmlHandle) {
  // Optional, default false (unchanged behavior when the xml does not set it)
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "ShouldSkipTracksThatCannotFormPfo", m_shouldSkipTracksThatCannotFormPfo));

  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
