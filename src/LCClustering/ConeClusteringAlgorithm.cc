/**
 *  @file   LCContent/src/LCClustering/ConeClusteringAlgorithm.cc
 * 
 *  @brief  Implementation of the clustering algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCClustering/ConeClusteringAlgorithm.h"

#include "LCHelpers/SortingHelper.h"

#include <list>

using namespace pandora;

namespace lc_content
{

ConeClusteringAlgorithm::ConeClusteringAlgorithm() :
    m_clusterSeedStrategy(2),
    m_shouldUseOnlyECalHits(false),
    m_shouldUseIsolatedHits(false),
    m_layersToStepBackFine(3),
    m_layersToStepBackCoarse(3),
    m_clusterFormationStrategy(0),
    m_genericDistanceCut(1.f),
    m_minHitTrackCosAngle(0.f),
    m_minHitClusterCosAngle(0.f),
    m_shouldUseTrackSeed(true),
    m_trackSeedCutOffLayer(0),
    m_shouldFollowInitialDirection(false),
    m_sameLayerPadWidthsFine(2.8f),
    m_sameLayerPadWidthsCoarse(1.8f),
    m_coneApproachMaxSeparation2(1000.f * 1000.f),
    m_tanConeAngleFine(0.3f),
    m_tanConeAngleCoarse(0.5f),
    m_additionalPadWidthsFine(2.5f),
    m_additionalPadWidthsCoarse(2.5f),
    m_maxClusterDirProjection(200.f),
    m_minClusterDirProjection(-10.f),
    m_trackPathWidth(2.f),
    m_maxTrackSeedSeparation2(250.f * 250.f),
    m_maxLayersToTrackSeed(3),
    m_maxLayersToTrackLikeHit(3),
    m_nLayersSpannedForFit(6),
    m_nLayersSpannedForApproxFit(10),
    m_nLayersToFit(8),
    m_nLayersToFitLowMipCut(0.5f),
    m_nLayersToFitLowMipMultiplier(2),
    m_fitSuccessDotProductCut1(0.75f),
    m_fitSuccessChi2Cut1(5.0f),
    m_fitSuccessDotProductCut2(0.50f),
    m_fitSuccessChi2Cut2(2.5f),
    m_mipTrackChi2Cut(2.5f),
    m_firstLayer(1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::Run()
{

    std::cout << "In Run() of ConeClusteringAlgorithm" << std::endl;

    m_firstLayer = (PandoraContentApi::GetPlugins(*this)->GetPseudoLayerPlugin()->GetPseudoLayerAtIp());

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    if (pCaloHitList->empty()) {
        std::cout << "GP ERROR: Calo hit list is empty! No clustering done in this call." << std::endl;
        return STATUS_CODE_SUCCESS;
    }
    const TrackList *pTrackList = nullptr;
    if (0 != m_clusterSeedStrategy) {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pTrackList));
    }

    this->InitializeKDTrees(pTrackList,pCaloHitList);

    OrderedCaloHitList orderedCaloHitList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(*pCaloHitList));

    ClusterVector clusterVector;
    m_tracksToClusters.clear();
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SeedClustersWithTracks(pTrackList,clusterVector));

    // do the clustering
    m_hitsToClusters.clear();
    // I think first ordered by pseudoLayer
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        const unsigned int pseudoLayer(iter->first);
        // std::cout << "PseudoLayer: " << pseudoLayer << std::endl;
        CaloHitVector relevantCaloHits;

        // then ordered by hits in the layer, I assume
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

            if ((m_shouldUseIsolatedHits || !pCaloHit->IsIsolated()) &&
                (!m_shouldUseOnlyECalHits || (ECAL == pCaloHit->GetHitType())) &&
                (PandoraContentApi::IsAvailable(*this, pCaloHit)))
            {
                relevantCaloHits.push_back(pCaloHit);
                // std::cout << "Found a relevant calo hit!" << std::endl;
            }
            // else {
            //     std::cout << "Discarding calo hit." << std::endl;
            // }
        }

        ClusterFitResultMap clusterFitResultMap;
        // below seems to be an intiial clustering fit, without any consideration of tracks, filled to clusterFitResultMap
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetCurrentClusterFitResults(clusterVector, clusterFitResultMap));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FindHitsInPreviousLayers(pseudoLayer, relevantCaloHits, clusterFitResultMap, clusterVector));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FindHitsInSameLayer(pseudoLayer, relevantCaloHits, clusterFitResultMap, clusterVector));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RemoveEmptyClusters(clusterVector));

    //reset our kd trees and maps if everything turned out well
    m_tracksKdTree.clear();
    m_hitsKdTree.clear();
    m_hitsToClusters.clear();
    m_tracksToClusters.clear();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::InitializeKDTrees(const TrackList *const pTrackList, const CaloHitList *const pCaloHitList)
{
    // load the kd-tree of tracks that we will use
    m_tracksKdTree.clear();
    m_trackNodes.clear();

    if (nullptr != pTrackList)
    {
        KDTreeCube tracksBoundingRegion = fill_and_bound_3d_kd_tree(this, *pTrackList, m_trackNodes);
        m_tracksKdTree.build(m_trackNodes,tracksBoundingRegion);
        m_trackNodes.clear();
    }

    // make sure the hit kd tree is ready
    m_hitsKdTree.clear();
    m_hitNodes.clear();
    KDTreeTesseract hitsBoundingRegion = fill_and_bound_4d_kd_tree(this,*pCaloHitList,m_hitNodes);
    m_hitsKdTree.build(m_hitNodes,hitsBoundingRegion);
    m_hitNodes.clear();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::SeedClustersWithTracks(const TrackList *const pTrackList, ClusterVector &clusterVector)
{

    // std::cout << "In SeedClustersWithTracks" << std::endl;

    if (0 == m_clusterSeedStrategy)
        return STATUS_CODE_SUCCESS;

    int canFormPFOInEvent = 0;

    // if we are known to be seeding with tracks we must have a track list
    if (nullptr == pTrackList)
        return STATUS_CODE_FAILURE;

    for (TrackList::const_iterator iter = pTrackList->begin(), iterEnd = pTrackList->end(); iter != iterEnd; ++iter)
    {
        const Track *const pTrack = *iter;

        if (!pTrack->CanFormPfo())
            continue;

        canFormPFOInEvent += 1;

        bool useTrack(false);

        if (2 == m_clusterSeedStrategy)
        {
            useTrack = true;
        }
        else if ((1 == m_clusterSeedStrategy) && pTrack->IsProjectedToEndCap())
        {
            useTrack = true;
        }

        if (useTrack)
        {
            const Cluster *pCluster = nullptr;
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_pTrack = pTrack;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
            clusterVector.push_back(pCluster);
            m_tracksToClusters.emplace(pTrack,pCluster);
        }
        else { 
            std::cout << "GP ERROR: useTrack = false! This shouldn't be happening with cluster seed strategy = 2" << std::endl; // I never see this. Good!
        }
    }

    std::cout << "Number of tracks that can form a PFO found: " << canFormPFOInEvent << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::GetCurrentClusterFitResults(const ClusterVector &clusterVector, ClusterFitResultMap &clusterFitResultMap) const
{

    // std::cout << "In GetCurrentClusterFitResults" << std::endl;

    if (!clusterFitResultMap.empty()) {
        std::cout << "GP ERROR: Found something in clusterFitResultMap! Not sure what will happen now." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        // std::cout << "Looping over cluster! This shouldn't be happening if there is no track that can form a PFO!" << std::endl; warning, this prints out many many times
        const Cluster *const pCluster = *iter;
        ClusterFitResult clusterFitResult;

        if (pCluster->GetNCaloHits() > 1)
        {
            const unsigned int innerLayer(pCluster->GetInnerPseudoLayer());
            const unsigned int outerLayer(pCluster->GetOuterPseudoLayer());
            const unsigned int nLayersSpanned(outerLayer - innerLayer);
            // std::cout << "number of layers spanned: " << nLayersSpanned << std::endl;

            if (nLayersSpanned > m_nLayersSpannedForFit) // default value: 6 layers. Not sure if this is small or large.
            {
                unsigned int nLayersToFit(m_nLayersToFit);

                if (pCluster->GetMipFraction() - m_nLayersToFitLowMipCut < std::numeric_limits<float>::epsilon())
                    nLayersToFit *= m_nLayersToFitLowMipMultiplier; // muliply the number of fit layers by two if the MIP fraction is low. Default cut is @ 0.5
                    // MIP Fraction: fraction of energy compared to some MIP value? Why would we fit more layers for lower MIP, rather than less?

                const unsigned int startLayer( (nLayersSpanned > nLayersToFit) ? (outerLayer - nLayersToFit) : innerLayer);
                (void) ClusterFitHelper::FitLayerCentroids(pCluster, startLayer, outerLayer, clusterFitResult);

                if (clusterFitResult.IsFitSuccessful())
                {
                    const float dotProduct(clusterFitResult.GetDirection().GetDotProduct(pCluster->GetInitialDirection()));
                    const float chi2(clusterFitResult.GetChi2());

                    if (((dotProduct < m_fitSuccessDotProductCut1) && (chi2 > m_fitSuccessChi2Cut1)) ||
                        ((dotProduct < m_fitSuccessDotProductCut2) && (chi2 > m_fitSuccessChi2Cut2)) )
                    {
                        clusterFitResult.SetSuccessFlag(false);
                    }
                }
            }
            else if (nLayersSpanned > m_nLayersSpannedForApproxFit)
            {
                const CartesianVector centroidChange(pCluster->GetCentroid(outerLayer) - pCluster->GetCentroid(innerLayer));
                clusterFitResult.Reset();
                clusterFitResult.SetDirection(centroidChange.GetUnitVector());
                clusterFitResult.SetSuccessFlag(true);
            }
        }

        if (!clusterFitResultMap.insert(ClusterFitResultMap::value_type(pCluster, clusterFitResult)).second)
            return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// inprogress
StatusCode ConeClusteringAlgorithm::FindHitsInPreviousLayers(unsigned int pseudoLayer, const CaloHitVector &relevantCaloHits,
    const ClusterFitResultMap &clusterFitResultMap, ClusterVector & /*clusterVector*/)
{
    // TODO: this seems to change when ConeClustering is recalled after a successful run. BUG?
    const float maxTrackSeedSeparation = std::sqrt(m_maxTrackSeedSeparation2);

    // std::cout << "maxTrackSeedSeparation: " << maxTrackSeedSeparation << std::endl; 

    std::vector<HitKDNode> found_hits;
    std::vector<TrackKDNode> found_tracks;
    ClusterSet nearby_clusters;

    for (const CaloHit *const pCaloHit : relevantCaloHits) // loop over relevant Calo hits. What happened to the initially fit clusters?? Starting from hits again?
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit)) {
            std::cout << "GP ERROR: Hit is not available! Not sure what will happen." << std::endl;
            continue;
        }
        const float additionalPadWidths = ((PandoraContentApi::GetGeometry(*this)->GetHitTypeGranularity(pCaloHit->GetHitType()) <= FINE) ?
            m_additionalPadWidthsFine * pCaloHit->GetCellLengthScale() : m_additionalPadWidthsCoarse * pCaloHit->GetCellLengthScale());
        const float largestAllowedDistanceForSearch = std::max(maxTrackSeedSeparation, m_maxClusterDirProjection + additionalPadWidths);
        // maxTrackSeedSeparation default: 250 mm (unit?), m_maxClusterDirProjection default: 200 (mm?) 

        if (largestAllowedDistanceForSearch != maxTrackSeedSeparation) {
            // std::cout << "GP MESSAGE: additional pad widths are relevant - the search cone is now " << largestAllowedDistanceForSearch << std::endl;
        }

        const Cluster *pBestCluster = nullptr;
        float bestClusterEnergy(0.f);
        float smallestGenericDistance(m_genericDistanceCut); // looks to be 1 mm?
        const unsigned int layersToStepBack((PandoraContentApi::GetGeometry(*this)->GetHitTypeGranularity(pCaloHit->GetHitType()) <= FINE) ?
            m_layersToStepBackFine : m_layersToStepBackCoarse); // both 3 layers by default

        // Associate with existing clusters in stepBack layers. If stepBackLayer == pseudoLayer, will examine track projections
        for (unsigned int stepBackLayer = 1; (stepBackLayer <= layersToStepBack) && (stepBackLayer <= pseudoLayer); ++stepBackLayer)
        {
            const unsigned int searchLayer(pseudoLayer - stepBackLayer);

            // need to reorganize this to use a kd-tree. On rechits comprising clusters we are mutating
            // goal -> determine search distances for KD-tree from cut values and associated scalings
            // search for tracks that would satisfy the search criteria in GetGenericDistanceToHit()
            KDTreeCube searchRegionTks = build_3d_kd_search_region(pCaloHit, largestAllowedDistanceForSearch, largestAllowedDistanceForSearch, largestAllowedDistanceForSearch);
            m_tracksKdTree.search(searchRegionTks,found_tracks); // I suppose this finds tracks
            int nFoundTrks = 0;
            for (auto &track : found_tracks )
            {
                nFoundTrks +=1 ;
                // std::cout << "Found a track. Number times it's been hit for this relevant calo hit: " << nFoundTrks << std::endl;
                auto assc_cluster = m_tracksToClusters.find(track.data);
                // below skips if there is no cluster associated to that track
                if (assc_cluster != m_tracksToClusters.end())
                {
                    nearby_clusters.insert(assc_cluster->second);
                }
            }
            found_tracks.clear();

            // now search for hits-in-clusters that would also satisfy the criteria
            KDTreeTesseract searchRegionHits = build_4d_kd_search_region(pCaloHit, largestAllowedDistanceForSearch, largestAllowedDistanceForSearch, largestAllowedDistanceForSearch, searchLayer);
            m_hitsKdTree.search(searchRegionHits,found_hits);
            for (auto &hit : found_hits)
            {
                auto assc_cluster = m_hitsToClusters.find(hit.data);
                if (assc_cluster != m_hitsToClusters.end())
                {
                    nearby_clusters.insert(assc_cluster->second);
                }
            }
            found_hits.clear();

            ClusterList nearbyClusterList(nearby_clusters.begin(), nearby_clusters.end());
            nearbyClusterList.sort(SortingHelper::SortClustersByNHits);
            nearby_clusters.clear();

            // Instead of using the full cluster list we use only those clusters that are found to be nearby according to the KD-tree
            // ---- This can be optimized further for sure. (for instance having a match by KD-tree qualifies a ton of the loops later
            // See if hit should be associated with any existing clusters
            for (const Cluster *const pCluster : nearbyClusterList)
            {
                float genericDistance(std::numeric_limits<float>::max());
                const float clusterEnergy(pCluster->GetHadronicEnergy());

                PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_UNCHANGED, !=, this->GetGenericDistanceToHit(pCluster,
                    pCaloHit, searchLayer, clusterFitResultMap, genericDistance));

                // std::cout << "This is the distance now: " << genericDistance << std::endl;

                if ((genericDistance < smallestGenericDistance) ||
                    ((genericDistance == smallestGenericDistance) && (clusterEnergy > bestClusterEnergy)))
                {
                    pBestCluster = pCluster;
                    bestClusterEnergy = clusterEnergy;
                    smallestGenericDistance = genericDistance;
                }
            }

            // Add best hit found after completing examination of a stepback layer
            if ((0 == m_clusterFormationStrategy) && (nullptr != pBestCluster))
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pBestCluster, pCaloHit));
                m_hitsToClusters.emplace(pCaloHit, pBestCluster);
                break;
            }
        }

        // Add best hit found after examining all stepback layers
        if ((1 == m_clusterFormationStrategy) && (nullptr != pBestCluster))
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pBestCluster, pCaloHit));
            m_hitsToClusters.emplace(pCaloHit, pBestCluster);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::FindHitsInSameLayer(unsigned int pseudoLayer, const CaloHitVector &relevantCaloHits,
    const ClusterFitResultMap &clusterFitResultMap, ClusterVector &clusterVector)
{
    const float maxTrackSeedSeparation = std::sqrt(m_maxTrackSeedSeparation2);

    // std::cout << "Printing maxTrackSeedSeparation: " << maxTrackSeedSeparation << std::endl;

    //keep a list of available hits with the most energetic available hit at the back
    std::list<unsigned> available_hits_in_layer;
    for (unsigned i = 0; i < relevantCaloHits.size(); ++i )
    {
        if (!PandoraContentApi::IsAvailable(*this, relevantCaloHits[i]))
            continue;
        available_hits_in_layer.push_back(i);
    }

    //tactical cache hits -> tracks
    std::unordered_multimap<const CaloHit*, const Track*> hitsToTracksLocal;

    //tactical cache hits -> hits
    std::unordered_multimap<const CaloHit*, const CaloHit*> hitsToHitsLocal;

    //pull out result caches to help keep memory locality
    std::vector<TrackKDNode> found_tracks;
    std::vector<HitKDNode> found_hits;

    ClusterSet nearby_clusters;

    while (!available_hits_in_layer.empty())
    {
        bool clustersModified = true;

        while (clustersModified)
        {
            clustersModified = false;

            for (auto iter = available_hits_in_layer.begin(); iter != available_hits_in_layer.end();)
            {
                // this his is assured to be usable by the lines above and algorithm course
                const CaloHit *const pCaloHit = relevantCaloHits[*iter];

                const float pad_search_width = ((PandoraContentApi::GetGeometry(*this)->GetHitTypeGranularity(pCaloHit->GetHitType()) <= FINE) ?
                    (m_sameLayerPadWidthsFine * pCaloHit->GetCellLengthScale()) :
                    (m_sameLayerPadWidthsCoarse * pCaloHit->GetCellLengthScale()) );

                const float track_search_width = maxTrackSeedSeparation;
                const float hit_search_width = pad_search_width;

                // std::cout << "Printing hit_search_width: " << hit_search_width << std::endl;

                const Cluster *pBestCluster = nullptr;
                float bestClusterEnergy(0.f);
                float smallestGenericDistance(m_genericDistanceCut);

                // search for tracks that would satisfy the search criteria in GetGenericDistanceToHit()
                auto track_assc_cache = hitsToTracksLocal.find(pCaloHit);
                if (track_assc_cache != hitsToTracksLocal.end())
                {
                    auto range = hitsToTracksLocal.equal_range(pCaloHit);
                    for (auto itr = range.first; itr != range.second; ++itr)
                    {
                        auto assc_cluster = m_tracksToClusters.find(itr->second);
                        if(assc_cluster != m_tracksToClusters.end())
                        {
                            nearby_clusters.insert(assc_cluster->second);
                        }
                    }
                }
                else
                {
                    KDTreeCube searchRegionTks = build_3d_kd_search_region(pCaloHit, track_search_width, track_search_width, track_search_width);
                    m_tracksKdTree.search(searchRegionTks,found_tracks);
                    for (auto &track : found_tracks)
                    {
                        hitsToTracksLocal.emplace(pCaloHit, track.data);
                        auto assc_cluster = m_tracksToClusters.find(track.data);
                        if (assc_cluster != m_tracksToClusters.end())
                        {
                            nearby_clusters.insert(assc_cluster->second);
                        }
                    }
                    found_tracks.clear();
                }
                // now search for hits-in-clusters that would also satisfy the criteria
                auto hits_assc_cache = hitsToHitsLocal.find(pCaloHit);
                if (hits_assc_cache != hitsToHitsLocal.end())
                {
                    auto range = hitsToHitsLocal.equal_range(pCaloHit);
                    for (auto itr = range.first; itr != range.second; ++itr)
                    {
                        auto assc_cluster = m_hitsToClusters.find(itr->second);
                        if( assc_cluster != m_hitsToClusters.end() )
                        {
                            nearby_clusters.insert(assc_cluster->second);
                        }
                    }
                }
                else
                {
                    KDTreeTesseract searchRegionHits = build_4d_kd_search_region(pCaloHit, hit_search_width, hit_search_width, hit_search_width, pseudoLayer);
                    m_hitsKdTree.search(searchRegionHits,found_hits);
                    for (auto &hit : found_hits)
                    {
                        hitsToHitsLocal.emplace(pCaloHit, hit.data);
                        auto assc_cluster = m_hitsToClusters.find(hit.data);
                        if (assc_cluster != m_hitsToClusters.end())
                        {
                            nearby_clusters.insert(assc_cluster->second);
                        }
                    }
                    found_hits.clear();
                }

                ClusterList nearbyClusterList(nearby_clusters.begin(), nearby_clusters.end());
                nearbyClusterList.sort(SortingHelper::SortClustersByNHits);
                nearby_clusters.clear();

                // See if hit should be associated with any existing clusters
                for (ClusterList::iterator clusterIter = nearbyClusterList.begin(), clusterIterEnd = nearbyClusterList.end();
                    clusterIter != clusterIterEnd; ++clusterIter)
                {
                    const Cluster *const pCluster = *clusterIter;
                    float genericDistance(std::numeric_limits<float>::max());
                    const float clusterEnergy(pCluster->GetHadronicEnergy());

                    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_UNCHANGED, !=, this->GetGenericDistanceToHit(pCluster,
                        pCaloHit, pseudoLayer, clusterFitResultMap, genericDistance));

                    if ((genericDistance < smallestGenericDistance) || ((genericDistance == smallestGenericDistance) && (clusterEnergy > bestClusterEnergy)))
                    {
                        pBestCluster = pCluster;
                        bestClusterEnergy = clusterEnergy;
                        smallestGenericDistance = genericDistance;
                    }
                }

                if (nullptr != pBestCluster)
                {
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pBestCluster, pCaloHit));
                    m_hitsToClusters.emplace(pCaloHit, pBestCluster);
                    // remove this hit and advance to the next
                    iter = available_hits_in_layer.erase(iter);
                    clustersModified = true;
                }
                else
                { // otherwise advance the iterator
                    ++iter;
                }
            }
        }

        // If there is no cluster within the search radius, seed a new cluster with this hit
        if (!available_hits_in_layer.empty())
        {
            unsigned index = *(available_hits_in_layer.begin());
            available_hits_in_layer.pop_front();
            const CaloHit *const  pCaloHit = relevantCaloHits[index];
            // hit is assured to be valid
            const Cluster *pCluster = nullptr;
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.push_back(pCaloHit);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
            clusterVector.push_back(pCluster);
            m_hitsToClusters.emplace(pCaloHit,pCluster);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::GetGenericDistanceToHit(const Cluster *const pCluster, const CaloHit *const pCaloHit, const unsigned int searchLayer,
    const ClusterFitResultMap &clusterFitResultMap, float &genericDistance) const
{
    const unsigned int firstLayer = m_firstLayer;

    // Use position of track projection at calorimeter. Proceed only if projection is reasonably compatible with calo hit
    if (((searchLayer == 0) || (searchLayer < firstLayer)) && pCluster->IsTrackSeeded())
    {
        const TrackState &trackState(pCluster->GetTrackSeed()->GetTrackStateAtCalorimeter());
        const CartesianVector trackDirection(trackState.GetMomentum().GetUnitVector());

        if (pCaloHit->GetExpectedDirection().GetCosOpeningAngle(trackDirection) < m_minHitTrackCosAngle) {
            std::cout << "GP ERROR: Cluster skipped due to cos opening angle of track being inconsistent. This shouldn't be happening (default minimum requirement is 0.)" << std::endl;
            return STATUS_CODE_UNCHANGED;
        }
        return this->GetConeApproachDistanceToHit(pCaloHit, trackState.GetPosition(), trackDirection, genericDistance);
    }

    // std::cout << "How often is this actually getting past the first if?" << std::endl;

    // Check that cluster is occupied in the searchlayer and is reasonably compatible with calo hit
    OrderedCaloHitList::const_iterator clusterHitListIter = pCluster->GetOrderedCaloHitList().find(searchLayer);

    if (pCluster->GetOrderedCaloHitList().end() == clusterHitListIter)
        return STATUS_CODE_UNCHANGED;

    ClusterFitResultMap::const_iterator fitResultIter(clusterFitResultMap.find(pCluster));
    const ClusterFitResult *const clusterFitResult((clusterFitResultMap.end() != fitResultIter) ? &(fitResultIter->second) : nullptr);
    const CartesianVector &clusterDirection( (clusterFitResult != nullptr && clusterFitResult->IsFitSuccessful()) ?  
         clusterFitResult->GetDirection() : pCluster->GetInitialDirection());

    if (pCaloHit->GetExpectedDirection().GetCosOpeningAngle(clusterDirection) < m_minHitClusterCosAngle)
        return STATUS_CODE_UNCHANGED;

    // Cone approach measurements
    float initialDirectionDistance(std::numeric_limits<float>::max());
    float currentDirectionDistance(std::numeric_limits<float>::max());
    float trackSeedDistance(std::numeric_limits<float>::max());

    const bool useTrackSeed(m_shouldUseTrackSeed && pCluster->IsTrackSeeded());
    const bool followInitialDirection(m_shouldFollowInitialDirection && pCluster->IsTrackSeeded() && (searchLayer > m_trackSeedCutOffLayer));

    if (!useTrackSeed || (searchLayer > m_trackSeedCutOffLayer))
    {
        const CaloHitList *pClusterCaloHitList = clusterHitListIter->second;

        if (searchLayer == pCaloHit->GetPseudoLayer())
        {
            return this->GetDistanceToHitInSameLayer(pCaloHit, pClusterCaloHitList, genericDistance);
        }

        // Measurement using initial cluster direction
        StatusCode initialStatusCode = this->GetConeApproachDistanceToHit(pCaloHit, pClusterCaloHitList, pCluster->GetInitialDirection(),
            initialDirectionDistance);

        if (STATUS_CODE_SUCCESS == initialStatusCode)
        {
            if (followInitialDirection)
                initialDirectionDistance /= 5.;
        }
        else if (STATUS_CODE_UNCHANGED != initialStatusCode)
        {
            return initialStatusCode;
        }

        // Measurement using current cluster direction
        if (clusterFitResult != nullptr && clusterFitResult->IsFitSuccessful())
        {
            StatusCode currentStatusCode = this->GetConeApproachDistanceToHit(pCaloHit, pClusterCaloHitList, clusterFitResult->GetDirection(),
                currentDirectionDistance);

            if (STATUS_CODE_SUCCESS == currentStatusCode)
            {
                if ((currentDirectionDistance < m_genericDistanceCut) && pCluster->IsTrackSeeded())
                {
                    try
                    {
                        if (clusterFitResult->GetChi2() < m_mipTrackChi2Cut)
                            currentDirectionDistance /= 5.;
                    }
                    catch (StatusCodeException &) {}
                }
            }
            else if (STATUS_CODE_UNCHANGED != currentStatusCode)
            {
                return currentStatusCode;
            }
        }
    }

    // Seed track distance measurements
    if (useTrackSeed && !followInitialDirection)
    {
        // std::cout << "How often are we getting to this track seed?" << std::endl;
        StatusCode trackStatusCode = this->GetDistanceToTrackSeed(pCluster, pCaloHit, searchLayer, trackSeedDistance);
        // std::cout << "Track seed distance: " << trackSeedDistance << std::endl;
        if (STATUS_CODE_SUCCESS == trackStatusCode)
        {
            if (trackSeedDistance < m_genericDistanceCut)
                trackSeedDistance /= 5.;
        }
        else if (STATUS_CODE_UNCHANGED != trackStatusCode)
        {
            return trackStatusCode;
        }
    }

    if (std::min(trackSeedDistance, std::min(initialDirectionDistance, currentDirectionDistance)) ==  trackSeedDistance) {
        // std::cout << "trackSeedDistance is being returned! This is probably good news?" << std::endl;
    }

    // Identify best measurement of generic distance
    const float smallestDistance(std::min(trackSeedDistance, std::min(initialDirectionDistance, currentDirectionDistance)));
    if (smallestDistance < genericDistance)
    {
        genericDistance = smallestDistance;

        if (std::numeric_limits<float>::max() != genericDistance)
            return STATUS_CODE_SUCCESS;
    }

    return STATUS_CODE_UNCHANGED;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::GetDistanceToHitInSameLayer(const CaloHit *const pCaloHit, const CaloHitList *const pCaloHitList,
    float &distance) const
{
    const float dCut ((PandoraContentApi::GetGeometry(*this)->GetHitTypeGranularity(pCaloHit->GetHitType()) <= FINE) ?
        (m_sameLayerPadWidthsFine * pCaloHit->GetCellLengthScale()) :
        (m_sameLayerPadWidthsCoarse * pCaloHit->GetCellLengthScale()) );

    if (dCut < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_FAILURE;

    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

    bool hitFound(false);
    float smallestDistanceSquared(std::numeric_limits<float>::max());
    const float rDCutSquared(1.f / (dCut * dCut));

    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        const CaloHit *const pHitInCluster = *iter;
        const CartesianVector &hitInClusterPosition(pHitInCluster->GetPositionVector());
        const float separationSquared((hitPosition - hitInClusterPosition).GetMagnitudeSquared());
        const float hitDistanceSquared(separationSquared * rDCutSquared);

        if (hitDistanceSquared < smallestDistanceSquared)
        {
            smallestDistanceSquared = hitDistanceSquared;
            hitFound = true;
        }
    }

    if (!hitFound)
        return STATUS_CODE_UNCHANGED;

    distance = std::sqrt(smallestDistanceSquared);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::GetConeApproachDistanceToHit(const CaloHit *const pCaloHit, const CaloHitList *const pCaloHitList,
    const CartesianVector &clusterDirection, float &distance) const
{
    bool hitFound(false);
    float smallestDistance(std::numeric_limits<float>::max());

    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        const CaloHit *const pHitInCluster = *iter;
        float hitDistance(std::numeric_limits<float>::max());

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_UNCHANGED, !=, this->GetConeApproachDistanceToHit(pCaloHit,
            pHitInCluster->GetPositionVector(), clusterDirection, hitDistance));

        if (hitDistance < smallestDistance)
        {
            smallestDistance = hitDistance;
            hitFound = true;
        }
    }

    if (!hitFound)
        return STATUS_CODE_UNCHANGED;

    distance = smallestDistance;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::GetConeApproachDistanceToHit(const CaloHit *const pCaloHit, const CartesianVector &clusterPosition,
    const CartesianVector &clusterDirection, float &distance) const
{
    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
    const CartesianVector positionDifference(hitPosition - clusterPosition);

    // just a reminder that clusterDirection is determined by the track seeding the cluster. So it is the track position @ ECal face.
    // std::cout << "hit position is: " << hitPosition << std::endl;
    // std::cout << "track position @ ECal face is: " << clusterPosition << std::endl;
    // std::cout << "track radial distance @ ECal face: " << std::sqrt(clusterPosition.GetX()*clusterPosition.GetX() + clusterPosition.GetY()*clusterPosition.GetY()) << std::endl;
    // std::cout << "position difference between hit and track @ ECal face is: " << positionDifference.GetMagnitudeSquared() << std::endl;

    if (positionDifference.GetMagnitudeSquared() > m_coneApproachMaxSeparation2) { 
        // std::cout << "GP ERROR: Track (@ ECal face) and hit position very inconsistent!" << std::endl;
        // std::cout << "position difference between hit and track @ ECal face was: " << positionDifference.GetMagnitudeSquared() << std::endl;
        return STATUS_CODE_UNCHANGED;
    }
    const float dAlong(clusterDirection.GetDotProduct(positionDifference));

    // if (dAlong > m_maxClusterDirProjection) {
    //     std::cout << "GP MESSAGE: dAlong is too large. Perhaps the track and hit are far away..?" << std::endl << "dAlong is: " << dAlong << std::endl;
    // }
    // if (dAlong < m_minClusterDirProjection) {
    //     std::cout << "GP ERROR: Track and hit are in opposite directions. Something might be wrong." << std::endl; 
    // }

    if ((dAlong < m_maxClusterDirProjection) && (dAlong > m_minClusterDirProjection))
    {
        const float dCut ((PandoraContentApi::GetGeometry(*this)->GetHitTypeGranularity(pCaloHit->GetHitType()) <= FINE) ?
            (std::fabs(dAlong) * m_tanConeAngleFine) + (m_additionalPadWidthsFine * pCaloHit->GetCellLengthScale()) :
            (std::fabs(dAlong) * m_tanConeAngleCoarse) + (m_additionalPadWidthsCoarse * pCaloHit->GetCellLengthScale()) );

        if (dCut < std::numeric_limits<float>::epsilon()) {
            std::cout << "GP ERROR: For some reason, the hit and found track are back-to-back. This shouldn't happen." << std::endl;
            return STATUS_CODE_FAILURE;
        }
        const float dPerp (clusterDirection.GetCrossProduct(positionDifference).GetMagnitude());

        distance = dPerp / dCut;
        // std::cout << "this is the distance calculated: " << distance << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    return STATUS_CODE_UNCHANGED;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::GetDistanceToTrackSeed(const Cluster *const pCluster, const CaloHit *const pCaloHit, unsigned int searchLayer,
    float &distance) const
{
    if (searchLayer < m_maxLayersToTrackSeed)
        return this->GetDistanceToTrackSeed(pCluster, pCaloHit, distance);

    const int searchLayerInt(static_cast<int>(searchLayer));
    const int startLayer(std::max(0, searchLayerInt - static_cast<int>(m_maxLayersToTrackLikeHit)));

    const OrderedCaloHitList &orderedCaloHitList = pCluster->GetOrderedCaloHitList();

    for (int iLayer = startLayer; iLayer < searchLayerInt; ++iLayer)
    {
        OrderedCaloHitList::const_iterator listIter = orderedCaloHitList.find(iLayer);
        if (orderedCaloHitList.end() == listIter)
            continue;

        for (CaloHitList::const_iterator iter = listIter->second->begin(), iterEnd = listIter->second->end(); iter != iterEnd; ++iter)
        {
            float tempDistance(std::numeric_limits<float>::max());
            PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_UNCHANGED, !=, this->GetDistanceToTrackSeed(pCluster, *iter,
                tempDistance));

            if (tempDistance < m_genericDistanceCut)
                return this->GetDistanceToTrackSeed(pCluster, pCaloHit, distance);
        }
    }

    return STATUS_CODE_UNCHANGED;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::GetDistanceToTrackSeed(const Cluster *const pCluster, const CaloHit *const pCaloHit, float &distance) const
{
    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
    const CartesianVector &trackSeedPosition(pCluster->GetTrackSeed()->GetTrackStateAtCalorimeter().GetPosition());

    const CartesianVector positionDifference(hitPosition - trackSeedPosition);
    const float separationSquared(positionDifference.GetMagnitudeSquared());

    if (separationSquared < m_maxTrackSeedSeparation2)
    {
        const float flexibility(1.f + (m_trackPathWidth * std::sqrt(separationSquared / m_maxTrackSeedSeparation2)));

        const float dCut ((PandoraContentApi::GetGeometry(*this)->GetHitTypeGranularity(pCaloHit->GetHitType()) <= FINE) ?
            flexibility * (m_additionalPadWidthsFine * pCaloHit->GetCellLengthScale()) :
            flexibility * (m_additionalPadWidthsCoarse * pCaloHit->GetCellLengthScale()) );

        if (dCut < std::numeric_limits<float>::epsilon())
            return STATUS_CODE_FAILURE;

        const float dPerp((pCluster->GetInitialDirection().GetCrossProduct(positionDifference)).GetMagnitude());

        distance = dPerp / dCut;
        return STATUS_CODE_SUCCESS;
    }

    return STATUS_CODE_UNCHANGED;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::RemoveEmptyClusters(const ClusterVector &clusterVector) const
{
    ClusterList clusterDeletionList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (0 == (*iter)->GetNCaloHits())
        {
            clusterDeletionList.push_back(*iter);
        }
    }

    if (!clusterDeletionList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &clusterDeletionList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Track seeding parameters
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterSeedStrategy", m_clusterSeedStrategy));

    // High level clustering parameters
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldUseOnlyECalHits", m_shouldUseOnlyECalHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldUseIsolatedHits", m_shouldUseIsolatedHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayersToStepBackFine", m_layersToStepBackFine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayersToStepBackCoarse", m_layersToStepBackCoarse));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterFormationStrategy", m_clusterFormationStrategy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GenericDistanceCut", m_genericDistanceCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitTrackCosAngle", m_minHitTrackCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitClusterCosAngle", m_minHitClusterCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldUseTrackSeed", m_shouldUseTrackSeed));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackSeedCutOffLayer", m_trackSeedCutOffLayer));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldFollowInitialDirection", m_shouldFollowInitialDirection));

    // Same layer distance parameters
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SameLayerPadWidthsFine", m_sameLayerPadWidthsFine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SameLayerPadWidthsCoarse", m_sameLayerPadWidthsCoarse));

    // Cone approach distance parameters
    float coneApproachMaxSeparation = std::sqrt(m_coneApproachMaxSeparation2);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ConeApproachMaxSeparation", coneApproachMaxSeparation));
    m_coneApproachMaxSeparation2 = coneApproachMaxSeparation * coneApproachMaxSeparation;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TanConeAngleFine", m_tanConeAngleFine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TanConeAngleCoarse", m_tanConeAngleCoarse));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AdditionalPadWidthsFine", m_additionalPadWidthsFine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AdditionalPadWidthsCoarse", m_additionalPadWidthsCoarse));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterDirProjection", m_maxClusterDirProjection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterDirProjection", m_minClusterDirProjection));

    // Track seed distance parameters
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackPathWidth", m_trackPathWidth));

    float maxTrackSeedSeparation = std::sqrt(m_maxTrackSeedSeparation2);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrackSeedSeparation", maxTrackSeedSeparation));
    m_maxTrackSeedSeparation2 = maxTrackSeedSeparation * maxTrackSeedSeparation;

    if (m_shouldUseTrackSeed && (m_maxTrackSeedSeparation2 < std::numeric_limits<float>::epsilon()))
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLayersToTrackSeed", m_maxLayersToTrackSeed));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLayersToTrackLikeHit", m_maxLayersToTrackLikeHit));

    // Cluster current direction and mip track parameters
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NLayersSpannedForFit", m_nLayersSpannedForFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NLayersSpannedForApproxFit", m_nLayersSpannedForApproxFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NLayersToFit", m_nLayersToFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NLayersToFitLowMipCut", m_nLayersToFitLowMipCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NLayersToFitLowMipMultiplier", m_nLayersToFitLowMipMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FitSuccessDotProductCut1", m_fitSuccessDotProductCut1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FitSuccessChi2Cut1", m_fitSuccessChi2Cut1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FitSuccessDotProductCut2", m_fitSuccessDotProductCut2));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FitSuccessChi2Cut2", m_fitSuccessChi2Cut2));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MipTrackChi2Cut", m_mipTrackChi2Cut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
