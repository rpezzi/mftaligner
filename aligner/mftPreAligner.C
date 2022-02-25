#include <vector>

class preAlignerMFT
{
 public:
  // Pre alignement functions
  void initialize(std::string geometryFileName = "", std::string alignParamFileName = "");                       // loads geometry, tracks, clusters, apply initial alignment if provided
  void computeResiduals(int layerA = -1, int layerB = -1);                                                       // Store residuals in TProfiles? ; Filter tracks if reference layers are provided
  void drawResiduals();                                                                                          // Draw TProfiles for residuals X and Y for each sensor
  std::vector<o2::detectors::AlignParam> computeCorrectionsFromResiduals();                                      // Create vector of alignment parameters from stored TProfiles
  void applyAlignment(std::vector<o2::detectors::AlignParam> alignment);                                         // Apply alignment corrections to loaded geometry
  void UpdateTrackParameters(int layerA = -1, int layerB = -1);                                                  // Recalculate track parameters mimicking the tracker initialization for alignment linear tracks
  void savePreAlignedGeometry(std::vector<o2::detectors::AlignParam> alignment, std::string alignParamFileName); // Saves ccdb-like alignment objects

  // Pre alignment step based on two reference layers
  // Calculate alignment correction to centralize residuals for all tracks defined by the referece layers
  void preAlignStep(int layer1, int layer2)
  {
    UpdateTrackParameters(layer1, layer2); // Ensure all track parameters are consistent with current geometry
    computeResiduals(layer1, layer2);
    auto thisCorrection = computeCorrectionsFromResiduals();
    applyAlignment(thisCorrection);
    savePreAlignedGeometry(thisCorrection, "o2sim-prealignedmft.root");
  }

  void printTracks(int nTracks = 5)
  {
    for (auto& track : mMFTTracks) {
      std::cout << "X = " << track.getX() << " ; Y = " << track.getY() << " ; Z = " << track.getZ() << " ; phi = " << track.getPhi() << " ; tanl = " << track.getTanl() << " ; => outparam: X = " << track.getOutParam().getX() << " ; Y = " << track.getOutParam().getY() << " ; Z = " << track.getOutParam().getZ() << " ; phi = " << track.getOutParam().getPhi() << " ; tanl = " << track.getOutParam().getTanl() << std::endl;
      if (!nTracks--) {
        break;
      }
    }
  }

  o2::math_utils::Point3D<float> toGlobalCoordinates(o2::itsmft::CompClusterExt cluster)
  {
    // Get the global position of the cluster
    auto chipID = cluster.getChipID();
    auto pattID = cluster.getPatternID();

    o2::math_utils::Point3D<float> locC;

    if (pattID != o2::itsmft::CompCluster::InvalidPatternID) {
      locC = mDict.getClusterCoordinates(cluster);
    }

    // Transformation to the local --> global
    auto gloC = gman->getMatrixL2G(chipID) * locC;
    return gloC;
  }

  auto getCluster(o2::mft::TrackMFT mftTrack, int nCluster)
  { // Retrieve the MFT cluster position in global coordinates
    auto offset = mftTrack.getExternalClusterIndexOffset();
    auto clsEntry = trackExtClsVec[offset + nCluster];
    return toGlobalCoordinates(mCompClusters[clsEntry]);
  };

 private:
  std::unique_ptr<TProfile> mXResiduals, mYResiduals;
  std::vector<o2::mft::TrackMFT> mMFTTracks, *trackMFTVecP = &mMFTTracks;
  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  std::vector<o2::itsmft::CompClusterExt> mCompClusters, *clsVecP = &mCompClusters;

  o2::itsmft::TopologyDictionary mDict;
  o2::itsmft::ChipMappingMFT mftChipMapper;
  TTree *mftTrackTree = nullptr, *clsTree = nullptr;
  o2::mft::GeometryTGeo* gman;
};

//_________________________________________________________________________________________________
inline void preAlignerMFT::initialize(std::string geometryFileName = "", std::string alignParamFileName = "")
{
  std::cout << "Initializing MFT pre aligner" << std::endl;

  o2::base::GeometryManager::loadGeometry("", true, true);
  gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));

  // MFT Tracks
  TFile* trkFileIn = new TFile("mfttracks.root");
  mftTrackTree = (TTree*)trkFileIn->Get("o2sim");
  mftTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);
  mftTrackTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);
  std::cout << "Loading MFT tracks file with " << mftTrackTree->GetEntries() << " entries\n";
  mftTrackTree->GetEntry(0);

  // MFT Clusters
  TFile clusterFile("mftclusters.root");
  clsTree = (TTree*)clusterFile.Get("o2sim");
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);
  std::cout << "Loading MFT clusters file with " << clsTree->GetEntries() << " entries\n";
  clsTree->GetEntry(0);

  // Cluster pattern dictionary
  std::string dictfile = "MFTdictionary.bin";

  std::ifstream file(dictfile.c_str());
  if (file.good()) {
    printf("Running with dictionary: %s \n", dictfile.c_str());
    mDict.readBinaryFile(dictfile);
  } else {
    printf("Can not run without dictionary !\n");
    exit;
  }
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::computeResiduals(int layerA = -1, int layerB = -1)
{
  // Compute residuals for each MFT sensor.

  std::cout << "Computing residuals" << std::endl;

  auto layerFilter = [this, layerA, layerB](auto track) {
    auto zLayerA = o2::mft::constants::LayerZPosition[layerA]; // TODO: make this alignment proof
    auto zLayerB = o2::mft::constants::LayerZPosition[layerB];
    return (std::abs(track.getZ() - zLayerA) < 0.2 or std::abs(track.getOutParam().getZ() - zLayerB) < 0.2);
  };

  auto fillResidualProfiles = [this](auto track) {
    // loop over all track clusters, compute and fill TProfile histograms
    for (auto icls = 0; icls < track.getNumberOfPoints(); icls++) {
      const auto cluster = getCluster(track, icls);
      track.propagateParamToZlinear(cluster.z());
      // Calculate residuals
      // Fill histograms
    }
  };

  if (layerA == -1 or layerB == -1) {
    // compute residuals for all tracks and respective clusters
    for (auto& track : mMFTTracks) {
      fillResidualProfiles(track);
    }
  } else {
    // compute unbiased residuals:
    //   1. reference layers are fixed;
    //   2. using only tracks with reference clusters in provided layers
    //   3. This method assumes layerA and layerB are reference for the track sample.
    for (auto& track : mMFTTracks) {
      if (layerFilter(track)) {
        fillResidualProfiles(track);
      }
    }
  }
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::drawResiduals()
{
  std::cout << "Drawing residuals" << std::endl;
}

//_________________________________________________________________________________________________
inline std::vector<o2::detectors::AlignParam> preAlignerMFT::computeCorrectionsFromResiduals()
{
  // Compute alignment parameters to centralize residuals in each MFT sensors
  std::cout << "Computing alignement corrections from residuals" << std::endl;
  std::vector<o2::detectors::AlignParam> params;
  return params;
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::applyAlignment(std::vector<o2::detectors::AlignParam> alignment)
{
  std::cout << "Applying alignment" << std::endl;
  o2::base::GeometryManager::loadGeometry("", true, true);
  gman = o2::mft::GeometryTGeo::Instance();

  for (auto& alpar : alignment) {
    alpar.applyToGeometry();
  }
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::UpdateTrackParameters(int layerA = -1, int layerB = -1)
{
  // This method updates track linear track parameters.
  // By default all tracks are updated using a linear model computed using first and last track-clusters.
  // If reference layers are provided, only those tracks with clusters in both layers are updated

  std::cout << "Updating track parameters for new geometry" << std::endl;

  auto getClusterIdLayer = [this](o2::mft::TrackMFT mftTrack, int layer) { // Returns the MFT
    auto offset = mftTrack.getExternalClusterIndexOffset();
    for (auto icls = 0; icls < mftTrack.getNumberOfPoints(); icls++) {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto clsLayer = mftChipMapper.chip2Layer(mCompClusters[clsEntry].getChipID());
      if (layer == clsLayer) {
        return icls;
      }
    }
    return -1;
  };

  auto updateTrackParam = [](o2::mft::TrackMFT& track, const auto clusterA, const auto clusterB) { // Set linear track parameters
    auto xA = clusterA.X();
    auto yA = clusterA.Y();
    auto zA = clusterA.Z();

    auto xB = clusterB.X();
    auto yB = clusterB.Y();
    auto zB = clusterB.Z();

    auto deltaX = xA - xB;
    auto deltaY = yA - yB;
    auto deltaZ = zA - zB;
    auto deltaR = TMath::Sqrt(deltaX * deltaX + deltaY * deltaY);
    auto tanl0 = -std::abs(deltaZ) / deltaR;
    double phi0 = TMath::ATan2(-deltaY, -deltaX);

    track.setX(xA);
    track.setY(yA);
    track.setZ(zA);
    track.setPhi(phi0);
    track.setTanl(tanl0);

    o2::track::TrackParCovFwd outParam;
    outParam.setX(xB);
    outParam.setY(yB);
    outParam.setZ(zB);
    outParam.setPhi(phi0);
    outParam.setTanl(tanl0);
    track.setOutParam(outParam);
  };

  if (layerA == -1 or layerB == -1) { // update all tracks using first and last clusters as reference
    for (auto& track : mMFTTracks) {
      auto clusterA = getCluster(track, 0);
      auto clusterB = getCluster(track, track.getNumberOfPoints() - 1);
      updateTrackParam(track, clusterA, clusterB);
    }
  } else { // Update tracks using layers A and B as reference
    for (auto& track : mMFTTracks) {
      auto clusterIdA = getClusterIdLayer(track, layerA);
      auto clusterIdB = getClusterIdLayer(track, layerB);
      if (clusterIdA != -1 and clusterIdB != -1) { // Ignore tracks with no clusters in reference layers
        auto clusterA = getCluster(track, clusterIdA);
        auto clusterB = getCluster(track, clusterIdB);
        updateTrackParam(track, clusterA, clusterB);
      }
    }
  }
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::savePreAlignedGeometry(std::vector<o2::detectors::AlignParam> alignment, std::string alignParamFileName)
{
  std::cout << "Saving something..." << std::endl;
}

//_________________________________________________________________________________________________
void mftPreAligner()
{
  preAlignerMFT pa;

  pa.initialize();
  pa.printTracks();
  int layerA = 0, layerB = 9;
  pa.preAlignStep(layerA, layerB); // Start with longest tracks
  pa.printTracks();

  pa.computeResiduals(layerA, layerB); // Filtering tracks -> unbiased  residuals should center at zero
  pa.drawResiduals();
  pa.computeResiduals(); // No filtering
  pa.drawResiduals();

  //  preAlignStep(1, 8); // Move to other layers once we are happy with previous step
  //  preAlignStep(2, 7);
  //  preAlignStep(3, 6);
}