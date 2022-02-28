class preAlignerMFT
{
 public:
  // Pre alignement functions
  void initialize(std::string geometryFileName = "", std::string alignParamFileName = "");                         // loads geometry, tracks, clusters, apply initial alignment if provided
  void computeResiduals(int layerA = -1, int layerB = -1);                                                         // Store residuals in TProfiles? ; Filter tracks if reference layers are provided
  void drawResiduals(std::string compareFile = "");                                                                // Draw TProfiles for residuals X and Y for each sensor
  std::vector<o2::detectors::AlignParam> computeCorrectionsFromResiduals();                                        // Create vector of alignment parameters from stored TProfiles
  void applyAlignment(std::vector<o2::detectors::AlignParam> alignment, std::string geomFile = "");                // Apply alignment corrections to loaded geometry
  void UpdateTrackParameters(int layerA = -1, int layerB = -1);                                                    // Recalculate track parameters mimicking the tracker initialization for alignment linear tracks
  void savePreAlignedParameters(std::vector<o2::detectors::AlignParam> alignment, std::string alignParamFileName); // Saves ccdb-like alignment objects
  void exportPreAlignedGeometry(std::string geomFile = "");                                                        // Export aligned geometry file
  void pushAlignmentToCCDB(std::vector<o2::detectors::AlignParam> alignment, const std::string& ccdbHost,
                           long tmin, long tmax, const std::string& objectPath); // Push alignment parameters to CCDB

  auto& getFinalAlignment() const { return mResultingAlignment; }
  // Pre alignment step based on two reference layers
  // Calculate alignment correction to centralize residuals for all tracks defined by the referece layers
  void preAlignStep(int layer1, int layer2, std::string geomFile = "")
  {
    std::cout << "Prealign step with reference layers: " << layer1 << " and " << layer2 << std::endl;
    UpdateTrackParameters(layer1, layer2); // Ensure all track parameters are consistent with current geometry
    computeResiduals(layer1, layer2);
    drawResiduals("");
    auto thisCorrection = computeCorrectionsFromResiduals();
    applyAlignment(thisCorrection, geomFile);
    exportPreAlignedGeometry("mftpre_geometry-aligned.root");
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

  o2::math_utils::Point3D<float> toGlobalCoordinates(const o2::BaseCluster<float> cluster)
  {
    // Transformation to the local --> global
    auto gloC = gman->getMatrixL2G(cluster.getSensorID()) * cluster.getXYZ();
    return gloC;
  }

  auto getCluster(o2::mft::TrackMFT mftTrack, int nCluster)
  { // Retrieve the MFT cluster position in global coordinates
    auto offset = mftTrack.getExternalClusterIndexOffset();
    auto clsEntry = trackExtClsVec[offset + nCluster];
    return toGlobalCoordinates(mUnCompClusters[clsEntry]);
  };

  auto getTrackClusterChipID(o2::mft::TrackMFT mftTrack, int nCluster)
  { // Retrieve chipID of cluster in track
    auto offset = mftTrack.getExternalClusterIndexOffset();
    auto clsEntry = trackExtClsVec[offset + nCluster];
    return mUnCompClusters[clsEntry].getSensorID();
  };

 private:
  TProfile* mXResiduals = nullptr;
  TProfile* mYResiduals = nullptr;
  std::vector<o2::mft::TrackMFT> mMFTTracks, *trackMFTVecP = &mMFTTracks;
  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  std::vector<o2::itsmft::CompClusterExt> mCompClusters, *clsVecP = &mCompClusters;
  std::vector<unsigned char> mPatterns, *mPatternsP = &mPatterns;

  std::vector<o2::BaseCluster<float>> mUnCompClusters; // MFT Clusters in local coordinate system

  std::vector<o2::detectors::AlignParam> mResultingAlignment;

  o2::itsmft::TopologyDictionary mDict;
  o2::itsmft::ChipMappingMFT mftChipMapper;
  TTree *mftTrackTree = nullptr, *clsTree = nullptr;
  o2::mft::GeometryTGeo* gman;
};

//_________________________________________________________________________________________________
inline void preAlignerMFT::initialize(std::string geometryFileName = "", std::string alignParamFileName = "")
{
  std::cout << "Initializing MFT pre aligner" << std::endl;

  o2::base::GeometryManager::loadGeometry(geometryFileName.c_str(), false, true);
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
  clsTree->SetBranchAddress("MFTClusterPatt", &mPatternsP);
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

  // Cache compact cluster positions (use topology dictionary and clusters patterns only once)
  if (!mUnCompClusters.size()) {
    mUnCompClusters.reserve(mCompClusters.size());
    std::vector<unsigned char>::iterator pattIt = mPatterns.begin();

    for (const auto& compCluster : mCompClusters) {

      auto chipID = compCluster.getChipID();
      auto pattID = compCluster.getPatternID();

      o2::math_utils::Point3D<float> locC;
      float sigmaX2 = 1.806335945e-06F, sigmaY2 = 2.137443971e-06F; // Dummy COG errors (about half pixel size) from o2::mft::ioutils::DefClusError2Row and o2::mft::ioutils::DefClusError2Col
      if (pattID != o2::itsmft::CompCluster::InvalidPatternID) {
        sigmaX2 = mDict.getErr2X(pattID); // ALPIDE local X coordinate => MFT global X coordinate (ALPIDE rows)
        sigmaY2 = mDict.getErr2Z(pattID); // ALPIDE local Z coordinate => MFT global Y coordinate (ALPIDE columns)
        if (!mDict.isGroup(pattID)) {
          locC = mDict.getClusterCoordinates(compCluster);
        } else {
          o2::itsmft::ClusterPattern patt(pattIt);
          locC = mDict.getClusterCoordinates(compCluster, patt);
        }
      } else {
        o2::itsmft::ClusterPattern patt(pattIt);
        locC = mDict.getClusterCoordinates(compCluster, patt, false);
      }
      mUnCompClusters.emplace_back(chipID, locC);
    }
  }

  // Initialize total alignment vector
  if (!mResultingAlignment.size()) {
    Int_t nChip = 0;
    mResultingAlignment.clear();
    TString sname = gman->composeSymNameMFT();
    Int_t nHalf = gman->getNumberOfHalfs();
    mResultingAlignment.emplace_back(sname, -1, 0, 0, 0, 0, 0, 0, true);

    for (Int_t hf = 0; hf < nHalf; hf++) {
      Int_t nDisks = gman->getNumberOfDisksPerHalf(hf);
      sname = gman->composeSymNameHalf(hf);
      mResultingAlignment.emplace_back(sname, -1, 0, 0, 0, 0, 0, 0, true);

      for (Int_t dk = 0; dk < nDisks; dk++) {
        sname = gman->composeSymNameDisk(hf, dk);
        mResultingAlignment.emplace_back(sname, -1, 0, 0, 0, 0, 0, 0, true);

        Int_t nLadders = 0;
        for (Int_t sensor = gman->getMinSensorsPerLadder();
             sensor < gman->getMaxSensorsPerLadder() + 1; sensor++) {
          nLadders += gman->getNumberOfLaddersPerDisk(hf, dk, sensor);
        }

        for (Int_t lr = 0; lr < nLadders; lr++) { // nLadders
          sname = gman->composeSymNameLadder(hf, dk, lr);
          Int_t nSensorsPerLadder = gman->getNumberOfSensorsPerLadder(hf, dk, lr);
          mResultingAlignment.emplace_back(sname, -1, 0, 0, 0, 0, 0, 0, true);

          for (Int_t sr = 0; sr < nSensorsPerLadder; sr++) {
            sname = gman->composeSymNameChip(hf, dk, lr, sr);
            // Follow the geometrical order for sensors construction
            Int_t chipID = o2::itsmft::ChipMappingMFT::mChipIDGeoToRO[nChip++];
            Int_t uid = o2::base::GeometryManager::getSensID(o2::detectors::DetID::MFT, chipID);
            mResultingAlignment.emplace_back(sname, uid, 0, 0, 0, 0, 0, 0, true);
            nChip++;
          }
        }
      }
    }
  }
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::computeResiduals(int layerA = -1, int layerB = -1)
{
  // Compute residuals for each MFT sensor.
  static int step = 0;

  std::cout << "Computing residuals" << std::endl;

  if (mXResiduals) {
    delete mXResiduals;
  }

  if (mYResiduals) {
    delete mYResiduals;
  }

  mXResiduals = new TProfile(Form("mXResiduals_%d", step), "x residuals vs chipID", 936, 0., 936);
  mXResiduals->SetXTitle("cluster.chipID ");
  mXResiduals->SetYTitle("#delta_{x} (cm)");
  mXResiduals->Sumw2();

  mYResiduals = new TProfile(Form("mYResiduals_%d", step), "y residuals vs chipID", 936, 0., 936);
  mYResiduals->SetXTitle("cluster.chipID ");
  mYResiduals->SetYTitle("#delta_{y} (cm)");
  mYResiduals->Sumw2();
  step++;

  auto layerFilter = [this, layerA, layerB](auto track) {
    auto zLayerA = o2::mft::constants::LayerZPosition[layerA]; // TODO: make this alignment proof
    auto zLayerB = o2::mft::constants::LayerZPosition[layerB];
    return (std::abs(track.getZ() - zLayerA) < 0.2 and std::abs(track.getOutParam().getZ() - zLayerB) < 0.2);
  };

  auto fillResidualProfiles = [this](auto track) {
    // loop over all track clusters, compute and fill TProfile histograms
    for (auto icls = 0; icls < track.getNumberOfPoints(); icls++) {
      const auto cluster = getCluster(track, icls);
      const auto chipID = getTrackClusterChipID(track, icls);
      track.propagateParamToZlinear(cluster.z());

      // Calculate residuals
      auto resX = track.getX() - cluster.x();
      auto resY = track.getY() - cluster.y();

      // Fill histograms
      if (resX * resY) {
        std::cout << " Residuals for icls " << icls << ": resX = " << resX << " ; resY = " << resY << " from chipID = " << chipID << std::endl;
        mXResiduals->Fill(1.0 * chipID, 1.0 * resX);
        mYResiduals->Fill(chipID, resY);
      }
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
inline void preAlignerMFT::drawResiduals(std::string compareFile = "")
{
  std::cout << "Drawing residuals" << std::endl;
  static int step = 0;
  auto c = new TCanvas(Form("c_%d", step), Form("c_%d", step), 1200, 800);
  step++;
  c->Divide(2, 1);

  c->cd(1);
  mXResiduals->DrawClone("EPZ");
  c->cd(2);
  mYResiduals->DrawClone("EPZ");

  /*if (!compareFile.empty()) {
    TFile* f = new TFile(compareFile.c_str(), "read");
    auto leg = new TLegend();
    TProfile* mXResidualsToCompare;
    TProfile* mYResidualsToCompare;
    mXResidualsToCompare = (TProfile*)f->Get("mXResiduals");
    mYResidualsToCompare = (TProfile*)f->Get("mYResiduals");

    c->cd(1);
    mXResidualsToCompare->Draw("EPZ same");
    mXResidualsToCompare->SetLineColor(kRed + 1);
    c->cd(2);
    mYResidualsToCompare->Draw("EPZ same");
    mYResidualsToCompare->SetLineColor(kRed + 1);

    leg->AddEntry(*mYResiduals, "prealign", "pl");
    leg->AddEntry(*mYResidualsToCompare, "file to compare", "pl");
    leg->SetBorderSize(0);
    leg->Draw();
  } */
}

//_________________________________________________________________________________________________
inline std::vector<o2::detectors::AlignParam> preAlignerMFT::computeCorrectionsFromResiduals()
{
  // Compute alignment parameters to centralize residuals in each MFT sensors
  std::cout << "Computing alignement corrections from residuals" << std::endl;
  std::vector<o2::detectors::AlignParam> params;

  double lPsi = 0, lTheta = 0, lPhi = 0.;
  Int_t nChip = 0;
  bool glo = true;

  TString sname = gman->composeSymNameMFT();
  Int_t nHalf = gman->getNumberOfHalfs();

  // Fill empty align param for each element which is not a sensor
  params.emplace_back(sname, -1, 0, 0, 0, 0, 0, 0, glo);

  for (Int_t hf = 0; hf < nHalf; hf++) {
    Int_t nDisks = gman->getNumberOfDisksPerHalf(hf);
    sname = gman->composeSymNameHalf(hf);
    params.emplace_back(sname, -1, 0, 0, 0, 0, 0, 0, glo);

    for (Int_t dk = 0; dk < nDisks; dk++) {
      sname = gman->composeSymNameDisk(hf, dk);
      params.emplace_back(sname, -1, 0, 0, 0, 0, 0, 0, glo);

      Int_t nLadders = 0;
      for (Int_t sensor = gman->getMinSensorsPerLadder();
           sensor < gman->getMaxSensorsPerLadder() + 1; sensor++) {
        nLadders += gman->getNumberOfLaddersPerDisk(hf, dk, sensor);
      }

      for (Int_t lr = 0; lr < nLadders; lr++) { // nLadders
        sname = gman->composeSymNameLadder(hf, dk, lr);
        Int_t nSensorsPerLadder = gman->getNumberOfSensorsPerLadder(hf, dk, lr);
        params.emplace_back(sname, -1, 0, 0, 0, 0, 0, 0, glo);

        for (Int_t sr = 0; sr < nSensorsPerLadder; sr++) {
          sname = gman->composeSymNameChip(hf, dk, lr, sr);
          // Follow the geometrical order for sensors construction
          Int_t chipID = o2::itsmft::ChipMappingMFT::mChipIDGeoToRO[nChip++];
          Int_t uid = o2::base::GeometryManager::getSensID(o2::detectors::DetID::MFT, chipID);
          auto resX = mXResiduals->GetBinContent(chipID + 1);
          auto resY = mYResiduals->GetBinContent(chipID + 1);
          if (resX * resY)
            printf("resX=%f, resY=%f for chipID %d \n ", resX, resY, chipID);
          params.emplace_back(sname, uid, resX, resY, 0, lPsi, lTheta, lPhi, glo);
        }
      }
    }
  }

  return params;
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::applyAlignment(std::vector<o2::detectors::AlignParam> alignment, std::string geomFile)
{
  std::cout << "Applying alignment" << std::endl;
  gman = o2::mft::GeometryTGeo::Instance();

  o2::base::GeometryManager::loadGeometry(geomFile.c_str(), false, true);
  o2::base::GeometryManager::applyAlignment(alignment);

  int index = 0;
  for (auto& alpar : alignment) {
    // alpar.applyToGeometry();
    // alpar.print();
    mResultingAlignment[index].setGlobalParams(mResultingAlignment[index].createMatrix() * alpar.createMatrix()); // Increment aligment params: is this true?
    index++;
  }
  gman->updateL2GMatrixCache();
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::UpdateTrackParameters(int layerA = -1, int layerB = -1)
{
  // This method updates track linear track parameters.
  // By default all tracks are updated using a linear model computed using first and last track-clusters.
  // If reference layers are provided, only those tracks with clusters in both layers are updated

  std::cout << "Updating track parameters for current geometry" << std::endl;

  auto getClusterIdLayer = [this](o2::mft::TrackMFT mftTrack, int layer) { // Returns the MFT
    auto offset = mftTrack.getExternalClusterIndexOffset();
    for (auto icls = 0; icls < mftTrack.getNumberOfPoints(); icls++) {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto clsLayer = mftChipMapper.chip2Layer(mUnCompClusters[clsEntry].getSensorID());
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
inline void preAlignerMFT::savePreAlignedParameters(std::vector<o2::detectors::AlignParam> alignment, std::string alignParamFileName)
{
  std::cout << "Saving alignment parameters to " << alignParamFileName << std::endl;

  if (!alignParamFileName.empty()) {
    printf("Storing MFT alignment in local file %s \n", alignParamFileName.c_str());
    TFile algFile(alignParamFileName.c_str(), "recreate");
    algFile.WriteObjectAny(&alignment, "std::vector<o2::detectors::AlignParam>", "alignment");
    algFile.Close();
  }
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::exportPreAlignedGeometry(std::string geomFile)
{
  std::cout << "Exporting geometry file to " << geomFile << std::endl;

  gGeoManager->Export(geomFile.c_str());
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::pushAlignmentToCCDB(std::vector<o2::detectors::AlignParam> alignment, const std::string& ccdbHost, long tmin, long tmax, const std::string& objectPath)
{
  /*
  std::cout << "Pushing alignment to CCDB..." << std::endl;
  if (!ccdbHost.empty() {
    std::string path = objectPath.empty() ? o2::base::DetectorNameConf::getAlignmentPath(o2::detectors::DetID::MFT) : objectPath;
    path = "MFT/test_CCDB/MFT"; // testing the ccdb
    printf("Storing alignment object on  %s/%s \n", ccdbHost, path);
    o2::ccdb::CcdbApi api;
    std::map<std::string, std::string> metadata; // can be empty
    api.init(ccdbHost.c_str());                  // or http://localhost:8080 for a local installation
    metadata["test"] = fmt::format("Alignment objects for DetID:{}", o2::detectors::DetID::MFT);
    // store abitrary user object in strongly typed manner
    api.storeAsTFileAny(&alignment, path, metadata, tmin, tmax);
   }
   */
}

//_________________________________________________________________________________________________
void mftPreAligner()
{
  preAlignerMFT pa;
  pa.initialize();
  pa.printTracks();
  int layerA = 0, layerB = 9;
  pa.preAlignStep(layerA, layerB);          // Start with longest tracks
  pa.UpdateTrackParameters(layerA, layerB); // Ensure all track parameters are consistent with current geometry
  pa.computeResiduals(layerA, layerB);
  pa.drawResiduals("");

  auto test = true;
  if (test) { // Check validity of file produced in previous step
    std::cout << " Test1: Check use of file produced in previous step\n";
    pa.initialize("mftpre"); // Reload, starting from geometry file produced by previous call to preAlignStep
    pa.UpdateTrackParameters(layerA, layerB);
    pa.computeResiduals(layerA, layerB); // Initial residuals should be identical as previous execution
    pa.drawResiduals("");

    layerA = 1, layerB = 8;
    std::cout << " Test1 contd.: Prealign with reference layers: " << layerA << " and " << layerB << std::endl;
    pa.preAlignStep(layerA, layerB, "mftpre"); // Prealign with different reference layers
    pa.UpdateTrackParameters(layerA, layerB);
    pa.computeResiduals(layerA, layerB);
    pa.drawResiduals("");
  }

  if (test) { // Are accumulated alignment parameters are equivalent to previous step?
    std::cout << " Test2: Test if accumulated alignment parameters are equivalent to previous step\n";
    int layerA = 1, layerB = 8;
    pa.initialize();                                     // Reload ideal geometry
    auto alignPar = pa.getFinalAlignment();              // Get accumulated deltas from all previous pre-alignment steps
    o2::base::GeometryManager::applyAlignment(alignPar); // Apply and check resisuals
    pa.UpdateTrackParameters(layerA, layerB);
    pa.computeResiduals(layerA, layerB);
    pa.drawResiduals("");
  }

  //  preAlignStep(1, 8); // Move to other layers once we are happy with previous step
  //  preAlignStep(2, 7);
  //  preAlignStep(3, 6);
}