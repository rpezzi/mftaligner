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
    UpdateTrackParameters(); // Ensure all track parameters are consistent with current geometry
    computeResiduals(layer1, layer2);
    auto thisCorrection = computeCorrectionsFromResiduals();
    applyAlignment(thisCorrection);
    savePreAlignedGeometry(thisCorrection, "o2sim-prealignedmft.root");
  }

  void runPrealigner()
  {
    initialize();
    int layerA = 0, layerB = 9;
    preAlignStep(layerA, layerB);     // Start with longest tracks
    computeResiduals(layerA, layerB); // Filtering tracks -> residuals should center at zero
    drawResiduals();
    computeResiduals(); // No filtering
    drawResiduals();

    //  preAlignStep(1, 8); // Move to other layers once we are happy with previous step
    //  preAlignStep(2, 7);
    //  preAlignStep(3, 6);
  }

 private:
  std::unique_ptr<TProfile> mXResiduals, mYResiduals;
  std::vector<o2::itsmft::CompClusterExt> mCompClusters; // Used MFT compact clusters
  std::vector<o2::mft::TrackMFTExt> mTracks;             // MFT Tracks
  o2::itsmft::TopologyDictionary mDict;
};

//_________________________________________________________________________________________________
inline void preAlignerMFT::initialize(std::string geometryFileName = "", std::string alignParamFileName = "")
{
  std::cout << "Initializing MFT pre aligner" << std::endl;
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::computeResiduals(int layerA = -1, int layerB = -1)
{
  std::cout << "Computing residuals" << std::endl;
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::drawResiduals()
{
  std::cout << "Drawing residuals" << std::endl;
}

//_________________________________________________________________________________________________
inline std::vector<o2::detectors::AlignParam> preAlignerMFT::computeCorrectionsFromResiduals()
{
  std::cout << "Comuting alignement corrections from residuals" << std::endl;

  std::vector<o2::detectors::AlignParam> params;
  return params;
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::applyAlignment(std::vector<o2::detectors::AlignParam> alignment)
{
  std::cout << "Applying alignment" << std::endl;
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::UpdateTrackParameters(int layerA = -1, int layerB = -1)
{
  std::cout << "Updating track parameters for new geometry" << std::endl;
}

//_________________________________________________________________________________________________
inline void preAlignerMFT::savePreAlignedGeometry(std::vector<o2::detectors::AlignParam> alignment, std::string alignParamFileName)
{
  std::cout << "Saving something..." << std::endl;
}

void mftPreAligner()
{
  preAlignerMFT prealigner;
  prealigner.runPrealigner();
}