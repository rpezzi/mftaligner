#include "mftPreAligner.h"

//_________________________________________________________________________________________________
void mftPreAligner()
{
  preAlignerMFT pa;
  pa.initialize();

  auto debugAlignStep = [&pa](int layerA, int layerB) { // Draws detailed residual status at each alignment step
    pa.computeResiduals();                              // Complete residuals
    pa.drawResiduals("");
    pa.preAlignStep(layerA, layerB); // Unbiased residuals
    pa.drawResiduals("");
    pa.computeResiduals(layerA, layerB); // Corrected residuals / centered at zero
    pa.drawResiduals("");
  };

  // Pre alignment steps
  debugAlignStep(0, 9); // Start with longest tracks
  debugAlignStep(1, 8);
  pa.preAlignStep(2, 7);
  pa.preAlignStep(3, 6);
  pa.preAlignStep(0, 9);
  pa.preAlignStep(1, 8);
  debugAlignStep(0, 9);

  // Compute final residuals and save
  pa.computeResiduals(); // Complete residuals
  pa.drawResiduals("");
  pa.saveResidualCanvases();

  // =====================================================================================================================
  // Sanity test checks
  auto test = !true;
  int layerA = 0;
  int layerB = 9;
  if (test) { // Check validity of file produced in previous step
    std::cout << " Test1: Check use of file produced in previous step\n";
    pa.initialize("mftpre");             // Reload, starting from geometry file produced by previous call to preAlignStep
    pa.computeResiduals(layerA, layerB); // Initial residuals should be identical as previous execution
    pa.drawResiduals("");

    layerA = 1, layerB = 8;
    std::cout << " Test1 contd.: Prealign with reference layers: " << layerA << " and " << layerB << std::endl;
    pa.preAlignStep(layerA, layerB, "mftpre"); // Prealign with different reference layers
    pa.drawResiduals("");
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
  // =====================================================================================================================
}