#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

typedef SingleObjectSelector<
  std::vector<pat::PackedCandidate>, 
  StringCutObjectSelector<pat::PackedCandidate> > PackedCandidateSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PackedCandidateSelector);
