#pragma once

#include "BambooFactory.hh"
#include "BambooAnalysis.hh"

class RelicsOpticalAnalysis : public BambooAnalysis {
  public:
    RelicsOpticalAnalysis(const BambooParameters &pars);

    ~RelicsOpticalAnalysis() = default;

    static AnalysisRegister<RelicsOpticalAnalysis> reg;
};
