// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/// \author Peter Hristov <Peter.Hristov@cern.ch>, CERN
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Henrique J C Zanoli <henrique.zanoli@cern.ch>, Utrecht University

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Analysis/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Analysis/MC.h"

#include <TH1F.h>
#include <TH2F.h>
#include <cmath>

namespace o2fw = o2::framework;
namespace o2exp = o2::framework::expressions;
namespace o2df = o2::dataformats;

namespace track_utils
{
/// Converts a angle phi to the -pi to pi range.
double ConvertPhiRange(double phi)
{
  if (phi > M_PI) {
    phi -= 2 * M_PI;
  } else if (phi < -M_PI) {
    phi += 2 * M_PI;
  }

  return phi;
}

/// Determines the impact parameter and its error for a given track.
/// \param track the track to get the impact parameter from.
/// \param primaryVertex the primary vertex of th collision.
/// \param impactParameterRPhi variable to save the impact parameter in micrometers.
/// \param impactParameterRPhiError variable to save the impact parameter error in micrometers.
template <typename Track>
bool GetImpactParameterAndError(const Track& track, const o2df::VertexBase& primaryVertex,
                                double& impactParameterRPhi, double& impactParameterRPhiError,
                                double& impactParameterZ, double& impactParameterErrorZ)
{
  impactParameterRPhi = -999.;
  impactParameterRPhiError = -999.;
  impactParameterZ = -999;
  impactParameterErrorZ = -999;

  o2df::DCA dca;
  // FIXME: get this from CCDB
  constexpr float magneticField{5.0}; // in kG
  auto trackParameter = getTrackParCov(track);
  bool propagate = trackParameter.propagateToDCA(primaryVertex, magneticField, &dca);

  constexpr float conversion_to_micrometer = 1000;
  if (propagate) {
    impactParameterRPhi = conversion_to_micrometer * dca.getY();
    impactParameterRPhiError = conversion_to_micrometer * std::sqrt(dca.getSigmaY2());
    impactParameterZ = conversion_to_micrometer * dca.getZ();
    impactParameterErrorZ = conversion_to_micrometer * std::sqrt(dca.getSigmaZ2());
  }
  return propagate;
}
} // namespace track_utils

namespace o2::qa
{
/// Class to abstract the naming of a particular feature. It can help you to build the histogram
/// labels in a consistent way, and can generate the titles. You can also pass the
/// plotWithMatplotlib when constructing it to make the labels friendly to matplotlib.
class Feature
{
 public:
  Feature() = default;
  Feature(std::string name, std::string unit, bool plotWithMatplotlib = false)
    : mPlotWithMatplotlib(plotWithMatplotlib), mName(std::move(name)), mUnit(std::move(unit)){};

  std::string NameNoUnit() const { return FormatString(mName); };
  std::string Name() const { return FormatString(mName) + UnitFormatted(); };
  std::string MC() const { return FormatString(MCRaw()) + UnitFormatted(); };
  std::string Rec() const { return FormatString(RecRaw()) + UnitFormatted(); };
  std::string MCRecDiff() const
  {
    return FormatString(MCRaw() + " - " + RecRaw()) + UnitFormatted();
  };

  std::string RelativeMCRecDiff() const
  {
    return FormatString("(" + MCRaw() + " - " + RecRaw() + ")/(" + RecRaw() + ")");
  };
  std::string Unit() const { return "[" + UnitRaw() + "]"; };

 private:
  std::string MCRaw() const { return mName + "^{MC}"; };
  std::string RecRaw() const { return mName + "^{Rec}"; };
  std::string UnitRaw() const { return mUnit; };
  std::string UnitFormatted() const
  {
    if (UnitRaw().empty()) {
      return "";
    }

    return " [" + UnitRaw() + "]";
  };

  std::string FormatString(const std::string& s) const
  {
    if (mPlotWithMatplotlib)
      return "$" + s + "$";
    return s;
  }

  const bool mPlotWithMatplotlib{false};
  const std::string mName;
  const std::string mUnit;
};

/// Handle consistent naming of the histograms in the task and adds the possibility to make the
/// histogram titles friendly to matplotlib.
class Features
{
 public:
  Features(bool plotWithMatplotlib = false)
    : mMathSymbol(plotWithMatplotlib ? "\\" : "#"),
      mTrackMultiplicity{"Track Multiplicity", "", plotWithMatplotlib},
      mEta{Math("eta"), "", plotWithMatplotlib},
      mPhi{Math("varphi"), "rad", plotWithMatplotlib},
      mPt{"p_{T}", "GeV/c", plotWithMatplotlib},
      mVertexX{"X", "cm", plotWithMatplotlib},
      mVertexY{"Y", "cm", plotWithMatplotlib},
      mVertexZ{"Z", "cm", plotWithMatplotlib},
      mImpactParameterRPhi{"Impact Parameter r" + Math("varphi"), Math("mu") + "m",
                           plotWithMatplotlib},
      mImpactParameterRPhiError{"Impact Parameter Error r" + Math("varphi"), Math("mu") + "m",
                                plotWithMatplotlib},
      mImpactParameterZ{"Impact Parameter Z", Math("mu") + "m", plotWithMatplotlib},
      mImpactParameterZError{"Impact Parameter Z Error ", Math("mu") + "m", plotWithMatplotlib} {};

  const Feature& Eta() const { return mEta; }
  const Feature& Phi() const { return mPhi; }
  const Feature& Pt() const { return mPt; }
  const Feature& VertexX() const { return mVertexX; }
  const Feature& VertexY() const { return mVertexY; }
  const Feature& VertexZ() const { return mVertexZ; }
  const Feature& TrackMultiplicity() const { return mTrackMultiplicity; }
  const Feature& ImpactParameterRPhi() const { return mImpactParameterRPhi; }
  const Feature& ImpactParameterRPhiError() const { return mImpactParameterRPhiError; }
  const Feature& ImpactParameterZ() const { return mImpactParameterZ; }
  const Feature& ImpactParameterZError() const { return mImpactParameterZError; }
  const std::string& Counts() const { return mCounts; }

  /// Given a string with the content of the X title, builds the title field for the histogram using
  /// an empty title and the string that represents the number of counts.
  std::string Title1D(const std::string& xAxis) const { return ";" + xAxis + ";" + Counts(); }

  /// Given a string with the content of the X and Y titles, builds the title field for the
  /// histogram using an empty title and the string that represents the number of counts.
  std::string Title2D(const std::string& xAxis, const std::string& yAxis) const
  {
    return ";" + xAxis + ";" + yAxis + ";" + Counts();
  }

 private:
  /// Given a text (for example "varphi") returns the mathematical representation, using the current
  /// value defined in mMathSymbol (for example "#varphi" for ROOT or "\varphi" for matplotlib).
  std::string Math(const std::string& command) { return mMathSymbol + command; }
  const std::string mMathSymbol;
  const Feature mTrackMultiplicity;
  const Feature mEta;
  const Feature mPhi;
  const Feature mPt;
  const Feature mVertexX;
  const Feature mVertexY;
  const Feature mVertexZ;
  const Feature mImpactParameterRPhi;
  const Feature mImpactParameterRPhiError;
  const Feature mImpactParameterZ;
  const Feature mImpactParameterZError;

  const std::string mCounts{"Counts"};
};
} // namespace o2::qa

/// Task to QA global observables of the event
struct QAGlobalObservables {
  o2fw::Configurable<bool> histogramAxisForMatplotlib{
    "histogramAxisForMatplotlib", false,
    "Sets the histograms title to be friendly to matplotlib instead of ROOT"};
  o2::qa::Features features = o2::qa::Features(histogramAxisForMatplotlib);

  std::array<float, 2> collisionPositionRange = {-20., 20.};
  std::array<float, 2> numberOfTracksRange = {0, 2000};

  o2fw::Configurable<int> nBinsNumberOfTracks{"nBinsNumberOfTracks", 2000,
                                              "Number of bins fot the Number of Tracks"};
  o2fw::Configurable<int> nBinsVertexPosition{"nBinsPt", 100,
                                              "Number of bins for the Vertex Position"};

  o2fw::OutputObj<TH1F> hCollisionX{
    TH1F("collisionX", features.Title1D(features.VertexX().Name()).c_str(), nBinsVertexPosition,
         collisionPositionRange[0], collisionPositionRange[1])};

  o2fw::OutputObj<TH1F> hCollisionY{
    TH1F("collisionY", features.Title1D(features.VertexY().Name()).c_str(), nBinsVertexPosition,
         collisionPositionRange[0], collisionPositionRange[1])};

  o2fw::OutputObj<TH1F> hCollisionZ{
    TH1F("collisionZ", features.Title1D(features.VertexZ().Name()).c_str(), nBinsVertexPosition,
         collisionPositionRange[0], collisionPositionRange[1])};

  o2fw::OutputObj<TH1F> hNumberOfTracks{
    TH1F("NumberOfTracks", features.Title1D(features.TrackMultiplicity().Name()).c_str(),
         nBinsNumberOfTracks, numberOfTracksRange[0], numberOfTracksRange[1])};

  void process(const o2::aod::Collision& collision, const o2::aod::Tracks& tracks)
  {
    hCollisionX->Fill(collision.posX());
    hCollisionY->Fill(collision.posY());
    hCollisionZ->Fill(collision.posZ());

    int nTracks(0);
    for (const auto& track : tracks) {
      nTracks++;
    }
    hNumberOfTracks->Fill(nTracks);
  }
};

/// Task to QA the kinematic properties of the tracks
struct QATrackingKine {
  o2fw::Configurable<bool> histogramAxisForMatplotlib{
    "histogramAxisForMatplotlib", false,
    "Sets the histograms title to be friendly to matplotlib instead of ROOT"};
  o2::qa::Features features = o2::qa::Features(histogramAxisForMatplotlib);
  o2fw::Configurable<int> nBinsPt{"nBinsPt", 100, "Number of bins for Pt"};
  std::array<double, 2> ptRange = {0, 10.};

  o2fw::Configurable<int> nBinsPhi{"nBinsPhi", 100, "Number of bins for Phi"};

  o2fw::Configurable<int> nBinsEta{"nBinsEta", 100, "Number of bins for the eta histogram."};
  std::array<double, 2> etaRange = {-6, 6};

  o2fw::OutputObj<TH1F> hPt{
    TH1F("pt", features.Title1D(features.Pt().Name()).c_str(), nBinsPt, ptRange[0], ptRange[1])};
  o2fw::OutputObj<TH1F> hEta{TH1F("eta", features.Title1D(features.Eta().NameNoUnit()).c_str(),
                                  nBinsEta, etaRange[0], etaRange[1])};
  o2fw::OutputObj<TH1F> hPhi{
    TH1F("phi", features.Title1D(features.Phi().Name()).c_str(), nBinsPhi, 0, 2 * M_PI)};

  void process(const o2::aod::Track& track)
  {
    hEta->Fill(track.eta());
    hPhi->Fill(track.phi());
    hPt->Fill(track.pt());
  }
};

/// Task to evaluate the tracking resolution (Pt, Eta, Phi and impact parameter)
struct QATrackingResolution {
  o2fw::Configurable<bool> setupHistogramAxisForMatplotlib{
    "setupHistogramAxisForMatplotlib", false,
    "Sets the histograms title to be friendly to matplotlib instead of ROOT"};
  o2::qa::Features features = o2::qa::Features(setupHistogramAxisForMatplotlib);

  o2fw::Configurable<int> nBinsPt{"nBinsPt", 100, "Number of bins for the transverse momentum"};
  std::array<double, 2> ptRange = {0, 10.};

  o2fw::Configurable<int> nBinsEta{"nBinsEta", 400, "Number of bins for the pseudorapidity"};
  std::array<double, 2> etaRange = {-3, 3};

  o2fw::Configurable<int> nBinsPhi{"nBinsPhi", 100, "Number of bins for Phi"};
  std::array<double, 2> phiRange = {0, 2 * M_PI};

  o2fw::Configurable<int> nBinsDeltaPt{"nBinsDeltaPt", 400,
                                       "Number of bins for the transverse momentum differences"};
  std::array<double, 2> deltaPtRange = {-1., 1.};

  o2fw::Configurable<int> nBinsDeltaPhi{"nBinsDeltaPhi", 100,
                                        "Number of bins for the azimuthal angle differences"};
  std::array<double, 2> deltaPhiRange = {-0.1, 0.1};

  o2fw::Configurable<int> nBinsDeltaEta{"nBinsDeltaEta", 100,
                                        "Number of bins for the pseudorapidity differences"};
  std::array<double, 2> deltaEtaRange = {-0.1, 0.1};

  o2fw::Configurable<int> nBinsImpactParameter{"nBinsImpactParameter", 1000,
                                               "Number of bins for the Impact parameter"};

  std::array<double, 2> impactParameterRange = {-1500, 1500};       // micrometer
  std::array<double, 2> impactParameterResolutionRange = {0, 1000}; // micrometer

  // Eta resolution
  o2fw::OutputObj<TH1F> etaDiffMCRec{TH1F("etaDiffMCReco",
                                          features.Title1D(features.Eta().MCRecDiff()).c_str(),
                                          nBinsDeltaEta, deltaEtaRange[0], deltaEtaRange[1])};

  o2fw::OutputObj<TH2F> etaDiffMCRecoVsEtaMC{
    TH2F("etaDiffMCRecoVsEtaMC",
         features.Title2D(features.Eta().MCRecDiff(), features.Eta().MC()).c_str(), nBinsDeltaEta,
         deltaEtaRange[0], deltaEtaRange[1], nBinsEta, etaRange[0], etaRange[1])};

  o2fw::OutputObj<TH2F> etaDiffMCRecoVsEtaReco{
    TH2F("etaDiffMCRecoVsEtaReco",
         features.Title2D(features.Eta().MCRecDiff(), features.Eta().Rec()).c_str(), nBinsDeltaEta,
         deltaEtaRange[0], deltaEtaRange[1], nBinsEta, etaRange[0], etaRange[1])};

  // Phi Resolution
  o2fw::OutputObj<TH1F> phiDiffMCRec{TH1F("phiDiffMCRec",
                                          features.Title1D(features.Phi().MCRecDiff()).c_str(),
                                          nBinsDeltaPhi, deltaPhiRange[0], deltaPhiRange[1])};

  // Pt Resolution
  o2fw::OutputObj<TH1F> ptDiffMCRec{TH1F("ptDiffMCRec",
                                         features.Title1D(features.Pt().MCRecDiff()).c_str(),
                                         nBinsDeltaPt, deltaPtRange[0], deltaPtRange[1])};

  o2fw::OutputObj<TH1F> ptResolution{
    TH1F("ptResolution", features.Title1D(features.Pt().RelativeMCRecDiff()).c_str(), nBinsDeltaPt,
         -1., 1.)};

  o2fw::OutputObj<TH2F> ptResolutionVsPt{
    TH2F("ptResolutionVsPt",
         features.Title2D(features.Pt().Rec(), features.Pt().RelativeMCRecDiff()).c_str(), nBinsPt,
         ptRange[0], ptRange[1], nBinsDeltaPt, 0., deltaPtRange[1])};

  o2fw::OutputObj<TH2F> ptResolutionVsEta{
    TH2F("ptResolutionVsEta",
         features.Title2D(features.Eta().Rec(), features.Pt().RelativeMCRecDiff()).c_str(),
         nBinsEta, etaRange[0], etaRange[1], nBinsDeltaPt, 0., deltaPtRange[1])};

  o2fw::OutputObj<TH2F> ptResolutionVsPhi{
    TH2F("ptResolutionVsPhi",
         features.Title2D(features.Phi().Rec(), features.Pt().RelativeMCRecDiff()).c_str(),
         nBinsPhi, phiRange[0], phiRange[1], nBinsDeltaPt, 0., deltaPtRange[1])};

  // Impact parameter RPhi
  o2fw::OutputObj<TH2F> impactParameterRPhiVsPt{
    TH2F("impactParameterRPhiVsPt",
         features.Title2D(features.Pt().Rec(), features.ImpactParameterRPhi().Name()).c_str(),
         nBinsPt, ptRange[0], ptRange[1], nBinsImpactParameter, impactParameterRange[0],
         impactParameterRange[1])};

  o2fw::OutputObj<TH2F> impactParameterRPhiVsEta{
    TH2F("impactParameterRPhiVsEta",
         features.Title2D(features.Eta().Rec(), features.ImpactParameterRPhi().Name()).c_str(),
         nBinsEta, etaRange[0], etaRange[1], nBinsImpactParameter, impactParameterRange[0],
         impactParameterRange[1])};

  o2fw::OutputObj<TH2F> impactParameterRPhiVsPhi{
    TH2F("impactParameterRPhiVsPhi",
         features.Title2D(features.Phi().Rec(), features.ImpactParameterRPhi().Name()).c_str(),
         nBinsPhi, phiRange[0], phiRange[1], nBinsImpactParameter, impactParameterRange[0],
         impactParameterRange[1])};

  o2fw::OutputObj<TH2F> impactParameterErrorRPhiVsPt{
    TH2F("impactParameterErrorRPhiVsPt",
         features.Title2D(features.Pt().Rec(), features.ImpactParameterRPhiError().Name()).c_str(),
         nBinsPt, ptRange[0], ptRange[1], nBinsImpactParameter, impactParameterResolutionRange[0],
         impactParameterResolutionRange[1])};

  o2fw::OutputObj<TH2F> impactParameterErrorRPhiVsEta{
    TH2F("impactParameterErrorRPhiVsEta",
         features.Title2D(features.Eta().Rec(), features.ImpactParameterRPhiError().Name()).c_str(),
         nBinsEta, etaRange[0], etaRange[1], nBinsImpactParameter,
         impactParameterResolutionRange[0], impactParameterResolutionRange[1])};

  o2fw::OutputObj<TH2F> impactParameterErrorRPhiVsPhi{
    TH2F("impactParameterErrorRPhiVsPhi",
         features.Title2D(features.Phi().Rec(), features.ImpactParameterRPhiError().Name()).c_str(),
         nBinsPhi, phiRange[0], phiRange[1], nBinsImpactParameter,
         impactParameterResolutionRange[0], impactParameterResolutionRange[1])};

  // Impact parameter Z
  o2fw::OutputObj<TH2F> impactParameterZVsPt{
    TH2F("impactParameterZVsPt",
         features.Title2D(features.Pt().Rec(), features.ImpactParameterZ().Name()).c_str(), nBinsPt,
         ptRange[0], ptRange[1], nBinsImpactParameter, impactParameterRange[0],
         impactParameterRange[1])};

  o2fw::OutputObj<TH2F> impactParameterZVsEta{
    TH2F("impactParameterZVsEta",
         features.Title2D(features.Eta().Rec(), features.ImpactParameterZ().Name()).c_str(),
         nBinsEta, etaRange[0], etaRange[1], nBinsImpactParameter, impactParameterRange[0],
         impactParameterRange[1])};
  o2fw::OutputObj<TH2F> impactParameterZVsPhi{
    TH2F("impactParameterZVsPhi",
         features.Title2D(features.Phi().Rec(), features.ImpactParameterZ().Name()).c_str(),
         nBinsPhi, phiRange[0], phiRange[1], nBinsImpactParameter, impactParameterRange[0],
         impactParameterRange[1])};

  o2fw::OutputObj<TH2F> impactParameterErrorZVsPt{
    TH2F("impactParameterErrorZVsPt",
         features.Title2D(features.Pt().Rec(), features.ImpactParameterZError().Name()).c_str(),
         nBinsPt, ptRange[0], ptRange[1], nBinsImpactParameter, impactParameterResolutionRange[0],
         impactParameterResolutionRange[1])};

  o2fw::OutputObj<TH2F> impactParameterErrorZVsEta{
    TH2F("impactParameterErrorVsEta",
         features.Title2D(features.Eta().Rec(), features.ImpactParameterZError().Name()).c_str(),
         nBinsEta, etaRange[0], etaRange[1], nBinsImpactParameter, impactParameterRange[0],
         impactParameterRange[1])};

  o2fw::OutputObj<TH2F> impactParameterErrorZVsPhi{
    TH2F("impactParameterErrorZVsPhi",
         features.Title2D(features.Phi().Rec(), features.ImpactParameterZError().Name()).c_str(),
         nBinsPhi, phiRange[0], phiRange[1], nBinsImpactParameter, impactParameterRange[0],
         impactParameterRange[1])};

  void process(
    const o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>::iterator& collision,
    const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::McTrackLabels>& tracks,
    const o2::aod::McParticles& mcParticles, const o2::aod::McCollisions& mcCollisions)
  {
    const o2df::VertexBase primaryVertex = getPrimaryVertex(collision);

    for (const auto& track : tracks) {
      const double deltaPt = track.label().pt() - track.pt();
      ptDiffMCRec->Fill(deltaPt);

      const double deltaPtOverPt = deltaPt / track.pt();
      ptResolution->Fill(deltaPtOverPt);
      ptResolutionVsPt->Fill(track.pt(), std::abs(deltaPtOverPt));
      ptResolutionVsEta->Fill(track.eta(), std::abs(deltaPtOverPt));
      ptResolutionVsPhi->Fill(track.phi(), std::abs(deltaPtOverPt));

      const double deltaEta = track.label().eta() - track.eta();
      etaDiffMCRec->Fill(deltaEta);
      etaDiffMCRecoVsEtaMC->Fill(deltaEta, track.label().eta());
      etaDiffMCRecoVsEtaReco->Fill(deltaEta, track.eta());

      const auto deltaPhi = track_utils::ConvertPhiRange(track.label().phi() - track.phi());
      phiDiffMCRec->Fill(deltaPhi);

      double impactParameterRPhi{-999.}, impactParameterRPhiError{-999.};
      double impactParameterZ{-999.}, impactParameterErrorZ{-999.};

      const bool propagate = track_utils::GetImpactParameterAndError(
        track, primaryVertex, impactParameterRPhi, impactParameterRPhiError, impactParameterZ,
        impactParameterErrorZ);

      if (propagate) {
        impactParameterRPhiVsPt->Fill(track.pt(), impactParameterRPhi);
        impactParameterRPhiVsEta->Fill(track.eta(), impactParameterRPhi);
        impactParameterRPhiVsPhi->Fill(track.phi(), impactParameterRPhi);

        impactParameterZVsPt->Fill(track.pt(), impactParameterZ);
        impactParameterZVsEta->Fill(track.eta(), impactParameterZ);
        impactParameterZVsPhi->Fill(track.phi(), impactParameterZ);

        impactParameterErrorRPhiVsPt->Fill(track.pt(), impactParameterRPhiError);
        impactParameterErrorRPhiVsEta->Fill(track.eta(), impactParameterRPhiError);
        impactParameterErrorRPhiVsPhi->Fill(track.phi(), impactParameterRPhiError);

        impactParameterErrorZVsPt->Fill(track.pt(), impactParameterErrorZ);
        impactParameterErrorZVsEta->Fill(track.eta(), impactParameterErrorZ);
        impactParameterErrorZVsPhi->Fill(track.phi(), impactParameterErrorZ);
      }
    }
  }
};

o2fw::WorkflowSpec defineDataProcessing(o2fw::ConfigContext const&)
{
  return o2fw::WorkflowSpec{
    o2fw::adaptAnalysisTask<QAGlobalObservables>("qa-global-observables"),
    o2fw::adaptAnalysisTask<QATrackingKine>("qa-tracking-kine"),
    o2fw::adaptAnalysisTask<QATrackingResolution>("qa-tracking-resolution")};
}
