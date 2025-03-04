//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "TVector3.h"

#include "Physics/Common/PrimaryLeptonUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
  
using namespace genie;
using namespace genie::utils;


void genie::utils::SetPrimaryLeptonPolarization( GHepRecord * ev )
{

  // get the final state primary lepton
  GHepParticle * fsl = ev->FinalStatePrimaryLepton();
  if ( !fsl ) {
    LOG("LeptonicVertex", pERROR)
      << "Final state lepton not set yet! \n" << *ev;
    return;
  }

  //-- Get the interaction
  Interaction * interaction = ev->Summary();

  //-- Access cross section algorithm for running thread    //TODO Would it be better here if the polarization was stored to Interaction/Kinematics object?
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  const XSecAlgorithmI * xsec_alg = evg->CrossSectionAlg();
  TVector3 fsl_pol = xsec_alg->FinalLeptonPolarization(interaction);
  fsl->SetPolarization(fsl_pol);

  LOG("LeptonicVertex", pINFO)
    << "Setting polarization angles for particle: " << fsl->Name();

  if ( fsl->PolzIsSet() ) {
    LOG("LeptonicVertex", pINFO)
      << "Polarization : "  << fsl->GetPolarization().x() << "," << fsl->GetPolarization().y() << "," << fsl->GetPolarization().y();
  }
  
}
//____________________________________________________________________________
TVector3 genie::utils::CalculatePolarizationVectorInTargetRestFrame(
  const TLorentzVector & neutrinoMomTRF,
  const TLorentzVector & leptonMomTRF, 
  bool isLeftPolarized,
  double M,
  double W1,
  double W2,
  double W3,
  double W4,
  double W5
)
{
  /*
    Calculate the primary lepton polarization vector in the target rest frame

    References:
      [1] https://arxiv.org/pdf/hep-ph/0305324
  */

  //TODO not the same as the version from the Russian PR - why?

  // Pre-compute some important quantities
  double Enu = neutrinoMomTRF.E(); // Neutrino energy
  double El = leptonMomTRF.E(); // Primary lepton energy
  double ml = leptonMomTRF.M(); // Primary lepton mass
  double pl = leptonMomTRF.P(); // Primary lepton momentum
  double theta = leptonMomTRF.Angle(neutrinoMomTRF.Vect()); // nu-lepton angle
  double costheta = cos(theta);

  // Some terms if the following calculation flip sign based on nu vs nubar
  int sign = isLeftPolarized ? -1 : 1; // -1 for neutrinos, +1 for antineutrinos

  // Calculate denomenator term
  double F = (
      (2*W1 + pow(ml,2) / pow(M,2) * W4) * (El - pl * costheta)
      + W2 * (El + pl * costheta)
      + -sign * (W3 / M) * ((Enu * El) + pow(pl,2) - (Enu + El) * pl * costheta)
      - pow(ml,2) / M * W5
  );

  // Calculate polarization transverse component, [1] eqn 25a (s_x)
  // Note that a factor 1/2 has been removed to give the polarization vector (max magnitude = 1) rather than spin vector (max magnitude = 1/2) 
  double pol_t = sign * ml*sin(theta) * (2*W1 - W2 + ( -sign * (Enu / M) * W3)  - ( (pow(ml,2)/pow(M,2)) * W4) + ((El/M) * W5)) / F;

  // Calculate polarization longitudinal component: [1] eqn 25c (s_z)
  // Note that a factor 1/2 has been removed to give the polarization vector (max magnitude = 1) rather than spin vector (max magnitude = 1/2) 
  double pol_l = sign * (
      (2*W1 - ((pow(ml,2) / pow(M,2)) * W4)) * (pl - (El * costheta))
      + W2 * (pl + (El * costheta))
      +  (-sign) * (W3 / M) * ((Enu + El) * pl - (Enu * El + pow(pl,2)) * costheta)
      - (pow(ml,2) / M) * W5 * costheta
  ) / F;

  // Return the vector object, following the convention from [1] 25
  return TVector3(pol_t, 0., pol_l);

}
//____________________________________________________________________________
