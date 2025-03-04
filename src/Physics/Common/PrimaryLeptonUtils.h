//____________________________________________________________________________
/*!

\namespace  genie::utils

\brief      Common functions used for handling generation of the primary
            lepton, regardless of whether the relevant class inherits from
            PrimaryLeptonGenerator or not.

\author     Steven Gardiner <gardiner \at fnal.gov>
            Fermi National Accelerator Laboratory

\created    May 01, 2020

\cpright    Copyright (c) 2003-2024, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PRIMARY_LEPTON_UTILS_H
#define _PRIMARY_LEPTON_UTILS_H

#include "TVector3.h"
#include "TLorentzVector.h"

namespace genie {

class GHepRecord;

namespace utils {

  void SetPrimaryLeptonPolarization( GHepRecord* ev );

  TVector3 CalculatePolarizationVectorInTargetRestFrame(
      const TLorentzVector & neutrinoMomTRF,
      const TLorentzVector & leptonMomTRF, 
      bool isLeftPolarized,
      double M,
      double W1,
      double W2,
      double W3,
      double W4,
      double W5
  );

} // utils   namespace
} // genie   namespace

#endif // _PRIMARY_LEPTON_UTILS_H
