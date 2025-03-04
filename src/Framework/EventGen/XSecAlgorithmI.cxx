//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI() :
Algorithm()
{

}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
XSecAlgorithmI::~XSecAlgorithmI()
{

}
//___________________________________________________________________________
bool XSecAlgorithmI::ValidKinematics(const Interaction* interaction) const
{
// can offer common implementation for all concrete x-section models because
// the input interaction is aware of its kinematic limits

  if ( interaction->TestBit(kISkipKinematicChk) ) return true;

  const KPhaseSpace& kps = interaction->PhaseSpace();

  if ( ! kps.IsAboveThreshold() ) {
     LOG("XSecBase", pINFO)  << "*** Below energy threshold";
     return false;
  }
  if ( ! kps.IsAllowed() ) {
     LOG("XSecBase", pINFO)  << "*** Not in allowed kinematical space";
     return false;
  }
  return true;
}
//___________________________________________________________________________
TVector3 XSecAlgorithmI::FinalLeptonPolarization (const Interaction* i) const
{
    /*
    A generic final lepton polarization calculation.
    It assumes a massless lepton, so gives incorect results for tau in particular.
    Individual models may override this with more accurate calculations.
    */

    // Not implementation at present for EM interactions
    if ( i->ProcInfo().IsEM() ) {
        LOG("XSecBase", pWARN) << "Polarization not yet implemented for EM processes. Setting it to zero.";
        return TVector3(0., 0., 0.);
    }

    // For weak interactions, simply set:
    //   - Left-Handed (spin and momentum anti-parallel for particles)
    //   - Right-Handed (spin and momentum parallel for anti-particles)
    const Kinematics & kinematics = i->Kine();
    const TLorentzVector lepton_4p = kinematics.FSLeptonP4();
    TVector3 pol = lepton_4p.Vect().Unit(); // Start from momentum direction
    int lepton_pdg = i->FSPrimLeptonPdg();
    if ( pdg::IsNeutrino(lepton_pdg) || pdg::IsElectron(lepton_pdg) || pdg::IsMuon(lepton_pdg) || pdg::IsTau(lepton_pdg) ) {
        pol = -1. * pol; //Flip direction for particles (anti-parallel)
    }

    return pol;
}
//___________________________________________________________________________
