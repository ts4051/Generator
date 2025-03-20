//____________________________________________________________________________
/*!
TODO
*/
//____________________________________________________________________________

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

// #if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
// #include <fenv.h> // for `feenableexcept`
// #endif

#include <TLorentzVector.h>

#include "Framework/Interaction/Interaction.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Interaction/InteractionType.h"
#include "Framework/Interaction/ScatteringType.h"
#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Utils/CmdLnArgParser.h"

#include "Physics/DeepInelastic/XSection/QPMDISPXSec.h"
#include "Physics/DeepInelastic/XSection/KNOTunedQPMDISPXSec.h"

#include "Physics/Resonance/XSection/BergerSehgalRESPXSec2014.h"
#include "Framework/ParticleData/BaryonResonance.h"


using namespace genie;
// using namespace genie::controls;


//
// Globals
//

double M_tau = 1.7769; // tau mass


//
// Helper functions
//

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

   // from https://stackoverflow.com/questions/27028226/python-linspace-in-c

   std::vector<double> linspaced;

   double start = static_cast<double>(start_in);
   double end = static_cast<double>(end_in);
   double num = static_cast<double>(num_in);

   if (num == 0) { return linspaced; }
   if (num == 1) {
      linspaced.push_back(start);
      return linspaced;
   }

   double delta = (end - start) / (num - 1);

   for(int i=0; i < num-1; ++i) {
      linspaced.push_back(start + delta * i);
   }
   linspaced.push_back(end); // I want to ensure that start and end are exactly the same as the input
   return linspaced;
}


//
// Main
//


//____________________________________________________________________________
int main(int argc, char ** argv)
{

   std::cout << "Start..." << std::endl;


   //
   // Parse args
   //

   // Boilerplate setup
   // RunOpt::Instance()->EnableBareXSecPreCalc(true); //TODO what is this?

   // Common command line parsing
   RunOpt::Instance()->ReadFromCommandLine(argc, argv);

   // Parsing specific to this app
   CmdLnArgParser parser(argc,argv);

   // Neutrino energy
   double Ev = -999.;
  if( parser.OptionExists('e') ) {
      Ev = parser.ArgAsDouble('e');
      std::cout << "Neutrino energy = " << Ev << " GeV" << std::endl;
   }
   else {
      std::cerr << "Must specify energy argument : -e <energy>" << std::endl;
      exit(-1);
   }


   //
   // Steering
   //

   // Target  //TODO arg?
   int nucleus_pdg = 1000080160 ; // O16=1000080160, d=1000010020

   // Define scan parameters
   std::vector<std::string> int_type_scan{"dis", "res"};
   std::vector<double> El_scan = linspace( std::max(M_tau, 0.1*Ev), 0.9*Ev, 200); // Don't go below the tau mass
   // std::vector<double> theta_deg_scan{1e-15, 5., 10.}; // Problems with theta = 0, so just use a very small number instead
   std::vector<double> theta_deg_scan = linspace(1e-15, 30., 61); // Problems with theta = 0, so just use a very small number instead
   std::vector<int> neutrino_pdg_scan{16, -16};
   std::vector<int> nucleon_pdg_scan{2212, 2112}; // proton, neutron
   std::vector<int> quark_pdg_scan{-2, -1, +1, +2}; // Seems to make no difference
   std::vector<bool> sea_quark_scan{false} ; //sea_quark_scan{false, true}; // Seems to make no difference
   std::vector<Resonance_t> resonance_scan{kP33_1232} ; //TODO try others

   // Open output file
   std::string outfile_name = "gscanvars_" + std::to_string((int)Ev) + "GeV.csv"; //TODO target
   std::ofstream outfile; 
   outfile.open(outfile_name); 

   // Write outfile header
   outfile << "int_type" <<  "," << "nucleon_pdg" << "," << "quark_pdg" << ","<< "sea_quark" << "," << "resonance" << ",";
   outfile << "nubar" << ",";
   outfile << "nu_px" << "," << "nu_py" << "," << "nu_pz" << "," << "nu_E" << ",";
   outfile << "l_px" << "," << "l_py" << "," << "l_pz" << "," << "l_E"  << ",";
   outfile << "theta_deg" << "," << "cos_theta_p" << ",";
   outfile << "xsec" << ",";
   outfile << "x" << "," << "y" << "," << "Q2" << "," << "W" << ",";
   outfile << "pol_x" << "," <<  "pol_y" << "," <<  "pol_z";
   outfile << std::endl;



   //
   // Create Xsec models
   //

   // DIS model
   KNOTunedQPMDISPXSec dis_xsec_model;
   dis_xsec_model.Configure("Default"); //TODO Get this from config?
   ProcessInfo dis_proc_info(kScDeepInelastic, kIntWeakCC);

   // RES model
   BergerSehgalRESPXSec2014 res_xsec_model;
   res_xsec_model.Configure("NoPauliBlock"); //TODO Get this from config?
   ProcessInfo res_proc_info(kScResonant, kIntWeakCC);

   // Get RES/DIS transition
   bool res_dis_joining = res_xsec_model.GetConfig().GetBool("UseDRJoinScheme");
   assert(res_dis_joining);
   double Wcut = res_xsec_model.GetConfig().GetDouble("Wcut");

   //TODO force Wcut to 1.4 GeV ?



   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?
   //TODO KLN RES ?


   //
   // Scan
   //

   // Scan integration type
   for (std::string int_type : int_type_scan) {

      // Choose interaction model
      XSecAlgorithmI* xsec_model = NULL;
      if(int_type == "dis") {
         xsec_model = &dis_xsec_model;
      }
      else if(int_type == "res") {
        xsec_model = &res_xsec_model; 
      }
      else {
         std::cerr << "Unknown model : " << int_type << std::endl;
         exit(-1);
      }

      // Loop over target nucelon
      for (double nucleon_pdg : nucleon_pdg_scan) {

         // Loop over target quark
         // This is only relevent for DIS
         auto quark_pdg_scan_to_use = (int_type == "dis") ? quark_pdg_scan : std::vector<int>{0};
         auto sea_quark_scan_to_use = (int_type == "dis") ? sea_quark_scan : std::vector<bool>{false};
         for (double quark_pdg : quark_pdg_scan_to_use) {
            for (double sea_quark : sea_quark_scan_to_use) {

               // Loop over baryon
               // This is over relevent for resonance events
               auto resonance_scan_to_use = (int_type == "res") ? resonance_scan : std::vector<Resonance_t>{kNoResonance};
               for (Resonance_t resonance : resonance_scan_to_use) {

                  // Nu/nubar loop
                  for (double neutrino_pdg : neutrino_pdg_scan) {

                     // Scan the interaction parameters
                     for (double theta_deg : theta_deg_scan) {
                        for (double El : El_scan) {

                           std::cout << std::endl << "---------------------------------------------------------" << std::endl << std::endl << std::endl;

                           // Derived variables
                           double theta = theta_deg * (3.14159265359 / 180.);
                           TLorentzVector nu_p4(0., 0., Ev, Ev); // px,py,pz,E (GeV)


                           //
                           // Create interaction
                           //

                           // Define initial state
                           InitialState ini_state(nucleus_pdg, neutrino_pdg);
                           ini_state.SetProbeP4(nu_p4);

                           // Define  target
                           Target* target = ini_state.TgtPtr();
                           target->SetHitNucPdg(nucleon_pdg);
                           if(int_type == "dis") {
                              target->SetHitQrkPdg(quark_pdg);
                              target->SetHitSeaQrk(sea_quark);
                           }

                           // Define interaction process
                           ProcessInfo* proc_info = NULL;
                           if(int_type == "dis") {
                              proc_info = &dis_proc_info;
                           }
                           else if(int_type == "res") {
                              proc_info = &res_proc_info;
                           }

                           // Create interaction object
                           Interaction interaction(ini_state, *proc_info);

                           // Set ExclTag (specify resonance)
                           XclsTag excl_tag;
                           excl_tag.SetResonance(resonance);
                           interaction.SetExclTag(excl_tag);


                           //
                           // Define kinematics
                           //

                           Kinematics kinematics; 

                           // Define lepton
                           int lepton_pdg = neutrino_pdg < 0 ? (neutrino_pdg+1) : (neutrino_pdg-1);
                           double pl = sqrt( pow(El, 2.) - pow(M_tau, 2.));
                           // TLorentzVector l_p4(0., 0., pl, El);
                           // l_p4.RotateX(theta);
                           TLorentzVector l_p4(pl*sin(theta), 0., pl*cos(theta), El); // [1] eqn 4
                           kinematics.SetFSLeptonP4(l_p4);

                           // Calc kinematic variables
                           double y = 1. - (El / Ev); //TODO This only works if the target is at rest (which it is currently in this test)
                           auto q = nu_p4 - l_p4;
                           double Q2 = -q.Mag2();
                           double M = target->HitNucMass();
                           std::cout << ", y = " << y << ", Q2 = " << Q2 << std::endl;
                           double x = utils::kinematics::Q2YtoX(Ev, M, Q2, y);
                           double W = utils::kinematics::XYtoW(Ev, M, x, y);

                           // double M = target->HitNucMass();
                           // double W = utils::kinematics::W(&interaction); //TODO doesn't work, probably because I haven't defined the hadronic states
                           // double Q2 = utils::kinematics::Q2(&interaction);
                           // double x; double y;
                           // utils::kinematics::WQ2toXY(Ev, M, W, Q2, x, y);

                           std::cout << "x = " << x << ", y = " << y << ", Q2 = " << Q2 << ", W = " << W << std::endl;

                           // Set kinematic variables
                           kinematics.Setx(x); //TODO selected=true ?
                           kinematics.Sety(y);
                           kinematics.SetQ2(Q2);
                           kinematics.SetW(W);

                           // Pass kinematics to interaction 
                           interaction.SetKine(kinematics);

                           // Skip scan point if not valid
                           if (!xsec_model->ValidProcess(&interaction)) continue;
                           if (!xsec_model->ValidKinematics(&interaction)) continue;

                           // Hnadle RES-DIS transition
                           if(int_type == "res") {
                              if(W >= Wcut) continue;
                           }
                           else if(int_type == "dis") {
                              if(W < Wcut) continue;
                           }


                           //
                           // Calculate polarization
                           //

                           // First need to run xsec cslc
                           double xsec = xsec_model->XSec(&interaction, kPSxyfE); //TODO which KinePhaseSpace_t ?

                           // Calculate polarization, but only if xsec calc was successful
                           // TVector3 pol(-9999., -9999., -9999.);
                           // if(xsec != 0.) {
                           //    pol = xsec_model->FinalLeptonPolarization(&interaction); 
                           // }
                           TVector3 pol = xsec_model->FinalLeptonPolarization(&interaction); 

                           // Get theta_p
                           double cos_theta_p = cos( l_p4.Angle(pol) );

                           // Write to file
                           outfile << int_type << "," << nucleon_pdg << "," << quark_pdg << "," << sea_quark << "," << (int)resonance << ",";
                           outfile << (neutrino_pdg < 0) << ",";
                           outfile << nu_p4.Px() << "," << nu_p4.Py() << "," << nu_p4.Pz() << "," << nu_p4.E() << ",";
                           outfile << l_p4.Px() << "," << l_p4.Py() << "," << l_p4.Pz() << "," << l_p4.E()  << ",";
                           outfile << theta_deg << "," << cos_theta_p << ",";
                           outfile << xsec << ",";
                           outfile << x << "," << y << "," << Q2 << "," << W << ",";
                           outfile << pol.x() << "," <<  pol.y() << "," <<  pol.z();
                           outfile << std::endl;

                        }
                     }
                  }
               }
            }
         }
      }
   }

   std::cout << std::endl << "Done! Results in " << outfile_name << std::endl << std::endl;
   outfile.close(); 


}
