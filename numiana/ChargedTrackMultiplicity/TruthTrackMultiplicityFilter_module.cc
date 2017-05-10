////////////////////////////////////////////////////////////////////////
// Class:       TruthTrackMultiplicityFilter
// Plugin Type: filter (art v2_05_00)
// File:        TruthTrackMultiplicityFilter_module.cc
//
// Generated at Tue May  9 07:25:47 2017 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>
#include <iostream>

#include "TTree.h"

#include "lardataobj/MCBase/MCTrack.h"
#include "nusimdata/SimulationBase/MCTruth.h"


namespace ana {
  class TruthTrackMultiplicityFilter;
}


class ana::TruthTrackMultiplicityFilter : public art::EDFilter {
public:
  explicit TruthTrackMultiplicityFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TruthTrackMultiplicityFilter(TruthTrackMultiplicityFilter const &) = delete;
  TruthTrackMultiplicityFilter(TruthTrackMultiplicityFilter &&) = delete;
  TruthTrackMultiplicityFilter & operator = (TruthTrackMultiplicityFilter const &) = delete;
  TruthTrackMultiplicityFilter & operator = (TruthTrackMultiplicityFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;

private:

  art::InputTag fMCTruthTag;
  art::InputTag fMCTrackTag;

  int           fMinNPrimaryMCTracks;
  float         fMinPrimaryMCTrackLength;

  //int           fMinNSecVtx;
  //float         fMinSecVtxDist;
  //float         fMaxSecVtxDist;

  bool          fVerbose;
  
  TTree* fAnaTree;
  unsigned int run;
  unsigned int event;
  int nu_pdg;
  int ccnc;
  int mode;
  float q2;
  float nu_energy;
  float lep_energy;
  float vtx_x;
  float vtx_y;
  float vtx_z;
  int ntrks_1cm;
  int ntrks_5cm;
  int ntrks_10cm;
  int ntrks_20cm;
  void VarInit();
};


ana::TruthTrackMultiplicityFilter::TruthTrackMultiplicityFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

bool ana::TruthTrackMultiplicityFilter::filter(art::Event & ev)
{
  VarInit();

  if(fVerbose){
    std::cout << "Processing "
	      << "Run " << ev.run() << ", "
	      << "Event " << ev.event() << ", " 
	      << "Time " << ev.time().timeHigh() << std::endl;
  }
  
  run = ev.run();
  event = ev.event();
  
  auto const& mctrack_handle = ev.getValidHandle< std::vector<sim::MCTrack> >(fMCTrackTag);
  auto const& mctruth_handle = ev.getValidHandle< std::vector<simb::MCTruth> >(fMCTruthTag);
  
  TVector3 nu_pos;
  int n_tracks_pass=0;

  for(auto const& mctruth : *mctruth_handle){

    if(mctruth.Origin()!=simb::Origin_t::kBeamNeutrino)
      continue;

    auto const& nu = mctruth.GetNeutrino();
    ccnc = nu.CCNC();
    mode = nu.Mode();
    q2 = nu.QSqr();
    nu_energy = nu.Nu().E();
    nu_pdg = nu.Nu().PdgCode();
    vtx_x = nu.Nu().Vx();
    vtx_y = nu.Nu().Vy();
    vtx_z = nu.Nu().Vz();
    
    if(fVerbose) std::cout << nu_pdg << " " << nu.Nu().Vx() << "," << nu.Nu().Vy() << "," << nu.Nu().Vz() << std::endl;
    nu_pos = nu.Nu().Position().Vect();

    if(mode==0) //CC
      lep_energy = nu.Lepton().E();
  }

  ntrks_1cm=0; ntrks_5cm=0; ntrks_10cm=0; ntrks_20cm=0;
  float length;
  for(auto const& trk : *mctrack_handle){
    
    if(fVerbose)
      std::cout << "\t" << trk.Origin() << " " << trk.MotherPdgCode() << " "
		<< trk.Start().X() << "," << trk.Start().Y() << "," << trk.Start().Z() << " " 
		<< (trk.End().Position().Vect() - trk.Start().Position().Vect()).Mag() << std::endl;
    
    if(trk.Origin()!=simb::Origin_t::kBeamNeutrino ||	 
       (trk.Start().Position().Vect() - nu_pos).Mag() > 0.3 ) //wire pitch/resolution cut
      continue;
    
    length = (trk.End().Position().Vect() - trk.Start().Position().Vect()).Mag();
    if(length>1)
      ntrks_1cm++;
    if(length>5)
      ntrks_5cm++;
    if(length>10)
      ntrks_10cm++;
    if(length>20)
      ntrks_20cm++;

    if(length > fMinPrimaryMCTrackLength)
      n_tracks_pass++;
  }

  if(n_tracks_pass>fMinNPrimaryMCTracks)
    return true;
  
  return false;
}

void ana::TruthTrackMultiplicityFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fMCTruthTag = p.get<art::InputTag>("MCTruthTag");
  fMCTrackTag = p.get<art::InputTag>("MCTrackTag");

  fMinNPrimaryMCTracks     = p.get<int>("MinNPrimaryMCTracks",-1);
  fMinPrimaryMCTrackLength = p.get<float>("MinPrimaryMCTrackLength",0.0);

  //fMinNSecVtx     = p.get<int>("MinNSecVtx",-1);
  //fMinNSecVtxDist = p.get<float>("MinNSecVtxDist",0.0);
  //fMaxNSecVtxDist = p.get<float>("MaxNSecVtxDist",999999.999);

  fVerbose = p.get<bool>("Verbose",false);
}

void ana::TruthTrackMultiplicityFilter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fAnaTree = tfs->make<TTree>("ana","Nu Interaction Tree");//,"run:ev:pdg:ccnc:mode:q2:energy:vtx_x:vtx_y:vtx_z:ntrks_1cm:ntrks_10cm:ntrks_20cm");
  fAnaTree->Branch("run",&run,"run/I");
  fAnaTree->Branch("event",&event,"event/I");
  fAnaTree->Branch("nu_pdg",&nu_pdg,"nu_pdg/I");
  fAnaTree->Branch("ccnc",&ccnc,"ccnc/I");
  fAnaTree->Branch("mode",&mode,"mode/I");
  fAnaTree->Branch("q2",&q2,"q2/F");
  fAnaTree->Branch("nu_energy",&nu_energy,"nu_energy/F");
  fAnaTree->Branch("lep_energy",&lep_energy,"lep_energy/F");
  fAnaTree->Branch("vtx_x",&vtx_x,"vtx_x/F");
  fAnaTree->Branch("vtx_y",&vtx_y,"vtx_y/F");
  fAnaTree->Branch("vtx_z",&vtx_z,"vtx_z/F");
  fAnaTree->Branch("ntrks_1cm",&ntrks_1cm,"ntrks_1cm/I");
  fAnaTree->Branch("ntrks_5cm",&ntrks_5cm,"ntrks_5cm/I");
  fAnaTree->Branch("ntrks_10cm",&ntrks_10cm,"ntrks_10cm/I");
  fAnaTree->Branch("ntrks_20cm",&ntrks_20cm,"ntrks_20cm/I");
}

void ana::TruthTrackMultiplicityFilter::VarInit()
{
 run=-1;
 event=-1;
 nu_pdg = -99999;
 ccnc = -1;
 mode = -1;
 q2 = -999;
 nu_energy = -999;
 lep_energy = -999;
 vtx_x = -999999;
 vtx_y = -999999;
 vtx_z = -999999;
 ntrks_1cm = -1;
 ntrks_5cm = -1;
 ntrks_10cm = -1;
 ntrks_20cm = -1;


}

DEFINE_ART_MODULE(ana::TruthTrackMultiplicityFilter)
