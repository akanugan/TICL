// -*- C++ -*-
//
// Package:    NTupler
// Class:      ClusterNtupler
//
/**\class RecHitTupler RecHitTupler.cc NTupler/plugins/ClusterNtupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Philipp Zehetner
// Based on example written by Thorben Quast
//
//


#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class Ntupler : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Ntupler(const edm::ParameterSet&);
  ~Ntupler();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void endJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // some options
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  const edm::EDGetTokenT<TICLGraph> ticl_graph_token_;
  const edm::EDGetTokenT<std::vector<TICLCandidate>> ticl_candidates_token_;
  const edm::EDGetTokenT<std::vector<reco::Track>> tracks_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_merged_token_;

  // Output tree
  TTree* tree_;

  void clearVariables();

  unsigned int event_index;

  // Variables for branches
  unsigned int ev_event_;
  unsigned int ntracksters_;
  unsigned int nclusters_;

  std::vector<float_t> trackster_time;
  std::vector<float_t> trackster_timeError;
  std::vector<float_t> trackster_regressed_energy;
  std::vector<float_t> trackster_raw_energy;
  std::vector<float_t> trackster_raw_em_energy;
  std::vector<float_t> trackster_raw_pt;
  std::vector<float_t> trackster_raw_em_pt;
  std::vector<float_t> trackster_barycenter_x;
  std::vector<float_t> trackster_barycenter_y;
  std::vector<float_t> trackster_barycenter_z;
  std::vector<float_t> trackster_EV1;
  std::vector<float_t> trackster_EV2;
  std::vector<float_t> trackster_EV3;
  std::vector<float_t> trackster_eVector0_x;
  std::vector<float_t> trackster_eVector0_y;
  std::vector<float_t> trackster_eVector0_z;
  std::vector<float_t> trackster_sigmaPCA1;
  std::vector<float_t> trackster_sigmaPCA2;
  std::vector<float_t> trackster_sigmaPCA3;
  std::vector<std::vector<float_t>> trackster_id_probabilities;
  std::vector<std::vector<uint16_t> > trackster_vertices_indexes;
  std::vector<std::vector<float_t> > trackster_vertices_x;
  std::vector<std::vector<float_t> > trackster_vertices_y;
  std::vector<std::vector<float_t> > trackster_vertices_z;
  std::vector<std::vector<float_t> > trackster_vertices_energy;
  std::vector<std::vector<float_t> > trackster_vertices_correctedEnergy;
  std::vector<std::vector<float_t> > trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float_t> > trackster_vertices_multiplicity;
  
  // from TICLGraph
  std::vector<std::vector<uint16_t>> node_linked_inners;
  std::vector<std::vector<uint16_t>> node_linked_outers;
  std::vector<bool> isRootTrackster;

  // from TICLCandidate, product of linking
  size_t nCandidates;
  std::vector<int> candidate_charge;
  std::vector<int> candidate_pdgId;
  std::vector<float_t> candidate_energy;
  std::vector<double> candidate_px;
  std::vector<double> candidate_py;
  std::vector<double> candidate_pz;
  std::vector<float_t> candidate_time;
  std::vector<float_t> candidate_time_err;
  std::vector<std::vector<float_t>> candidate_id_probabilities;
  std::vector<std::vector<uint16_t>> tracksters_in_candidate;
  std::vector<uint16_t> track_in_candidate;

  // merged tracksters
  size_t nTrackstersMerged;
  std::vector<float_t> tracksters_merged_barycenter_x;
  std::vector<float_t> tracksters_merged_barycenter_y;
  std::vector<float_t> tracksters_merged_barycenter_z;
  std::vector<float_t> tracksters_merged_EV1;
  std::vector<float_t> tracksters_merged_EV2;
  std::vector<float_t> tracksters_merged_EV3;
  std::vector<float_t> tracksters_merged_eVector0_x;
  std::vector<float_t> tracksters_merged_eVector0_y;
  std::vector<float_t> tracksters_merged_eVector0_z;
  std::vector<float_t> tracksters_merged_sigmaPCA1;
  std::vector<float_t> tracksters_merged_sigmaPCA2;
  std::vector<float_t> tracksters_merged_sigmaPCA3;
  std::vector<std::vector<float_t>> tracksters_merged_id_probabilities;

  std::vector<uint32_t> cluster_seedID;
  std::vector<float_t> cluster_energy;
  std::vector<float_t> cluster_correctedEnergy;
  std::vector<float_t> cluster_correctedEnergyUncertainty;
  std::vector<float_t> cluster_position_x;
  std::vector<float_t> cluster_position_y;
  std::vector<float_t> cluster_position_z;
  std::vector<float_t> cluster_position_eta;
  std::vector<float_t> cluster_position_phi;
  std::vector<uint32_t> cluster_number_of_hits;
  std::vector<uint32_t> cluster_hitsAndFractions_ID;
  std::vector<uint32_t> cluster_hitsAndFractions_value;

  TTree* trackster_tree_;
  TTree* cluster_tree_;
  TTree* graph_tree_;
  TTree* candidate_tree_;
  TTree* tracksters_merged_tree_;
};

void Ntupler::clearVariables() {
  // event info
  ev_event_ = 0;
  ntracksters_ = 0;
  nclusters_ = 0;

  trackster_time.clear();
  trackster_timeError.clear();
  trackster_regressed_energy.clear();
  trackster_raw_energy.clear();
  trackster_raw_em_energy.clear();
  trackster_raw_pt.clear();
  trackster_raw_em_pt.clear();
  trackster_barycenter_x.clear();
  trackster_barycenter_y.clear();
  trackster_barycenter_z.clear();
  trackster_EV1.clear();
  trackster_EV2.clear();
  trackster_EV3.clear();
  trackster_eVector0_x.clear();
  trackster_eVector0_y.clear();
  trackster_eVector0_z.clear();
  trackster_sigmaPCA1.clear();
  trackster_sigmaPCA2.clear();
  trackster_sigmaPCA3.clear();
  trackster_id_probabilities.clear();
  trackster_vertices_indexes.clear();
  trackster_vertices_x.clear();
  trackster_vertices_y.clear();
  trackster_vertices_z.clear();
  trackster_vertices_energy.clear();
  trackster_vertices_correctedEnergy.clear();
  trackster_vertices_correctedEnergyUncertainty.clear();

  node_linked_inners.clear();
  node_linked_outers.clear();
  isRootTrackster.clear();
  
  nCandidates = 0;
  candidate_charge.clear();
  candidate_pdgId.clear();
  candidate_energy.clear();
  candidate_px.clear();
  candidate_py.clear();
  candidate_pz.clear();
  candidate_time.clear();
  candidate_time_err.clear();
  candidate_id_probabilities.clear();
  tracksters_in_candidate.clear();
  track_in_candidate.clear();

  nTrackstersMerged = 0;
  tracksters_merged_barycenter_x.clear();
  tracksters_merged_barycenter_y.clear();
  tracksters_merged_barycenter_z.clear();
  tracksters_merged_EV1.clear();
  tracksters_merged_EV2.clear();
  tracksters_merged_EV3.clear();
  tracksters_merged_eVector0_x.clear();
  tracksters_merged_eVector0_y.clear();
  tracksters_merged_eVector0_z.clear();
  tracksters_merged_sigmaPCA1.clear();
  tracksters_merged_sigmaPCA2.clear();
  tracksters_merged_sigmaPCA3.clear();
  tracksters_merged_id_probabilities.clear();

  cluster_seedID.clear();
  cluster_energy.clear();
  cluster_correctedEnergy.clear();
  cluster_correctedEnergyUncertainty.clear();
  cluster_position_x.clear();
  cluster_position_y.clear();
  cluster_position_z.clear();
  cluster_position_eta.clear();
  cluster_position_phi.clear();
  cluster_number_of_hits.clear();
  cluster_hitsAndFractions_ID.clear();
  cluster_hitsAndFractions_value.clear();
};

Ntupler::Ntupler(const edm::ParameterSet& ps)
    : tracksters_token_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("trackstersclue3d"))),
      layer_clusters_token_(
          consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerClusters"))),
      ticl_graph_token_(consumes<TICLGraph>(ps.getParameter<edm::InputTag>("ticlgraph"))),
      ticl_candidates_token_(consumes<std::vector<TICLCandidate>>(ps.getParameter<edm::InputTag>("ticlcandidates"))),
      tracks_token_(consumes<std::vector<reco::Track>>(ps.getParameter<edm::InputTag>("tracks"))),
      tracksters_merged_token_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("trackstersmerged"))){
      };

Ntupler::~Ntupler() { clearVariables(); };

void Ntupler::beginJob() {
  // Define tree and branches
  edm::Service<TFileService> fs;
  trackster_tree_ = fs->make<TTree>("tracksters", "TICL tracksters");
  cluster_tree_ = fs->make<TTree>("clusters", "TICL tracksters");
  graph_tree_ = fs->make<TTree>("graph", "TICL graph");
  candidate_tree_ = fs->make<TTree>("candidates", "TICL candidates");
  tracksters_merged_tree_ = fs->make<TTree>("trackstersMerged", "TICL tracksters merged");

  trackster_tree_->Branch("event", &ev_event_);
  trackster_tree_->Branch("NClusters", &nclusters_);
  trackster_tree_->Branch("NTracksters", &ntracksters_);
  trackster_tree_->Branch("time", &trackster_time);
  trackster_tree_->Branch("timeError", &trackster_timeError);
  trackster_tree_->Branch("regressed_energy", &trackster_regressed_energy);
  trackster_tree_->Branch("raw_energy", &trackster_raw_energy);
  trackster_tree_->Branch("raw_em_energy", &trackster_raw_em_energy);
  trackster_tree_->Branch("raw_pt", &trackster_raw_pt);
  trackster_tree_->Branch("raw_em_pt", &trackster_raw_em_pt);
  trackster_tree_->Branch("barycenter_x", &trackster_barycenter_x);
  trackster_tree_->Branch("barycenter_y", &trackster_barycenter_y);
  trackster_tree_->Branch("barycenter_z", &trackster_barycenter_z);
  trackster_tree_->Branch("EV1", &trackster_EV1);
  trackster_tree_->Branch("EV2", &trackster_EV2);
  trackster_tree_->Branch("EV3", &trackster_EV3);
  trackster_tree_->Branch("eVector0_x", &trackster_eVector0_x);
  trackster_tree_->Branch("eVector0_y", &trackster_eVector0_y);
  trackster_tree_->Branch("eVector0_z", &trackster_eVector0_z);
  trackster_tree_->Branch("sigmaPCA1", &trackster_sigmaPCA1);
  trackster_tree_->Branch("sigmaPCA2", &trackster_sigmaPCA2);
  trackster_tree_->Branch("sigmaPCA3", &trackster_sigmaPCA3);
  trackster_tree_->Branch("id_probabilities", &trackster_id_probabilities);
  trackster_tree_->Branch("vertices_indexes", &trackster_vertices_indexes);
  trackster_tree_->Branch("vertices_x", &trackster_vertices_x);
  trackster_tree_->Branch("vertices_y", &trackster_vertices_y);
  trackster_tree_->Branch("vertices_z", &trackster_vertices_z);
  trackster_tree_->Branch("vertices_energy", &trackster_vertices_energy);
  trackster_tree_->Branch("vertices_correctedEnergy", &trackster_vertices_correctedEnergy);
  trackster_tree_->Branch("vertices_correctedEnergyUncertainty", &trackster_vertices_correctedEnergyUncertainty);
  trackster_tree_->Branch("vertices_multiplicity", &trackster_vertices_multiplicity); //NEW

  graph_tree_->Branch("linked_inners", &node_linked_inners);
  graph_tree_->Branch("linked_outers", &node_linked_outers);
  graph_tree_->Branch("isRootTrackster", &isRootTrackster);

  candidate_tree_->Branch("NCandidates", &nCandidates);
  candidate_tree_->Branch("candidate_charge", &candidate_charge);
  candidate_tree_->Branch("candidate_pdgId", &candidate_pdgId);
  candidate_tree_->Branch("candidate_id_probabilities", &candidate_id_probabilities);
  candidate_tree_->Branch("candidate_time", &candidate_time);
  candidate_tree_->Branch("candidate_timeErr", &candidate_time_err);
  candidate_tree_->Branch("candidate_energy", &candidate_energy);
  candidate_tree_->Branch("candidate_px", &candidate_px);
  candidate_tree_->Branch("candidate_py", &candidate_py);
  candidate_tree_->Branch("candidate_pz", &candidate_pz);
  candidate_tree_->Branch("track_in_candidate", &track_in_candidate);
  candidate_tree_->Branch("tracksters_in_candidate", &tracksters_in_candidate);

  tracksters_merged_tree_->Branch("NTrackstersMerged", &nTrackstersMerged);
  tracksters_merged_tree_->Branch("barycenter_x", &tracksters_merged_barycenter_x);
  tracksters_merged_tree_->Branch("barycenter_y", &tracksters_merged_barycenter_y);
  tracksters_merged_tree_->Branch("barycenter_z", &tracksters_merged_barycenter_z);
  tracksters_merged_tree_->Branch("EV1", &tracksters_merged_EV1);
  tracksters_merged_tree_->Branch("EV2", &tracksters_merged_EV2);
  tracksters_merged_tree_->Branch("EV3", &tracksters_merged_EV3);
  tracksters_merged_tree_->Branch("eVector0_x", &tracksters_merged_eVector0_x);
  tracksters_merged_tree_->Branch("eVector0_y", &tracksters_merged_eVector0_y);
  tracksters_merged_tree_->Branch("eVector0_z", &tracksters_merged_eVector0_z);
  tracksters_merged_tree_->Branch("sigmaPCA1", &tracksters_merged_sigmaPCA1);
  tracksters_merged_tree_->Branch("sigmaPCA2", &tracksters_merged_sigmaPCA2);
  tracksters_merged_tree_->Branch("sigmaPCA3", &tracksters_merged_sigmaPCA3);
  tracksters_merged_tree_->Branch("id_probabilities", &tracksters_merged_id_probabilities);
  
  cluster_tree_->Branch("seedID", &cluster_seedID);
  cluster_tree_->Branch("energy", &cluster_energy);
  cluster_tree_->Branch("correctedEnergy", &cluster_correctedEnergy);
  cluster_tree_->Branch("correctedEnergyUncertainty", &cluster_correctedEnergyUncertainty);
  cluster_tree_->Branch("position_x", &cluster_position_x);
  cluster_tree_->Branch("position_y", &cluster_position_y);
  cluster_tree_->Branch("position_z", &cluster_position_z);
  cluster_tree_->Branch("position_eta", &cluster_position_eta);
  cluster_tree_->Branch("position_phi", &cluster_position_phi);
  cluster_tree_->Branch("cluster_number_of_hits", &cluster_number_of_hits);
  cluster_tree_->Branch("cluster_hitsAndFractions_ID", &cluster_hitsAndFractions_ID);
  cluster_tree_->Branch("cluster_hitsAndFractions_value", &cluster_hitsAndFractions_value);
  event_index = 0;
}

void Ntupler::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  event_index++;
  clearVariables();

  //get all the tracksters
  edm::Handle<std::vector<ticl::Trackster>> tracksters_handle;
  event.getByToken(tracksters_token_, tracksters_handle);
  const auto& tracksters = *tracksters_handle;

  //get all the layer clusters
  edm::Handle<std::vector<reco::CaloCluster>> layer_clusters_h;
  event.getByToken(layer_clusters_token_, layer_clusters_h);
  const auto& clusters = *layer_clusters_h;

  //TICL Graph
  edm::Handle<TICLGraph> ticl_graph_h;
  event.getByToken(ticl_graph_token_, ticl_graph_h);
  const auto& graph = *ticl_graph_h;

  //TICL Candidate
  edm::Handle<std::vector<TICLCandidate>> candidates_h;
  event.getByToken(ticl_candidates_token_, candidates_h);
  const auto& ticlcandidates = *candidates_h;

  //Track
  edm::Handle<std::vector<reco::Track>> tracks_h;
  event.getByToken(tracks_token_, tracks_h);
  const auto& tracks = *tracks_h;

  //Tracksters merged
  edm::Handle<std::vector<ticl::Trackster>> tracksters_merged_h;
  event.getByToken(tracksters_merged_token_, tracksters_merged_h);
  const auto& trackstersmerged = *tracksters_merged_h;

  ev_event_ = event_index;
  ntracksters_ = tracksters.size();
  nclusters_ = clusters.size();

  for (auto trackster_iterator = tracksters.begin(); trackster_iterator != tracksters.end(); ++trackster_iterator) {
      //per-trackster analysis
    trackster_time.push_back(trackster_iterator->time());
    trackster_timeError.push_back(trackster_iterator->timeError());
    trackster_regressed_energy.push_back(trackster_iterator->regressed_energy());
    trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    trackster_raw_em_energy.push_back(trackster_iterator->raw_em_energy());
    trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
    trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
    trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
    trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    trackster_eVector0_x.push_back((trackster_iterator->eigenvectors()[0]).x());
    trackster_eVector0_y.push_back((trackster_iterator->eigenvectors()[0]).y());
    trackster_eVector0_z.push_back((trackster_iterator->eigenvectors()[0]).z());
    trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    std::vector<float_t> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    trackster_id_probabilities.push_back(id_probs);
    
    // Clusters
    std::vector<uint16_t> vertices_indexes;
    std::vector<float_t> vertices_x;
    std::vector<float_t> vertices_y;
    std::vector<float_t> vertices_z;
    std::vector<float_t> vertices_energy;
    std::vector<float_t> vertices_correctedEnergy;
    std::vector<float_t> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
        vertices_indexes.push_back(idx);
        auto associated_cluster = (*layer_clusters_h)[idx];
        vertices_x.push_back(associated_cluster.x());
        vertices_y.push_back(associated_cluster.y());
        vertices_z.push_back(associated_cluster.z());
        vertices_energy.push_back(associated_cluster.energy());
        vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
        vertices_correctedEnergyUncertainty.push_back(associated_cluster.correctedEnergyUncertainty());
    }
    trackster_vertices_indexes.push_back(vertices_indexes);
    trackster_vertices_x.push_back(vertices_x);
    trackster_vertices_y.push_back(vertices_y);
    trackster_vertices_z.push_back(vertices_z);
    trackster_vertices_energy.push_back(vertices_energy);
    trackster_vertices_correctedEnergy.push_back(vertices_correctedEnergy);
    trackster_vertices_correctedEnergyUncertainty.push_back(vertices_correctedEnergyUncertainty);
    
    // Multiplicity
    std::vector<float_t> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    trackster_vertices_multiplicity.push_back(vertices_multiplicity);
  }

  node_linked_inners.resize(tracksters.size());
  node_linked_outers.resize(tracksters.size());
  isRootTrackster.resize(tracksters.size(),false);
  for (size_t i = 0; i < tracksters.size(); ++i) {
    const auto &node = graph.getNode((int) i);
    node_linked_inners[i].insert(node_linked_inners[i].end(), node.getInner().begin(), node.getInner().end());
    node_linked_outers[i].insert(node_linked_outers[i].end(), node.getOuter().begin(), node.getOuter().end());
    
    if (node.getInner().empty()) isRootTrackster[i] = true;
  }

  for (auto cluster_iterator = clusters.begin(); cluster_iterator != clusters.end(); ++cluster_iterator) {
    auto lc_seed = cluster_iterator->seed();
    cluster_seedID.push_back(lc_seed);
    cluster_energy.push_back(cluster_iterator->energy());
    cluster_correctedEnergy.push_back(cluster_iterator->correctedEnergy());
    cluster_correctedEnergyUncertainty.push_back(cluster_iterator->correctedEnergyUncertainty());
    cluster_position_x.push_back(cluster_iterator->x());
    cluster_position_y.push_back(cluster_iterator->y());
    cluster_position_z.push_back(cluster_iterator->z());
    cluster_position_eta.push_back(cluster_iterator->eta());
    cluster_position_phi.push_back(cluster_iterator->phi());
    uint32_t number_of_hits = cluster_iterator->hitsAndFractions().size();
    cluster_number_of_hits.push_back(number_of_hits);    
    
    // auto thickness = rhtools_.getSiThickness(lc_seed);
    // auto thicknessIndex = rhtools_.getSiThickIndex(lc_seed);
    // std::cout << "THICKNESS " << thickness <<  " THICKNESS INDEX " << thicknessIndex << std::endl;
  }

  tracksters_in_candidate.resize(ticlcandidates.size());
  track_in_candidate.resize(ticlcandidates.size());
  nCandidates = ticlcandidates.size();
  for (size_t i = 0; i < ticlcandidates.size(); ++i) {
    const auto& candidate = ticlcandidates[i];
    candidate_charge.push_back(candidate.charge());
    candidate_pdgId.push_back(candidate.pdgId());
    candidate_energy.push_back(candidate.energy());
    candidate_px.push_back(candidate.px());
    candidate_py.push_back(candidate.py());
    candidate_pz.push_back(candidate.pz());
    candidate_time.push_back(candidate.time());
    candidate_time_err.push_back(candidate.timeError());
    std::vector<float_t> id_probs;
    for (int i = 0; i < 8; i++) {
      ticl::Trackster::ParticleType type = static_cast<ticl::Trackster::ParticleType>(i);
      id_probs.push_back(candidate.id_probability(type));
    }
    candidate_id_probabilities.push_back(id_probs);

    auto trackster_ptrs = candidate.tracksters();
    auto track_ptr = candidate.trackPtr();
    for (auto ts_ptr : trackster_ptrs) {
      uint16_t ts_idx = ts_ptr.get() - (edm::Ptr<ticl::Trackster>(tracksters_handle, 0)).get();
      tracksters_in_candidate[i].push_back(ts_idx);
    }
    
    if (track_ptr.isNull()) continue;
    uint16_t tk_idx = track_ptr.get() - (edm::Ptr<reco::Track>(tracks_h, 0)).get();
    track_in_candidate[i] = tk_idx;
  }

  nTrackstersMerged = trackstersmerged.size();
  for (size_t i = 0; i < trackstersmerged.size(); ++i) {
    const auto& tsm = trackstersmerged[i];
    tracksters_merged_barycenter_x.push_back(tsm.barycenter().x());
    tracksters_merged_barycenter_y.push_back(tsm.barycenter().y());
    tracksters_merged_barycenter_z.push_back(tsm.barycenter().z());
    tracksters_merged_EV1.push_back(tsm.eigenvalues()[0]);
    tracksters_merged_EV2.push_back(tsm.eigenvalues()[1]);
    tracksters_merged_EV3.push_back(tsm.eigenvalues()[2]);
    tracksters_merged_eVector0_x.push_back(tsm.eigenvectors(0).x());
    tracksters_merged_eVector0_y.push_back(tsm.eigenvectors(0).y());
    tracksters_merged_eVector0_z.push_back(tsm.eigenvectors(0).z());
    tracksters_merged_sigmaPCA1.push_back(tsm.sigmasPCA()[0]);
    tracksters_merged_sigmaPCA2.push_back(tsm.sigmasPCA()[1]);
    tracksters_merged_sigmaPCA3.push_back(tsm.sigmasPCA()[2]);
    
    std::vector<float_t> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(tsm.id_probabilities(i));
    tracksters_merged_id_probabilities.push_back(id_probs);
  }

  trackster_tree_->Fill();
  cluster_tree_->Fill();
  graph_tree_->Fill();
  candidate_tree_->Fill();
  tracksters_merged_tree_->Fill();
}

void Ntupler::endJob() {}

void Ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("ticlgraph", edm::InputTag("ticlGraph"));
  desc.add<edm::InputTag>("ticlcandidates", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Ntupler);
