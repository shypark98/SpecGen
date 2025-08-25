// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2019-2025, The OpenROAD Authors

#include "par/PartitionMgr.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_set>
#include <omp.h>

#include "TritonPart.h"
#include "ArtNetSpec.h"
#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"
#include "odb/db.h"
#include "sta/ArcDelayCalc.hh"
#include "sta/Bfs.hh"
#include "sta/ConcreteNetwork.hh"
#include "sta/Corner.hh"
#include "sta/DcalcAnalysisPt.hh"
#include "sta/ExceptionPath.hh"
#include "sta/FuncExpr.hh"
#include "sta/Graph.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/Liberty.hh"
#include "sta/MakeConcreteNetwork.hh"
#include "sta/Network.hh"
#include "sta/NetworkClass.hh"
#include "sta/ParseBus.hh"
#include "sta/PathAnalysisPt.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/PatternMatch.hh"
#include "sta/PortDirection.hh"
#include "sta/Sdc.hh"
#include "sta/Search.hh"
#include "sta/SearchPred.hh"
#include "sta/Sequential.hh"
#include "sta/Sta.hh"
#include "sta/Units.hh"
#include "sta/VerilogWriter.hh"
#include "utl/Logger.h"

using odb::dbBlock;
using odb::dbInst;
using odb::dbITerm;
using odb::dbBTerm;
using odb::dbMaster;
using odb::dbMasterType;
using odb::dbIntProperty;
using odb::dbIoType;

using sta::Cell;
using sta::CellPortBitIterator;
using sta::ConcreteNetwork;
using sta::dbNetwork;
using sta::Instance;
using sta::InstancePinIterator;
using sta::isBusName;
using sta::Library;
using sta::Net;
using sta::NetPinIterator;
using sta::NetTermIterator;
using sta::NetworkReader;
using sta::parseBusName;
using sta::Pin;
using sta::Port;
using sta::PortDirection;
using sta::Term;
using sta::writeVerilog;
using utl::PAR;

namespace par {

int Cluster::next_id_ = 0;

// find the correct brackets used in the liberty libraries.
static void determineLibraryBrackets(const dbNetwork* db_network,
                                     char* left,
                                     char* right)
{
  *left = '[';
  *right = ']';

  sta::LibertyLibraryIterator* lib_iter = db_network->libertyLibraryIterator();
  while (lib_iter->hasNext()) {
    const sta::LibertyLibrary* lib = lib_iter->next();
    *left = lib->busBrktLeft();
    *right = lib->busBrktRight();
  }
  delete lib_iter;
}

// determine the required direction of a port.
static PortDirection* determinePortDirection(const Net* net,
                                             const std::set<Instance*>* insts,
                                             const dbNetwork* db_network)
{
  bool local_only = true;
  bool locally_driven = false;
  bool externally_driven = false;

  NetTermIterator* term_iter = db_network->termIterator(net);
  while (term_iter->hasNext()) {
    Term* term = term_iter->next();
    PortDirection* dir = db_network->direction(db_network->pin(term));
    if (dir->isAnyInput()) {
      externally_driven = true;
    }
    local_only = false;
  }
  delete term_iter;

  if (insts != nullptr) {
    NetPinIterator* pin_iter = db_network->pinIterator(net);
    while (pin_iter->hasNext()) {
      const Pin* pin = pin_iter->next();
      const PortDirection* dir = db_network->direction(pin);
      Instance* inst = db_network->instance(pin);

      if (insts->find(inst) == insts->end()) {
        local_only = false;
        if (dir->isAnyOutput()) {
          externally_driven = true;
        }
      } else {
        if (dir->isAnyOutput()) {
          locally_driven = true;
        }
      }
    }
    delete pin_iter;
  }

  // no port is needed
  if (local_only) {
    return nullptr;
  }

  if (locally_driven && externally_driven) {
    return PortDirection::bidirect();
  }

  if (externally_driven) {
    return PortDirection::input();
  }

  return PortDirection::output();
}

void PartitionMgr::writeArtNetSpec(const char* fileName) {
    
    std::unordered_map<std::string, std::pair<int, bool>> onlyUseMasters;
    std::string top_name;
    int numInsts = 0;
    int numPIs = 0;
    int numPOs = 0;
    int numSeq = 0;
    int Dmax = -1;
    int MDmax = -1;
    float Rratio;
    float p;
    float q;

    getFromODB(onlyUseMasters, top_name, numInsts, numPIs, numPOs, numSeq);
    std::cout << "getFromODB done" << std::endl;
    getFromSTA(Dmax, MDmax);
    std::cout << "getFromSTA done" << std::endl;
    getFromPAR(Rratio, p, q);
    std::cout << "getFromPAR done" << std::endl;
    std::cout << "Rratio: " << Rratio << " " << "rentP " << p << " " << "rentQ " << q << std::endl;
    writeFile(onlyUseMasters, top_name, numInsts, numPIs, numPOs, 
              numSeq, Dmax, MDmax, Rratio, p, q, fileName);

}

void PartitionMgr::getFromODB(std::unordered_map<std::string, std::pair<int, bool>>& onlyUseMasters,
                              std::string& top_name,
                              int& numInsts,
                              int& numPIs,
                              int& numPOs,
                              int& numSeq) 
{
    auto block = getDbBlock(); 
    odb::dbSet<dbInst> insts = block->getInsts();
    odb::dbSet<dbBTerm> bterms = block->getBTerms();
    numInsts = insts.size();

    for(auto bterm_itr = bterms.begin(); bterm_itr != bterms.end(); ++bterm_itr) {
        dbBTerm* bterm = *bterm_itr;
        if(bterm->getIoType() == odb::dbIoType::INPUT)
            numPIs++;
        if(bterm->getIoType() == odb::dbIoType::OUTPUT)
            numPOs++;
    }

    for(auto inst_itr = insts.begin(); inst_itr != insts.end(); ++inst_itr) {
        dbInst* inst = *inst_itr;
        dbMaster* master = inst->getMaster();
        bool isMacro = (master->getType() == dbMasterType::BLOCK ? 1 : 0);
        if (master->isSequential())
            numSeq++;
        if(onlyUseMasters.find(master->getName()) == onlyUseMasters.end())
            onlyUseMasters[master->getName()] = std::make_pair(0, isMacro);
        onlyUseMasters[master->getName()].first++;
    }
}

void PartitionMgr::getFromSTA(int& Dmax,
                              int& MDmax) 
{
    BuildTimingPath(Dmax, MDmax);
}

void PartitionMgr::BuildTimingPath(int& Dmax,
                                   int& MDmax) 
{
  sta_->ensureGraph();     // Ensure that the timing graph has been built
  sta_->searchPreamble();  // Make graph and find delays
  sta_->ensureLevelized();
  // Step 1:  find the top_n critical timing paths
  sta::ExceptionFrom* e_from = nullptr;
  sta::ExceptionThruSeq* e_thrus = nullptr;
  sta::ExceptionTo* e_to = nullptr;
  bool include_unconstrained = false;
  bool get_max = true;  // max for setup check, min for hold check
  // Timing paths are grouped into path groups according to the clock
  // associated with the endpoint of the path, for example, path group for clk
  //int group_count = top_n_;
  int group_count = 1000;
  int endpoint_count = 1;  // The number of paths to report for each endpoint.
  // Definition for findPathEnds function in Search.hh
  // PathEndSeq *findPathEnds(ExceptionFrom *from,
  //              ExceptionThruSeq *thrus,
  //              ExceptionTo *to,
  //              bool unconstrained,
  //              const Corner *corner,
  //              const MinMaxAll *min_max,
  //              int group_count,
  //              int endpoint_count,
  //              bool unique_pins,
  //              float slack_min,
  //              float slack_max,
  //              bool sort_by_slack,
  //              PathGroupNameSet *group_names,
  //              bool setup,
  //              bool hold,
  //              bool recovery,
  //              bool removal,
  //              bool clk_gating_setup,
  //              bool clk_gating_hold);
  // PathEnds represent search endpoints that are either unconstrained or
  // constrained by a timing check, output delay, data check, or path delay.
  sta::PathEndSeq path_ends = sta_->search()->findPathEnds(  // from, thrus, to,
                                                             // unconstrained
      e_from,   // return paths from a list of clocks/instances/ports/register
                // clock pins or latch data pins
      e_thrus,  // return paths through a list of instances/ports/nets
      e_to,     // return paths to a list of clocks/instances/ports or pins
      include_unconstrained,  // return unconstrained paths
      // corner, min_max,
      sta_->cmdCorner(),  // return paths for a process corner
      get_max ? sta::MinMaxAll::max()
              : sta::MinMaxAll::min(),  // return max/min paths checks
      // group_count, endpoint_count, unique_pins
      group_count,     // number of paths in total
      endpoint_count,  // number of paths for each endpoint
      true,
      -sta::INF,
      sta::INF,  // slack_min, slack_max,
      true,      // sort_by_slack
      nullptr,   // group_names
      // setup, hold, recovery, removal,
      get_max,
      !get_max,
      false,
      false,
      // clk_gating_setup, clk_gating_hold
      false,
      false);
  
  auto block = getDbBlock();  
  std::map<std::string, int> pathDepthMap;
    
  // check all the timing paths
  for (auto& path_end : path_ends) {
    // Printing timing paths to logger
    // sta_->reportPathEnd(path_end);
    auto* path = path_end->path();
    
    int depth = 0;
    std::string endPointName;
    std::unordered_set<std::string> visitedInstances;
    std::unordered_set<std::string> visitedBterms;

    sta::PathExpanded expand(path, sta_);
    for (size_t i = 0; i < expand.size(); i++) {
      const sta::Path* ref = expand.path(i);
      sta::Pin* pin = ref->vertex(sta_)->pin();
      // Nets connect pins at a level of the hierarchy
      auto net = db_network_->net(pin);  // sta::Net*
      // Check if the pin is connected to a net
      if (net == nullptr)
        continue;  // check if the net exists
        
      std::string name;

      if (db_network_->isTopLevelPort(pin) == true) {
        auto bterm = block->findBTerm(db_network_->pathName(pin));
        name = bterm->getName();
        if (visitedBterms.insert(name).second)
            depth++;

      } else {
        auto inst = db_network_->instance(pin);
        auto db_inst = block->findInst(db_network_->pathName(inst));
        name = db_inst->getName();
        if (visitedInstances.insert(name).second)
            depth++;
      }
    
      if (i == expand.size() - 1) {
        endPointName = name;
        pathDepthMap[endPointName] = depth;
      }
    }
  } // path_end
  
  int ff_max = 0;
  int mac_max = 0;
  for (auto path : pathDepthMap) {
     auto inst = block->findInst((path.first).c_str());
     if (inst) {
        if (inst->getMaster()->isBlock()) {
          mac_max = std::max(path.second, mac_max);   
        } else { 
          ff_max = std::max(path.second, ff_max);   
        }
     }
  }
  Dmax = ff_max;
  MDmax = mac_max;
}

void PartitionMgr::getFromPAR(float& Rratio,
                              float& p,
                              float& q) 
{
    getRents(Rratio, p, q);
}

void PartitionMgr::getRents(float& Rratio,
                            float& p,
                            float& q)
{
    auto block = getDbBlock();
    auto tree = ClusterTree(block); 
    if (!tree.getRoot()) {
        Rratio = 0.0f;
        p = 0.0f;
        q = 0.0f;
        return;
    }

    std::vector<Cluster*> current_level_clusters;
    current_level_clusters.push_back(tree.getRoot());
    
    double totPins = 0;
    for (dbInst* inst : block->getInsts()) {
        for (dbITerm* inst_iterm : inst->getITerms()) {
            if(inst_iterm->getIoType() == dbIoType::INPUT ||
               inst_iterm->getIoType() == dbIoType::OUTPUT)
                totPins++;
        }
    }
    double avgK = totPins / block->getInsts().size(); 
    
    double numPIs = 0, numPOs = 0;
    auto bterms = block->getBTerms();
    for(auto bterm_itr = bterms.begin(); bterm_itr != bterms.end(); ++bterm_itr) {
        dbBTerm* bterm = *bterm_itr;
        if(bterm->getIoType() == odb::dbIoType::INPUT)
            numPIs++;
        if(bterm->getIoType() == odb::dbIoType::OUTPUT)
            numPOs++;
    }
    
    auto m = std::make_shared<Module>(tree.getNumModules());
    m->setAvgK(avgK);
    m->setAvgInsts(block->getInsts().size());
    m->setAvgT(block->getBTerms().size());
    m->setInT(numPIs);
    m->setOutT(numPOs);
    m->setSigmaT(0.0);
    tree.addModule(m);
    
    int level = 0;
    while (true) {
        bool earlyStop = false;
        std::vector<Cluster*> next_level_clusters;
        std::cout << "start partitioning level " << level << std::endl; 
        #pragma omp parallel for shared(current_level_clusters, next_level_clusters, earlyStop)
        for (Cluster* current_cluster : current_level_clusters) {
            auto [child_0, child_1] = partitionCluster(current_cluster, tree);
            
            if (child_1 != nullptr) {
                #pragma omp critical
                {
                    next_level_clusters.push_back(child_0);
                    next_level_clusters.push_back(child_1);
                }
            } else { // child_1 == nullptr
                #pragma omp atomic
                earlyStop = true;
            }
            if (earlyStop)
                break;
        }
        //if (!earlyStop) {
            current_level_clusters = next_level_clusters;
            int sampleNum = current_level_clusters.size(); 
            std::vector<double> numInsts_vec, sumT_vec, inT_vec, outT_vec;
        
            for (auto cluster : current_level_clusters) {
                double sumT = 0, inT = 0, outT = 0; 
                getClusterIONum(cluster, sumT, inT, outT);
                numInsts_vec.push_back(cluster->getNumInsts());
                sumT_vec.push_back(sumT);
                inT_vec.push_back(inT);
                outT_vec.push_back(outT);
            }
            
            double total_insts = std::accumulate(numInsts_vec.begin(), numInsts_vec.end(), 0.0);
            double total_sumT = std::accumulate(sumT_vec.begin(), sumT_vec.end(), 0.0);
            double total_inT = std::accumulate(inT_vec.begin(), inT_vec.end(), 0.0);
            double total_outT = std::accumulate(outT_vec.begin(), outT_vec.end(), 0.0);
        
            double avgInsts = total_insts / sampleNum;
            double avgT = total_sumT / sampleNum;
            double avgInT = total_inT / sampleNum;
            double avgOutT = total_outT / sampleNum;
            
            double varT = 0.0;
            if (sampleNum > 1) {
                double sqdiffs = 0.0;
                for (double val : sumT_vec) {
                    sqdiffs += (val - avgT) * (val - avgT);
                }
                varT = sqdiffs / (sampleNum - 1);
            }
            double stdevT = sqrt(varT);
        
        
            if (avgInsts >= 1 && avgT >= 1) {
                // in this situation, module means partition level
                auto m = std::make_shared<Module>(tree.getNumModules());
                m->setAvgInsts(avgInsts);
                m->setAvgT(avgT);
                m->setInT(avgInT);
                m->setOutT(avgOutT);
                m->setSigmaT(stdevT);
                tree.addModule(m);
            }
        //}
        if (earlyStop) {
            break;
        }
        level++;
    }
    linCurvFit(tree, Rratio, p, q);
}

std::pair<Cluster*, Cluster*>
PartitionMgr::partitionCluster(Cluster* parent,
                               ClusterTree& tree) 
{   
    // Recursion termination condition 
    if (parent->getLeafInsts().size() <= 5) {
        std::cout << "Cluster '" << parent->getName() << "' is small enough ("
                  << parent->getLeafInsts().size() << " insts). Stopping recursion.\n";

        return {parent, nullptr};
    }
    
    auto block = getDbBlock();
    std::map<int, int> cluster_vertex_id_map;
    std::vector<float> vertex_weight;
    int vertex_id = 0;

    // Register external clusters as vertices
    for (auto& [cluster_id, cluster] : tree.maps.id_to_cluster) {
        cluster_vertex_id_map[cluster_id] = vertex_id++;
        vertex_weight.push_back(0.0f);
    }
    const int num_other_cluster_vertices = vertex_id;

    std::vector<odb::dbInst*> insts;
    std::map<odb::dbInst*, int> inst_vertex_id_map;
    
    for (auto& std_cell : parent->getLeafInsts()) {
        inst_vertex_id_map[std_cell] = vertex_id++;
        vertex_weight.push_back(computeMicronArea(std_cell));
        insts.push_back(std_cell);
    }

    std::vector<std::vector<int>> hyperedges;
    for (odb::dbNet* net : block->getNets()) {
        if (net->getSigType().isSupply()) {
            continue;
        }

        int driver_id = -1;
        std::set<int> loads_id;
        bool ignore = false;
        for (odb::dbITerm* iterm : net->getITerms()) {
            odb::dbInst* inst = iterm->getInst();
            if (isIgnoredInst(inst)) {
                ignore = true;
                break;
            }

            const int cluster_id = tree.maps.inst_to_cluster_id.at(inst);
            int vertex_id = (cluster_id != parent->getId())
                            ? cluster_vertex_id_map[cluster_id]
                            : inst_vertex_id_map[inst];
            if (iterm->getIoType() == odb::dbIoType::OUTPUT) {
                driver_id = vertex_id;
            } else {
                loads_id.insert(vertex_id);
            }
        }

        if (ignore) {
            continue;
        }

        for (odb::dbBTerm* bterm : net->getBTerms()) {
            const int cluster_id = tree.maps.bterm_to_cluster_id.at(bterm);
            if (bterm->getIoType() == odb::dbIoType::INPUT) {
                driver_id = cluster_vertex_id_map[cluster_id];
            } else {
                loads_id.insert(cluster_vertex_id_map[cluster_id]);
            }
        }
        loads_id.insert(driver_id);
        if (driver_id != -1 && loads_id.size() > 1
            && loads_id.size() < tree.large_net_threshold) {
            std::vector<int> hyperedge;
            hyperedge.insert(hyperedge.end(), loads_id.begin(), loads_id.end());
            hyperedges.push_back(hyperedge);
        }
    }
    const int seed = 0;
    constexpr float default_balance_constraint = 1.0f;
    float balance_constraint = default_balance_constraint;
    const int num_parts = 2;  // We use two-way partitioning here
    const int num_vertices = static_cast<int>(vertex_weight.size());
    std::vector<float> hyperedge_weights(hyperedges.size(), 1.0f);
    
    // Due to the discrepancy that may exist between the weight of vertices
    // that represent macros/std cells, the partitioner may fail to meet the
    // balance constraint. This may cause the output to be completely unbalanced
    // and lead to infinite partitioning recursion. To handle that, we relax
    // the constraint until we find a reasonable split.
    constexpr float balance_constraint_relaxation_factor = 10.0f;
    std::vector<int> solution;
    do {
        if (balance_constraint >= 90) {
            std::cout << "Cannot find a balanced partitioning for the clusters." << std::endl;
        }

        solution = PartitionKWaySimpleMode(num_parts,
                                           balance_constraint,
                                           seed,
                                           hyperedges,
                                           vertex_weight,
                                           hyperedge_weights);

        balance_constraint += balance_constraint_relaxation_factor;
    } while (partitionerSolutionIsFullyUnbalanced(solution, 
                                                  num_other_cluster_vertices));
    
    parent->clearLeafInsts();

    const std::string cluster_name = parent->getName();
    parent->setName(cluster_name + std::string("_0"));
    auto cluster_part_1 = std::make_unique<Cluster>(
            cluster_name + std::string("_1"));

    for (int i = num_other_cluster_vertices; i < num_vertices; i++) {
        odb::dbInst* inst = insts[i - num_other_cluster_vertices];
        if (solution[i] == 0) {
            parent->addLeafInst(inst);
        } else {
            cluster_part_1->addLeafInst(inst);
        }
    }

    Cluster* raw_part_1 = cluster_part_1.get();
    cluster_part_1->setParent(parent);
    parent->addChild(std::move(cluster_part_1));
    tree.registerCluster(raw_part_1);
    // Recursive break the cluster
    // until the size of the cluster is less than max_num_inst_
    return {parent, raw_part_1};
}

bool PartitionMgr::partitionerSolutionIsFullyUnbalanced(
    const std::vector<int>& solution,
    const int num_other_cluster_vertices)
{
  // The partition of the first vertex which represents
  // an actual macro or std cell.
  const int first_vertex_partition = solution[num_other_cluster_vertices];
  const int number_of_vertices = static_cast<int>(solution.size());

  // Skip all the vertices that represent other clusters.
  for (int vertex_id = num_other_cluster_vertices;
       vertex_id < number_of_vertices;
       ++vertex_id) {
    if (solution[vertex_id] != first_vertex_partition) {
      return false;
    }
  }

  return true;
}

float PartitionMgr::computeMicronArea(odb::dbInst* inst)
{
  auto block = getDbBlock();
  const float width = static_cast<float>(
      block->dbuToMicrons(inst->getBBox()->getBox().dx()));
  const float height = static_cast<float>(
      block->dbuToMicrons(inst->getBBox()->getBox().dy()));

  return width * height;
}

bool PartitionMgr::isIgnoredInst(odb::dbInst* inst)
{
  odb::dbMaster* master = inst->getMaster();

  return master->isPad() || master->isCover() || master->isEndCap()
         || inst->isFixed();
}

void PartitionMgr::getClusterIONum(Cluster* parent,
                                   double &sumT,
                                   double &numPIs,
                                   double &numPOs) 
{

    // Build submodule partitions
    std::map<sta::Net*, sta::Port*> sta_port_map;
    std::set<Instance*> instance_set;
    for (odb::dbInst* inst : parent->getLeafInsts())
        instance_set.insert(db_network_->dbToSta(inst));

    // Build submodule partitions
    const std::string subModule_name = parent->getName();

    // Create new network and library
    sta::NetworkReader* sub_network = sta::makeConcreteNetwork();
    sta::Library* sub_library = sub_network->makeLibrary("Partitions", nullptr);

    sta::Instance* sub_inst = buildSubInst(subModule_name.c_str(),
                                           "sub",
                                           sub_library,
                                           sub_network,
                                           nullptr, // No top instance
                                           &instance_set,
                                           &sta_port_map);

    if (!sub_inst)
        logger_->report("Failed to create instance ", subModule_name);

    reinterpret_cast<sta::ConcreteNetwork*>(sub_network)->setTopInstance(sub_inst);
    
    auto network_ = sta_->network();
    for (const auto& [net, port] : sta_port_map) {
        sta::PortDirection* dir = sub_network->direction(port);
        if (!sub_network->direction(port)->isPowerGround()) {
            if (dir->isAnyInput()) {
                if (sub_network->isBus(port))
                    std::cout << network_->fromIndex(port) << " " << network_->toIndex(port)
                        << " " << sub_network->size(port) << std::endl;                
                numPIs += sub_network->size(port);
            } else {
                if (sub_network->isBus(port))
                    std::cout << network_->fromIndex(port) << " " << network_->toIndex(port)
                        << " " << sub_network->size(port) << std::endl;
                numPOs += sub_network->size(port);
            }
        }
    }
    sumT = numPIs + numPOs;
    delete sub_network; // Clean up memory

}

sta::Instance* PartitionMgr::buildSubInst(const char* name,
                                     const char* port_prefix,
                                     sta::Library* library,
                                     sta::NetworkReader* network,
                                     sta::Instance* parent,
                                     const std::set<Instance*>* insts,
                                     std::map<Net*, Port*>* port_map)
{
  // build cell
  Cell* cell = network->makeCell(library, name, false, nullptr);

  // add global ports
  auto pin_iter = db_network_->pinIterator(db_network_->topInstance());
  while (pin_iter->hasNext()) {
    const Pin* pin = pin_iter->next();

    bool add_port = false;
    Net* net = db_network_->net(db_network_->term(pin));
    if (net) {
      NetPinIterator* net_pin_iter = db_network_->pinIterator(net);
      while (net_pin_iter->hasNext()) {
        // check if port is connected to instance in this partition
        if (insts->find(db_network_->instance(net_pin_iter->next()))
            != insts->end()) {
          add_port = true;
          break;
        }
      }
      delete net_pin_iter;
    }

    if (add_port) {
      const char* portname = db_network_->name(pin);

      Port* port = network->makePort(cell, portname);
      // copy exactly the parent port direction
      network->setDirection(port, db_network_->direction(pin));
      PortDirection* sub_module_dir
          = determinePortDirection(net, insts, db_network_);
      if (sub_module_dir != nullptr) {
        network->setDirection(port, sub_module_dir);
      }
      port_map->insert({net, port});
    }
  }
  delete pin_iter;

  // make internal ports for partitions and if port is not needed.
  std::set<Net*> local_nets;
  for (Instance* inst : *insts) {
    InstancePinIterator* pin_iter = db_network_->pinIterator(inst);
    while (pin_iter->hasNext()) {
      Net* net = db_network_->net(pin_iter->next());
      if (net != nullptr &&                          // connected
          port_map->find(net) == port_map->end() &&  // port not present
          local_nets.find(net) == local_nets.end()) {
        // check if connected to anything in a different partition
        NetPinIterator* net_pin_iter = db_network_->pinIterator(net);
        while (net_pin_iter->hasNext()) {
          Net* net = db_network_->net(net_pin_iter->next());
          PortDirection* port_dir
              = determinePortDirection(net, insts, db_network_);
          if (port_dir == nullptr) {
            local_nets.insert(net);
            continue;
          }
          std::string port_name = port_prefix;
          port_name += db_network_->name(net);

          Port* port = network->makePort(cell, port_name.c_str());
          network->setDirection(port, port_dir);

          port_map->insert({net, port});
          break;
        }
        delete net_pin_iter;
      }
    }
    delete pin_iter;
  }
  // loop over buses and to ensure all bit ports are created, only needed for
  // partitioned modules
  char path_escape = db_network_->pathEscape();
  char left_bracket;
  char right_bracket;
  determineLibraryBrackets(db_network_, &left_bracket, &right_bracket);
  std::map<std::string, std::vector<Port*>> port_buses;
  for (auto& [net, port] : *port_map) {
      std::string portname = network->name(port);

    // check if bus and get name
    if (isBusName(portname.c_str(), left_bracket, right_bracket, path_escape)) {
      std::string bus_name;
      bool is_bus;
      int idx;
      parseBusName(portname.c_str(),
                   left_bracket,
                   right_bracket,
                   path_escape,
                   is_bus,
                   bus_name,
                   idx);

      portname = bus_name;
      port_buses[portname].push_back(port);
    }
  }
  for (auto& [bus, ports] : port_buses) {
    std::set<int> port_idx;
    std::set<PortDirection*> port_dirs;
    for (Port* port : ports) {
      std::string bus_name;
      bool is_bus;
      int idx;
      parseBusName(network->name(port),
                   left_bracket,
                   right_bracket,
                   path_escape,
                   is_bus,
                   bus_name,
                   idx);

      port_idx.insert(idx);
      port_dirs.insert(network->direction(port));
    }

    // determine real direction of port
    PortDirection* overall_direction = nullptr;
    if (port_dirs.size() == 1) {  // only one direction is used.
      overall_direction = *port_dirs.begin();
    } else {
      overall_direction = PortDirection::bidirect();
    }

    // set port direction to match
    for (Port* port : ports) {
      network->setDirection(port, overall_direction);
    }

    // fill in missing ports in bus
    const auto [min_idx, max_idx]
        = std::minmax_element(port_idx.begin(), port_idx.end());
    for (int idx = *min_idx; idx <= *max_idx; idx++) {
      if (port_idx.find(idx) == port_idx.end()) {
        // build missing port
        std::string portname = bus;
        portname += left_bracket + std::to_string(idx) + right_bracket;
        Port* port = network->makePort(cell, portname.c_str());
        network->setDirection(port, overall_direction);
      }
    }
  }

  network->groupBusPorts(cell, [](const char*) { return true; });

  // build instance
  std::string instname = name;
  instname += "_inst";
  Instance* inst = network->makeInstance(cell, instname.c_str(), parent);

  // create nets for ports in cell
  for (auto& [db_net, port] : *port_map) {
    Net* net = network->makeNet(network->name(port), inst);
    Pin* pin = network->makePin(inst, port, nullptr);
    network->makeTerm(pin, net);
  }

 // create and connect instances
  for (Instance* instance : *insts) {
    Instance* leaf_inst = network->makeInstance(
        db_network_->cell(instance), db_network_->name(instance), inst);

    InstancePinIterator* pin_iter = db_network_->pinIterator(instance);
    while (pin_iter->hasNext()) {
      Pin* pin = pin_iter->next();
      Net* net = db_network_->net(pin);
      if (net != nullptr) {  // connected
        Port* port = db_network_->port(pin);

        // check if connected to a port
        auto port_find = port_map->find(net);
        if (port_find != port_map->end()) {
          Net* new_net
              = network->findNet(inst, network->name(port_find->second));
          network->connect(leaf_inst, port, new_net);
        } else {
          Net* new_net = network->findNet(inst, db_network_->name(net));
          if (new_net == nullptr) {
            new_net = network->makeNet(db_network_->name(net), inst);
          }
          network->connect(leaf_inst, port, new_net);
        }
      }
    }
    delete pin_iter;
  }

  return inst;
}

//from RentCon
void PartitionMgr::linCurvFit(ClusterTree& tree,
                              float& Rratio,
                              float& p,
                              float& q)
{
    
    int n = tree.getNumModules();
    double *x = new double[n];
    double *y = new double[n];
    
    auto modules = tree.getModules();
    double b = log(modules[0]->getAvgK());
    for (int i = 0; i < n; i++)
    {
        auto m = modules[i];
        x[i] = log(m->getAvgInsts());
        y[i] = log(m->getAvgT()) - b;
    }
    
    for (int i = 0; i < n; i++)
    {
        auto m = modules[i];
        double numInsts = m->getAvgInsts();
        double I = m->getInT();
        double O = m->getOutT();
        double T = m->getAvgT();
        double sigT = m->getSigmaT();
        std::cout << "Level: " << i << " G: " << numInsts << " In: " << I << " Out: " << O << std::endl;
    }

    auto [ratio, rentP, std_dev] = fitRent(x, y, n);
    delete x;
    delete y;
    Rratio = ratio;
    p = rentP;
    q = std_dev;
}

//from RentCon
std::tuple<double, double, double> 
PartitionMgr::fitRent(double* x, double* y, int n)
{
    int minPntNum = (int)(n * 0.75);
    double bestRent;
    int totPoints = n;
    int bestN = n;
    double rentP, cov11, sumsq;
    
    fit_mul(x, 1, y, 1, n, &rentP, &cov11, &sumsq);
    bestRent = rentP;

    double oldDev = sqrt(sumsq / n);

    while (n > minPntNum)
    {
        n--;
        fit_mul(x, 1, y, 1, n, &rentP, &cov11, &sumsq);
        // compute the standard deviation of the residuals
        double newDev = sqrt(sumsq / n);
        if (newDev > oldDev)
            break;
        else
        {
            oldDev = newDev;
            bestN = n;
            bestRent = rentP;
        }
    }
    double Rratio;
    if (bestN == totPoints)
        Rratio = 0.9;
    else
        Rratio = double(bestN) / totPoints;

    return std::make_tuple(Rratio, bestRent, oldDev);
}

//from gsl library
void PartitionMgr::fit_mul(const double *x, 
                           const size_t xstride,
                           const double *y, 
                           const size_t ystride,
                           const size_t n,
                           double *c1, 
                           double *cov_11, 
                           double *sumsq)
{
  double m_x = 0, m_y = 0, m_dx2 = 0, m_dxdy = 0;
  size_t i;
  for (i = 0; i < n; i++)
    {
      m_x += (x[i * xstride] - m_x) / (i + 1.0);
      m_y += (y[i * ystride] - m_y) / (i + 1.0);
    }
  for (i = 0; i < n; i++)
    {
      const double dx = x[i * xstride] - m_x;
      const double dy = y[i * ystride] - m_y;
      m_dx2 += (dx * dx - m_dx2) / (i + 1.0);
      m_dxdy += (dx * dy - m_dxdy) / (i + 1.0);
    }
  /* In terms of y =  b x */
  {
    double s2 = 0, d2 = 0;
    double b = (m_x * m_y + m_dxdy) / (m_x * m_x + m_dx2);
    *c1 = b;
    /* Compute chi^2 = \sum (y_i -  b * x_i)^2 */
    for (i = 0; i < n; i++)
      {
        const double dx = x[i * xstride] - m_x;
        const double dy = y[i * ystride] - m_y;
        const double d = (m_y - b * m_x) + dy - b * dx;
        d2 += d * d;
      }
    s2 = d2 / (n - 1.0);        /* chisq per degree of freedom */
    *cov_11 = s2 * 1.0 / (n * (m_x * m_x + m_dx2));
    *sumsq = d2;
  }
}

void PartitionMgr::writeFile(std::unordered_map<std::string, std::pair<int, bool>>& onlyUseMasters,
                             std::string& top_name,
                             int& numInsts,
                             int& numPIs,
                             int& numPOs,
                             int& numSeq,
                             int& Dmax,
                             int& MDmax,
                             float& Rratio,
                             float& p,
                             float& q,
                             const char* fileName) 
{

    std::ofstream outFile(fileName);
    if (!outFile.good()) {
        std::cout << "Error: cannot open file " << fileName << std::endl;
        exit(0);
    }

    outFile << "LIBRARY" << std::endl;
    outFile << "NAME lib" << std::endl;

    // unordered_map<string, int> --> cellName / isMacro
    for (auto it : onlyUseMasters) {
        if (!it.second.second)
            outFile << "STD_CELL " << it.first  << std::endl;
        else
            outFile << "MACRO_CELL " << it.first  << std::endl;
    }
    outFile << std::endl;

    outFile << "CIRCUIT" << std::endl;
    outFile << "NAME " << top_name << std::endl;
    outFile << "LIBRARIES lib" << std::endl;
    outFile << "DISTRIBUTION ";
    for (auto it : onlyUseMasters) {
        outFile << it.second.first << " ";
    }
    outFile << std::endl;
    outFile << "SIZE " << int(numInsts * Rratio)  << std::endl;
    outFile << "p " << p << std::endl;
    outFile << "q " << q << std::endl;
    outFile << "END" << std::endl;
    outFile << "SIZE " << numInsts << std::endl;
    outFile << "I " << numPIs << std::endl;
    outFile << "O " << numPOs << std::endl;
    outFile << "END" << std::endl;
    outFile.close();
    outFile << std::endl;
    std::cout << "Dmax " << Dmax << " " << "MDmax " << MDmax << std::endl;
}

}  // namespace par
