// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2021-2025, The OpenROAD Authors

#pragma once

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "odb/db.h"
#include "odb/dbTypes.h"
#include "odb/odb.h"

namespace odb {
class Rect;
class dbInst;
class dbModule;
class dbDatabase;
class dbITerm;
class dbTechLayer;
class dbBox;
class dbTrackGrid;
}  // namespace odb

namespace utl {
class Logger;
}

namespace par {

class Cluster;
using UniqueClusterVector = std::vector<std::unique_ptr<Cluster>>;

class Module;
using SharedModuleVector = std::vector<std::shared_ptr<Module>>;

class Module {
  public:
    //constructors
    Module(int i): id(i)
    { }

    //modifiers
    void setId(int i) { id = i; }
    void setAvgK(double aK) { avgK = aK; }
    void setAvgInsts(double aIn) { avgInsts = aIn; }
    void setAvgT(double aT) { avgT = aT; }
    void setInT(double inT) { avgInT = inT; }
    void setOutT(double outT) { avgOutT = outT; }
    void setSigmaT(double sigmaT) { SigmaT = sigmaT; }

    //accossers
    int getId() { return id; }
    double getAvgK() { return avgK; }
    double getAvgInsts() { return avgInsts; }
    double getAvgT() { return avgT; }
    double getInT() { return avgInT; }
    double getOutT() { return avgOutT; }
    double getSigmaT() { return SigmaT; }

  private:
    int id;
    double avgK = 0.0;
    double avgInsts = 0.0;
    double avgT = 0.0;
    double avgInT = 0.0;
    double avgOutT = 0.0;
    double SigmaT = 0.0;
};

class Cluster
{
  public:
    Cluster(const std::string& name) : id_(next_id_++), name_(name) {}

    int getId() const
    {
        return id_; 
    }

    // cluster name can be updated
    const std::string& getName() const
    {
        return name_;
    }
    
    const std::vector<odb::dbInst*> getLeafInsts() const
    {
        return leaf_insts_;
    }
    
    Cluster* getParent() const
    {
        return parent_;
    }

    const UniqueClusterVector& getChildren() const
    {
        return children_;
    }

    void setName(const std::string& name)
    {
        name_ = name; 
    }
    
    void setParent(Cluster* parent)
    {
        parent_ = parent;
    }

    void addLeafInst(odb::dbInst* leaf_std_cell)
    {
        leaf_insts_.push_back(leaf_std_cell);
    }

    int getNumInsts() const
    {
        return leaf_insts_.size();
    }
    void clearLeafInsts()
    {
        leaf_insts_.clear();
    }
    
    // Hierarchy Support
    void addChild(std::unique_ptr<Cluster> child)
    {
        children_.push_back(std::move(child));
    }
    std::unique_ptr<Cluster> releaseChild(const Cluster* candidate)
    {
        auto it = std::find_if(
        children_.begin(), children_.end(), [candidate](const auto& child) {
            return child.get() == candidate;
        });

        if (it != children_.end()) {
            std::unique_ptr<Cluster> released_child = std::move(*it);
            children_.erase(it);
            return released_child;
        }

        return nullptr;
    }
    void addChildren(UniqueClusterVector children)
    {
        std::move(children.begin(), children.end(), std::back_inserter(children_)); 
    }
    UniqueClusterVector releaseChildren()
    {
        UniqueClusterVector released_children = std::move(children_);
        children_.clear();
        
        return released_children;
    }

  private:
    int id_ = -1; // cluster id (a valid cluster id should be nonnegative)
    std::string name_;
    std::vector<odb::dbInst*> leaf_insts_;
    Cluster* parent_ = nullptr;  // parent of current cluster
    static int next_id_;
    UniqueClusterVector children_; 
};

class ClusterTree
{
  public:
    struct Maps {
        std::map<int, Cluster*> id_to_cluster;
        std::map<odb::dbInst*, int> inst_to_cluster_id;
        std::map<odb::dbBTerm*, int> bterm_to_cluster_id;
    };

    Maps maps;
    int large_net_threshold = 50; //same with RTLMP
    
    ClusterTree(odb::dbBlock* block) {
        root_ = std::make_unique<Cluster>("__cluster__");
        maps.id_to_cluster[root_->getId()] = root_.get();

        for (const auto inst : block->getInsts()) {
            root_->addLeafInst(inst);
            maps.inst_to_cluster_id[inst] = root_->getId();
        }
        for (const auto bterm : block->getBTerms()) {
            maps.bterm_to_cluster_id[bterm] = root_->getId();
        }
    }
    
    Cluster* getRoot() const { return root_.get(); }
    
    void registerCluster(Cluster* cluster) {
        maps.id_to_cluster[cluster->getId()] = cluster;
    }
    
    void addModule(std::shared_ptr<Module> module) 
    {
        modules_.push_back(module);
    }
    
    int getNumModules() const { return modules_.size(); }
    
    SharedModuleVector getModules() const { return modules_; }

private:
    std::unique_ptr<Cluster> root_; 
    SharedModuleVector modules_;
};
  

} //namespace
