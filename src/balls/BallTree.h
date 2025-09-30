/*
 * Copyright 2025 Pedro Maciel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#pragma once

#include <cassert>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "balls/Metric.h"
#include "balls/Vec3.h"


namespace balls {


// Utilities


struct Node {
    Vec3 center;      // center of the ball
    double radius{};  // radius under current metric
    int left{-1};
    int right{-1};
    uint32_t start{0};
    uint32_t end{0};  // leaf range in index array [start, end)
    bool is_leaf() const { return left < 0 && right < 0; }
};


// Ball Tree


class BallTree {
public:
    // Build constructor
    BallTree(std::vector<Vec3> pts, Metric* metric, uint32_t leaf_size = 16);

    // Load constructor (takes prebuilt arrays)
    BallTree(std::vector<Vec3> pts, std::vector<uint32_t> index, std::vector<Node> nodes, Metric* metric,
             uint32_t leaf_size = 16);

    // k nearest neighbors: returns vector of (distance, point index)
    std::vector<std::pair<double, uint32_t>> knn(const Vec3& q, uint32_t k) const;

    // radius query: return indices within distance <= r
    std::vector<uint32_t> radius_search(const Vec3& q, double r) const;

    // Serialization: save
    bool save(const std::string& path) const;

    // Serialization: load (caller supplies the id instance)
    static BallTree* load(const std::string& path, Metric::MetricId);

    const std::vector<Vec3>& points() const { return points_; }
    const std::vector<Node>& nodes() const { return nodes_; }

private:
    std::vector<Vec3> points_;
    std::vector<uint32_t> index_;  // permutation of indices into points_
    std::vector<Node> nodes_;
    Metric* metric_;  // non-owning; must outlive the BallTree
    uint32_t leaf_size_{16};

    void build();

    int build_recurse(uint32_t start, uint32_t end);

    uint32_t farthest_index(uint32_t pivot_idx, uint32_t start, uint32_t end) const;

    template <class Heap>
    void knn_recurse(int node_id, const Vec3& q, uint32_t k, Heap& best) const;

    void radius_recurse(int node_id, const Vec3& q, double r, std::vector<uint32_t>& out) const;
};


}  // namespace balls
