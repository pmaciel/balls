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


#include "balls/BallTree.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>


namespace balls {


namespace {


inline std::streamsize to_streamsize(uint64_t v) {
    assert(v <= static_cast<uint64_t>(std::numeric_limits<std::streamsize>::max()));
    return static_cast<std::streamsize>(v);
}


struct Header {
    uint32_t magic{0xBA115EED};  // "BALLSEED"
    uint32_t version{1};
    Metric::MetricId metric_id{Metric::GreatCircle};
    uint64_t n_points{0};
    uint64_t n_nodes{0};
};


}  // namespace


BallTree::BallTree(std::vector<Vec3> pts, Metric* metric, uint32_t leaf_size) :
    points_(std::move(pts)), metric_(metric), leaf_size_(leaf_size) {
    if (metric_ == nullptr) {
        throw std::invalid_argument("metric pointer must not be null");
    }
    build();
}


BallTree::BallTree(std::vector<Vec3> pts, std::vector<uint32_t> index, std::vector<Node> nodes, Metric* metric,
                   uint32_t leaf_size) :
    points_(std::move(pts)),
    index_(std::move(index)),
    nodes_(std::move(nodes)),
    metric_(metric),
    leaf_size_(leaf_size) {
    if (metric_ == nullptr) {
        throw std::invalid_argument("metric pointer must not be null");
    }
}


std::vector<std::pair<double, uint32_t> > BallTree::knn(const Vec3& q, uint32_t k) const {
    if (points_.empty() || k == 0) {
        return {};
    }
    using Pair = std::pair<double, uint32_t>;
    auto cmp   = [](const Pair& a, const Pair& b) { return a.first < b.first; };  // max-heap
    std::priority_queue<Pair, std::vector<Pair>, decltype(cmp)> best(cmp);
    knn_recurse(0, q, k, best);
    std::vector<Pair> out;
    out.reserve(best.size());
    while (!best.empty()) {
        out.push_back(best.top());
        best.pop();
    }
    std::reverse(out.begin(), out.end());
    return out;
}


std::vector<uint32_t> BallTree::radius_search(const Vec3& q, double r) const {
    std::vector<uint32_t> out;
    radius_recurse(0, q, r, out);
    return out;
}


bool BallTree::save(const std::string& path) const {
    std::ofstream os(path, std::ios::binary);
    if (!os) {
        return false;
    }
    Header h;
    h.metric_id = metric_->id();
    h.n_points  = points_.size();
    h.n_nodes   = nodes_.size();
    os.write(reinterpret_cast<const char*>(&h), to_streamsize(sizeof(h)));
    os.write(reinterpret_cast<const char*>(points_.data()), to_streamsize(sizeof(Vec3) * points_.size()));
    os.write(reinterpret_cast<const char*>(index_.data()), to_streamsize(sizeof(uint32_t) * index_.size()));
    os.write(reinterpret_cast<const char*>(nodes_.data()), to_streamsize(sizeof(Node) * nodes_.size()));
    return os.good();
}


BallTree* BallTree::load(const std::string& path, Metric::MetricId metric_id) {
    std::ifstream is(path, std::ios::binary);
    if (!is) {
        return nullptr;
    }

    Header h{};
    is.read(reinterpret_cast<char*>(&h), sizeof(h));
    if (!is || h.magic != 0xBA115EED || h.version != 1 || h.metric_id != metric_id) {
        return nullptr;
    }

    std::vector<Vec3> pts(h.n_points);
    std::vector<uint32_t> idx(h.n_points);
    std::vector<Node> nodes(h.n_nodes);
    is.read(reinterpret_cast<char*>(pts.data()), to_streamsize(sizeof(Vec3) * pts.size()));
    is.read(reinterpret_cast<char*>(idx.data()), to_streamsize(sizeof(uint32_t) * idx.size()));
    is.read(reinterpret_cast<char*>(nodes.data()), to_streamsize(sizeof(Node) * nodes.size()));
    if (!is) {
        return nullptr;
    }

    return new BallTree(std::move(pts), std::move(idx), std::move(nodes), Metric::make(metric_id), /*leaf_size=*/16);
}


void BallTree::build() {
    index_.resize(points_.size());
    std::iota(index_.begin(), index_.end(), 0U);
    nodes_.clear();
    nodes_.reserve(points_.size() * 2);
    build_recurse(0, points_.size());
}


int BallTree::build_recurse(uint32_t start, uint32_t end) {
    Node node{};
    node.start = start;
    node.end   = end;
    // Center: normalized mean of unit vectors
    Vec3 mean{0, 0, 0};
    for (uint32_t i = start; i < end; ++i) {
        mean = mean + points_[index_[i]];
    }
    mean        = mean * (1.0 / std::max<uint32_t>(1, end - start));
    node.center = normalize(mean);

    // Radius under metric
    double r = 0.0;
    for (uint32_t i = start; i < end; ++i) {
        r = std::max(r, metric_->dist(node.center, points_[index_[i]]));
    }
    node.radius = r;

    auto my_id = static_cast<int>(nodes_.size());
    nodes_.push_back(node);

    auto count = end - start;
    if (count <= leaf_size_) {
        return my_id;
    }

    // Two-means style split: choose farthest pair (A,B), then assign
    uint32_t i0 = start + (count / 2);
    auto a      = farthest_index(index_[i0], start, end);
    auto b      = farthest_index(index_[a], start, end);
    if (a == b) {
        return my_id;  // degenerate
    }

    auto ca = points_[index_[a]];
    auto cb = points_[index_[b]];

    auto it_mid = std::partition(index_.begin() + start, index_.begin() + end, [&](uint32_t idx) {
        auto da = metric_->dist(points_[idx], ca);
        auto db = metric_->dist(points_[idx], cb);
        return da < db;
    });

    auto mid = static_cast<uint32_t>(std::distance(index_.begin(), it_mid));
    if (mid == start || mid == end) {
        mid = start + count / 2;
    }  // fallback

    int left            = build_recurse(start, mid);
    int right           = build_recurse(mid, end);
    nodes_[my_id].left  = left;
    nodes_[my_id].right = right;
    return my_id;
}


uint32_t BallTree::farthest_index(uint32_t pivot_idx, uint32_t start, uint32_t end) const {
    const auto& p = points_[pivot_idx];

    double best_d = -1.0;
    auto best_i   = start;
    for (auto i = start; i < end; ++i) {
        if (auto d = metric_->dist(points_[index_[i]], p); d > best_d) {
            best_d = d;
            best_i = i;
        }
    }
    return best_i;
}


void BallTree::radius_recurse(int node_id, const Vec3& q, double r, std::vector<uint32_t>& out) const {
    if (node_id < 0) {
        return;
    }

    const auto& nd = nodes_[node_id];
    if (double dc = metric_->dist(q, nd.center); dc - nd.radius > r) {
        return;
    }

    if (nd.is_leaf()) {
        for (uint32_t i = nd.start; i < nd.end; ++i) {
            uint32_t idx = index_[i];
            if (metric_->dist(q, points_[idx]) <= r) {
                out.push_back(idx);
            }
        }
        return;
    }

    radius_recurse(nd.left, q, r, out);
    radius_recurse(nd.right, q, r, out);
}


template <class Heap>
void BallTree::knn_recurse(int node_id, const Vec3& q, uint32_t k, Heap& best) const {
    if (node_id < 0) {
        return;
    }

    const auto& nd = nodes_[node_id];
    auto dc        = metric_->dist(q, nd.center);

    double worst = best.empty() ? std::numeric_limits<double>::infinity() : best.top().first;
    if (dc - nd.radius > worst) {
        return;
    }

    if (nd.is_leaf()) {
        for (uint32_t i = nd.start; i < nd.end; ++i) {
            auto idx = index_[i];
            auto d   = metric_->dist(q, points_[idx]);
            if (best.size() < k) {
                best.emplace(d, idx);
            }
            else if (d < best.top().first) {
                best.pop();
                best.emplace(d, idx);
            }
        }
        return;
    }
    const auto& L = nodes_[nd.left];
    const auto& R = nodes_[nd.right];

    auto dl = std::max(0., metric_->dist(q, L.center) - L.radius);
    auto dr = std::max(0., metric_->dist(q, R.center) - R.radius);

    auto first  = nd.left;
    auto second = nd.right;
    if (dr < dl) {
        std::swap(first, second);
    }

    knn_recurse(first, q, k, best);
    knn_recurse(second, q, k, best);
}


}  // namespace balls
