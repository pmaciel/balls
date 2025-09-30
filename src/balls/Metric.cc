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


#include "balls/Metric.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "balls/Vec3.h"


namespace balls {


// (Haversine-equivalent)
struct GreatCircleMetric final : Metric {
    GreatCircleMetric() : Metric(MetricId::GreatCircle) {}

    double dist(const Vec3& a, const Vec3& b) const override {
        auto c = std::max(-1., std::min(1., dot(a, b)));
        return /* R * */ std::acos(c);
    }
};


struct EuclideanMetric final : Metric {
    EuclideanMetric() : Metric(MetricId::Euclidean) {}

    double dist(const Vec3& a, const Vec3& b) const override { return (a - b).norm(); }
};


Metric* Metric::make(MetricId id) {
    return id == MetricId::GreatCircle ? static_cast<Metric*>(new GreatCircleMetric)
           : id == MetricId::Euclidean ? static_cast<Metric*>(new EuclideanMetric)
                                       : throw std::invalid_argument("Metric: unknown");
}


}  // namespace balls
