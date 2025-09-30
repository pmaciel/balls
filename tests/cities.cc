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


#include <cassert>
#include <iostream>
#include <memory>

#include "balls/BallTree.h"
#include "balls/Metric.h"


int main() {
    using namespace balls;

    auto print_knn = [](const std::vector<std::pair<double, uint32_t>>& pairs, const std::vector<Vec3>& pts) {
        std::cout << "kNN results (dist, idx):\n";
        for (const auto& p : pairs) {
            std::cout << "  d=" << p.first << ", idx=" << p.second << "  xyz=(" << pts[p.second].x << ","
                      << pts[p.second].y << "," << pts[p.second].z << ")\n";
        }
    };

    auto latlon_to_unit = [](double lat, double lon) {
        double clat = std::cos(lat);
        return Vec3{clat * std::cos(lon), clat * std::sin(lon), std::sin(lat)};
    };

    auto deg2rad = [](double d) { return d * M_PI / 180.; };

    struct City {
        const char* name;
        double lat_deg, lon_deg;
    };

    std::vector<City> cities = {
        {"London", 51.5074, -0.1278},           //
        {"Paris", 48.8566, 2.3522},             //
        {"New_York", 40.7128, -74.0060},        //
        {"San_Francisco", 37.7749, -122.4194},  //
        {"Tokyo", 35.6762, 139.6503},           //
        {"Sydney", -33.8688, 151.2093},         //
        {"Cape_Town", -33.9249, 18.4241},       //
        {"Rio", -22.9068, -43.1729},            //
        {"Singapore", 1.3521, 103.8198},        //
        {"Toronto", 43.6532, -79.3832},         //
    };

    Vec3 berlin = latlon_to_unit(deg2rad(52.5200), deg2rad(13.4050));


    std::vector<Vec3> P;
    P.reserve(cities.size());
    for (auto& c : cities) {
        P.push_back(latlon_to_unit(deg2rad(c.lat_deg), deg2rad(c.lon_deg)));
    }

    // ---------- Great-circle

    {
        std::unique_ptr<Metric> metric(Metric::make(Metric::GreatCircle));
        BallTree gc_tree(P, metric.get(), /*leaf_size=*/4);

        auto gc_knn = gc_tree.knn(berlin, 3);
        std::cout << "\nGreat-circle metric (km): kNN from Berlin\n";
        print_knn(gc_knn, gc_tree.points());

        std::cout << "Within (great-circle):";
        for (auto i : gc_tree.radius_search(berlin, 1000.0)) {
            std::cout << " " << cities[i].name;
        }
        std::cout << "\n";

        const std::string gc_file = "cities_gc_poly.btree";
        if (!gc_tree.save(gc_file)) {
            std::cerr << "Failed to save " << gc_file << "\n";
            return 1;
        }

        // Caller supplies matching metric when loading
        std::unique_ptr<Metric> metric_for_load(Metric::make(Metric::GreatCircle));
        std::unique_ptr<BallTree> btree(BallTree::load(gc_file, metric_for_load->id()));
        assert(btree);

        std::cout << "Reloaded GC tree; points=" << btree->points().size() << ", nodes=" << btree->nodes().size()
                  << "\n";
    }

    // ---------- Euclidean (chord) ----------
    std::unique_ptr<Metric> eu_metric(Metric::make(Metric::Euclidean));

    BallTree eu_tree(P, eu_metric.get(), /*leaf_size=*/4);

    auto eu_knn = eu_tree.knn(berlin, 3);
    std::cout << "\nEuclidean (chord) metric: kNN from Berlin\n";
    print_knn(eu_knn, eu_tree.points());

    double arc           = 0.1;                      // radians on unit sphere
    double chord_on_unit = 2.0 * std::sin(arc / 2);  // unit sphere chord

    std::cout << "Within (chord):";
    for (auto i : eu_tree.radius_search(berlin, chord_on_unit)) {
        std::cout << " " << cities[i].name;
    }
    std::cout << "\n";

    const std::string eu_file = "cities_eu_poly.btree";
    if (!eu_tree.save(eu_file)) {
        std::cerr << "Failed to save " << eu_file << "\n";
        return 1;
    }

    std::unique_ptr<Metric> eu_metric_for_load(Metric::make(Metric::Euclidean));
    std::unique_ptr<BallTree> eu_loaded(BallTree::load(eu_file, eu_metric_for_load->id()));
    assert(eu_loaded);
    std::cout << "Reloaded EU tree; points=" << eu_loaded->points().size() << ", nodes=" << eu_loaded->nodes().size()
              << "\n";

    return 0;
}
