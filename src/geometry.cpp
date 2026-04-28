//============================================================
//
//      Type:        FWH source file
//
//      Author:      Tommaso Bellosta on 14/04/26.
//                   Dipartimento di Scienze e Tecnologie Aerospaziali
//                   Politecnico di Milano
//                   Via La Masa 34, 20156 Milano, ITALY
//                   e-mail: tommaso.bellosta@polimi.it
//
//      Copyright:   2026, authors above and the FWH contributors.
//                   This software is distributed under the MIT license, see LICENSE.txt
//
//============================================================
#include <cmath>
#include <map>

#include "geometry.h"

double dot(const Vect3& v1, const Vect3& v2) {
    double out = 0;
    for (int i = 0; i < 3; i++) out += v1[i] * v2[i];
    return out;
}


Vect3 cross(const Vect3& v1, const Vect3& v2) {
    Vect3 out;
    out[0] = v1[1]*v2[2] - v1[2]*v2[1];
    out[1] = v1[2]*v2[0] - v1[0]*v2[2];
    out[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return out;
}

double norm(const Vect3& v) {
    double out = 0;
    for (int i = 0; i < 3; i++) out += v[i] * v[i];
    return std::sqrt(out);
}


void generate_icosphere(std::vector<Vect3>& nodes,
                        std::vector<std::vector<size_t>>& elements,
                        int n_refine,
                        double radius) {
    double t = (1.0 + std::sqrt(5.0)) / 2.0;

    nodes = {
        Vect3(-1,  t, 0), Vect3( 1,  t, 0), Vect3(-1, -t, 0), Vect3( 1, -t, 0),
        Vect3( 0, -1,  t), Vect3( 0,  1,  t), Vect3( 0, -1, -t), Vect3( 0,  1, -t),
        Vect3( t,  0, -1), Vect3( t,  0,  1), Vect3(-t,  0, -1), Vect3(-t,  0,  1)
    };

    for (auto& v : nodes) {
        double r = norm(v);
        for (int d = 0; d < 3; ++d)
            v[d] /= r;
    }

    elements = {
        {0,11,5},  {0,5,1},   {0,1,7},   {0,7,10},  {0,10,11},
        {1,5,9},   {5,11,4},  {11,10,2},  {10,7,6},  {7,1,8},
        {3,9,4},   {3,4,2},   {3,2,6},    {3,6,8},   {3,8,9},
        {4,9,5},   {2,4,11},  {6,2,10},   {8,6,7},   {9,8,1}
    };

    for (int ref = 0; ref < n_refine; ++ref) {

        std::map<std::pair<size_t,size_t>, size_t> mid_cache;
        std::vector<std::vector<size_t>> new_elems;

        auto midpoint = [&](size_t a, size_t b) -> size_t {

            auto key = std::make_pair(std::min(a,b), std::max(a,b));
            auto it = mid_cache.find(key);
            if (it != mid_cache.end())
                return it->second;

            Vect3 m;
            for (int d = 0; d < 3; ++d)
                m[d] = 0.5 * (nodes[a][d] + nodes[b][d]);

            double r = norm(m);
            for (int d = 0; d < 3; ++d)
                m[d] /= r;

            size_t idx = nodes.size();
            nodes.push_back(m);
            mid_cache[key] = idx;
            return idx;
        };

        for (auto& e : elements) {
            size_t a = e[0], b = e[1], c = e[2];
            size_t ab = midpoint(a, b);
            size_t bc = midpoint(b, c);
            size_t ca = midpoint(c, a);

            new_elems.push_back({a,  ab, ca});
            new_elems.push_back({b,  bc, ab});
            new_elems.push_back({c,  ca, bc});
            new_elems.push_back({ab, bc, ca});
        }

        elements = new_elems;
    }

    for (auto& e : elements) {
        const Vect3& xa = nodes[e[0]];
        const Vect3& xb = nodes[e[1]];
        const Vect3& xc = nodes[e[2]];

        Vect3 centroid;
        for (int d = 0; d < 3; ++d)
            centroid[d] = (xa[d] + xb[d] + xc[d]) / 3.0;

        Vect3 v1, v2;
        for (int d = 0; d < 3; ++d) {
            v1[d] = xb[d] - xa[d];
            v2[d] = xc[d] - xa[d];
        }

        Vect3 n = cross(v1, v2);

        if (dot(n, centroid) < 0.0)
            std::swap(e[1], e[2]);
    }

    for (auto& v : nodes) {
        for (int d = 0; d < 3; ++d)
            v[d] *= radius;
    }
}
