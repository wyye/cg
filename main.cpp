#include <iostream>
#include <cstdint>
#include <set>
#include <vector>
#include <memory>
#include <fstream>
#include <tuple>
#include <map>
#include <cassert>

#include "fixedreal.h"

struct Segment;

struct Point {
    explicit Point (FixedReal<6> x = makeFixedReal(0),
                    FixedReal<6> y = makeFixedReal(0)) :
        x(x), y(y)
    {}
    FixedReal<6> x, y;
};
struct SegPoint {

    std::set<Segment*> L;
    std::set<Segment*> U;
    std::set<Segment*> C;
};

bool
operator < (const Point & p1, const Point & p2) {
    return (p1.y > p2.y) || ((p1.y == p2.y) && (p1.x < p2.x));
}

bool
operator <= (const Point & p1, const Point & p2) {
    return !(p2 < p1);
}

bool
operator >= (const Point & p1, const Point & p2) {
    return !(p1 < p2);
}

bool
operator > (const Point & p1, const Point & p2) {
    return p2 < p1;
}

bool
operator == (const Point & p1, const Point & p2) {
    return !(p1 < p2) && !(p2 < p1);
}

bool
operator != (const Point & p1, const Point & p2) {
    return !(p1 == p2);
}

bool
collinear(const Point & a, const Point & b, const Point & c) {
    return (b.x - a.x) * (c.y - a.y) == (b.y - a.y) * (c.x - a.x);
}

// point a in between
bool
between(FixedReal<6> a, FixedReal<6> b, FixedReal<6> c) {
    return (b <= a) && (a <= c);
}

struct Segment {
    Segment(const Point* p1, const Point* p2) :
        u(p1), d(p2)
    {}
    const Point* u;
    const Point* d;

    bool contains(const Point & point) {
        return collinear(point, *u, *d) &&
                (u->x != d->x ?
                     between(point.x, u->x, d->x) :
                     between(point.y, u->y, d->y));
    }
};


bool
compareSegments(const Segment & a, const Segment & b, const FixedReal<6> y) {
    // find point of intersection with horisontal line y
    // (x - x1) * (y2 - y1) = (x2 - x1) * (y - y1)
    FixedReal<6> a_x = (a.u->y == a.d->y) ? a.u->x :
                 (y - a.u->y) * (a.d->x - a.u->x) / (a.d->y - a.u->y) + a.u->x;
    FixedReal<6> b_x = (b.u->y == b.d->y) ? b.u->x :
                 (y - b.u->y) * (b.d->x - b.u->x) / (b.d->y - b.u->y) + b.u->x;
    if (a_x == b_x) { // a and b from common point on y
        // compare angles
        if (a.d->y == y) { // if a horisontal, b can't be horisontal, so its angle less, b < a
            return false;
        } else if (b.d->y == y) { // if b horisontal, a can't be horisontal, so its angle less, a < b
            return true;
        }
        FixedReal<6> a_tg = (a.d->x - a_x) / (a.d->y - y);
        FixedReal<6> b_tg = (b.d->x - b_x) / (b.d->y - y);
        return a_tg < b_tg;
    } else { // a and b from different points on y
        return a_x < b_x;
    }
}

using Points = std::map<Point, std::unique_ptr<SegPoint>>;
using Segments = std::vector<std::unique_ptr<Segment>>;
struct SegmentComparator {
    bool operator() (const Segment* a, const Segment* b) {
        return compareSegments(*a, *b, y);
    }
    FixedReal<6> y = makeFixedReal(0);
};
using State = std::set<Segment*, SegmentComparator>;

bool
findIntersection(const Segment & sa, const Segment & sb, Point & p) {
    FixedReal<6> cross11 = (sa.u->x - sb.u->x)*(sb.d->y - sb.u->y) - (sa.u->y - sb.u->y)*(sb.d->x - sb.u->x);
    FixedReal<6> cross12 = (sa.d->x - sb.u->x)*(sb.d->y - sb.u->y) - (sa.d->y - sb.u->y)*(sb.d->x - sb.u->x);
    FixedReal<6> cross22 = (sb.u->x - sa.u->x)*(sa.u->y - sa.d->y) - (sb.u->y - sa.u->y)*(sa.u->x - sa.d->x);
    FixedReal<6> cross21 = (sb.d->x - sa.u->x)*(sa.u->y - sa.d->y) - (sb.d->y - sa.u->y)*(sa.u->x - sa.d->x);
    if (cross11*cross12 <= makeFixedReal(0) && cross21*cross22 <= makeFixedReal(0)) {
//         solution of
//
//         x == xua + (xda - xua) ta && y == yua + (yda - yua) ta &&
//         x == xub + (xdb - xub) tb && y == yub + (ydb - yub) tb
//
//         is:
//
//           ta -> -((-xub yua + xdb yua + xua yub - xdb yub - xua ydb + xub ydb)/
//                  (xub yua - xdb yua - xua yub + xda yub - xub yda + xdb yda + xua ydb - xda ydb))
//           tb -> -((xub yua - xda yua - xua yub + xda yub + xua yda - xub yda)/
//                  (-xub yua + xdb yua + xua yub - xda yub + xub yda - xdb yda - xua ydb + xda ydb))
//           x -> -((xub xda yua - xda xdb yua - xua xdb yub + xda xdb yub -
//                   xua xub yda + xua xdb yda + xua xub ydb - xub xda ydb) /
//                  (-xub yua + xdb yua + xua yub - xda yub + xub yda - xdb yda - xua ydb + xda ydb))
//           y -> -((-xda yua yub + xdb yua yub + xua yub yda - xdb yub yda -
//                   xub yua ydb + xda yua ydb - xua yda ydb + xub yda ydb) /
//                  (xub yua - xdb yua - xua yub + xda yub - xub yda + xdb yda + xua ydb - xda ydb))
        p.x = - ((sb.u->x * sa.d->x * sa.u->y -
                  sa.d->x * sb.d->x * sa.u->y -
                  sa.u->x * sb.d->x * sb.u->y +
                  sa.d->x * sb.d->x * sb.u->y -
                  sa.u->x * sb.u->x * sa.d->y +
                  sa.u->x * sb.d->x * sa.d->y +
                  sa.u->x * sb.u->x * sb.d->y -
                  sb.u->x * sa.d->x * sb.d->y) /
                 (-sb.u->x * sa.u->y +
                  sb.d->x * sa.u->y +
                  sa.u->x * sb.u->y -
                  sa.d->x * sb.u->y +
                  sb.u->x * sa.d->y -
                  sb.d->x * sa.d->y -
                  sa.u->x * sb.d->y +
                  sa.d->x * sb.d->y));
        p.y = - ((-sa.d->x * sa.u->y * sb.u->y +
                  sb.d->x * sa.u->y * sb.u->y +
                  sa.u->x * sb.u->y * sa.d->y -
                  sb.d->x * sb.u->y * sa.d->y -
                  sb.u->x * sa.u->y * sb.d->y +
                  sa.d->x * sa.u->y * sb.d->y -
                  sa.u->x * sa.d->y * sb.d->y +
                  sb.u->x * sa.d->y * sb.d->y) /
                 (sb.u->x * sa.u->y -
                  sb.d->x * sa.u->y -
                  sa.u->x * sb.u->y +
                  sa.d->x * sb.u->y -
                  sb.u->x * sa.d->y +
                  sb.d->x * sa.d->y +
                  sa.u->x * sb.d->y -
                  sa.d->x * sb.d->y));
        return true;
    }
    return false;
}

void
determineLineIntersection(Points& points,
                          State& state,
                          SegmentComparator& scmp) {
    for(Points::const_iterator it = points.begin(); it != points.end(); ++it) {

        scmp.y = it->first.y;
        std::set<Segment*> & L = it->second->L;
        std::set<Segment*> & U = it->second->U;
        std::set<Segment*> & C = it->second->C;
        State::iterator up_left = state.end();
        State::iterator up_right = state.end();

        for (auto it_l = L.begin(); it_l != L.end(); ++it_l) {
            up_right = state.erase(state.find(*it_l));
        }
        for (auto it_c = C.begin(); it_c != C.end(); ++it_c) {
            up_right = state.erase(state.find(*it_c));
        }
        if (!state.empty() && up_right != state.begin()) {
            up_left = std::prev(up_right);
        }

        State::iterator it_max, it_min;
        if (U.empty() && C.empty()) {
            if (up_left != state.end() && up_right != state.end()) {
                Point p;
                if (findIntersection(**up_left, **up_right, p)) {
                    auto inserted = points.insert(std::make_pair(p, std::make_unique<SegPoint>()));
                    if (p > it->first) {
                        inserted.first->second->C.insert(*up_left);
                        inserted.first->second->C.insert(*up_right);
                    }
                }
            }
        } else {
            Segment* s_min;
            Segment* s_max;
            if (!U.empty()) {
                s_max = s_min = *U.begin();
                it_max = it_min = state.insert(*U.begin()).first;
            } else {
                s_max = s_min = *C.begin();
                it_max = it_min = state.insert(*C.begin()).first;
            }

            for (auto it_u = U.begin(); it_u != U.end(); ++it_u) {
                auto cur_it = state.insert(*it_u);
                if (scmp(*it_u, s_min)) {
                    it_min = cur_it.first;
                    s_min = *it_u;
                }
                if (scmp(s_max, *it_u)) {
                    it_max = cur_it.first;
                    s_max = *it_u;
                }

            }
            for (auto it_c = C.begin(); it_c != C.end(); ++it_c) {
                auto cur_it = state.insert(*it_c);
                if (scmp(*it_c, s_min)) {
                    it_min = cur_it.first;
                    s_min = *it_c;
                }
                if (scmp(s_max, *it_c)) {
                    it_max = cur_it.first;
                    s_max = *it_c;
                }
            }
            State::iterator it_left = state.end();
            if (!state.empty() && it_min != state.begin()) {
                it_left = std::prev(it_min);
            }

            State::iterator it_right = state.end();
            if (it_max != state.end()) {
                it_right = std::next(it_max);
            }

            if (it_left != state.end() && it_min != state.end()) {
                Point p;
                if (findIntersection(**it_left, **it_min, p)) {
                    auto inserted = points.insert(std::make_pair(p, std::make_unique<SegPoint>()));
                    if (p > it->first) {
                        inserted.first->second->C.insert(*it_left);
                        inserted.first->second->C.insert(*it_min);
                    }
                }
            }

            if (it_right != state.end() && it_max != state.end()) {
                Point p;
                if (findIntersection(**it_right, **it_max, p)) {
                    auto inserted = points.insert(std::make_pair(p, std::make_unique<SegPoint>()));
                    if (p > it->first) {
                        inserted.first->second->C.insert(*it_right);
                        inserted.first->second->C.insert(*it_max);
                    }
                }
            }
        }
    }
}

void
readData(std::ifstream & ifs,
         std::map<Point, std::unique_ptr<SegPoint>> & points,
         std::vector<std::unique_ptr<Segment>> & segments) {

    while (1) {
        FixedReal<6> x1, y1, x2, y2;
        ifs >> x1 >> y1 >> x2 >> y2;
        if (ifs.eof()) break;

        Point p1(x1, y1);
        Point p2(x2, y2);
        assert(p1 != p2);

        if (p1 > p2) {
            std::swap(p1, p2);
        }
        auto it1 = points.insert(std::make_pair(p1, std::make_unique<SegPoint>()));
        auto it2 = points.insert(std::make_pair(p2, std::make_unique<SegPoint>()));



        segments.emplace_back(std::make_unique<Segment>(&it1.first->first, &it2.first->first));
        it1.first->second->U.insert(segments.back().get());
        it2.first->second->L.insert(segments.back().get());
    }
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Launch : " << argv[0] << " file.txt" << std::endl;
        exit(1);
    }
    std::ifstream ifs(argv[1]);

    if (!ifs.is_open()) {
        std::cout << "Can't open file" << std::endl;
        exit(1);
    }

    Points points;
    Segments segments;
    readData(ifs, points, segments);

    SegmentComparator scmp;
    State state(scmp);

    determineLineIntersection(points, state, scmp);

    for(Points::const_iterator it = points.begin(); it != points.end(); ++it) {
        if (!it->second->C.empty()) {
            std::cout << it->first.x << " " << it->first.y << std::endl;
        }
    }

    test();
}
