#include <iostream>
#include <cstdint>
#include <set>
#include <vector>
#include <memory>
#include <fstream>
#include <tuple>
#include <map>
#include <algorithm>
#include <cassert>
#include <cmath>

#include "fixedreal.h"
#include "logger.h"

struct Segment;

struct Point {
    explicit Point (double x = 0,
                    double y = 0) :
        x(x), y(y)
    {}
    double x, y;
};

std::ostream& operator<<(std::ostream & stream, const Point & point) {
    stream << "(" << point.x << ", " << point.y << ")";
    return stream;
}


struct SegPoint {

    std::set<Segment*> L;
    std::set<Segment*> U;
    std::set<Segment*> C;
};


bool
operator == (const Point & p1, const Point & p2) {
    return std::abs(p1.x - p2.x) <= 0.001 && std::abs(p1.y - p2.y) <= 0.001;
}


bool
operator != (const Point & p1, const Point & p2) {
    return !(p1 == p2);
}

bool
operator < (const Point & p1, const Point & p2) {
    if (p1 == p2) return false;
    if (std::abs(p1.y - p2.y) <= 0.001) {
        return p1.x < p2.x;
    } else {
        return p1.y < p2.y;
    }
}

bool
operator > (const Point & p1, const Point & p2) {
    if (p1 == p2) return false;
    if (std::abs(p1.y - p2.y) <= 0.001) {
        return p1.x > p2.x;
    } else {
        return p1.y > p2.y;
    }
}

bool
operator <= (const Point & p1, const Point & p2) {
    return p2 < p1 || p1 == p2;
}

bool
operator >= (const Point & p1, const Point & p2) {
    return p1 > p2 || p1 == p2;
}

bool
collinear(const Point & a, const Point & b, const Point & c) {
    return (b.x - a.x) * (c.y - a.y) == (b.y - a.y) * (c.x - a.x);
}

// point a in between
bool
between(double a, double b, double c) {
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

std::ostream& operator<<(std::ostream & stream, const Segment & segment) {
    stream << "[" << *segment.u << " " << *segment.d << "]";
    return stream;
}

/*
 *    ------------------------------------
 *                   |\
 *                   | \
 *                   |  \
 *                   |  *\
 *                   |**  \
 *                   |     \
 *                   |      \
 *                   |       \
 *                   |        \
 *                   | angle   \
 *                   |          \
 */

int
compareSegments(const Segment & a, const Segment & b, const double x, const double y) {
    // find point of intersection with horisontal line y
    // (x - x1) * (y2 - y1) = (x2 - x1) * (y - y1)
    double a_x = std::abs(a.u->y - a.d->y) < 0.001 ? a.u->x :
                 (y - a.u->y) * (a.d->x - a.u->x) / (a.d->y - a.u->y) + a.u->x;
    double b_x = std::abs(b.u->y - b.d->y) < 0.001 ? b.u->x :
                 (y - b.u->y) * (b.d->x - b.u->x) / (b.d->y - b.u->y) + b.u->x;
    if (std::abs(a_x - b_x) < 0.001) { // a and b from common point on y
        // compare angles
        if (std::abs(y - a.d->y) < 0.001 && std::abs(y - b.d->y) < 0.001) {
            return 0;
        } else if (std::abs(y - a.d->y) < 0.001) { // if a horisontal, b can't be horisontal, so its angle less, b < a
            return 1;
        } else if (std::abs(y - b.d->y) < 0.001) { // if b horisontal, a can't be horisontal, so its angle less, a < b
            return -1;
        }

        int res;
        if (a_x <= x) {
            double a_tg = (a.d->x - a_x) / (y - a.d->y);
            double b_tg = (b.d->x - b_x) / (y - b.d->y); // delta y will be always positive
            res = std::abs(a_tg - b_tg) < 0.001 ? 0 : a_tg < b_tg ? -1 : 1;
        } else {
            double a_tg = (a.u->x - a_x) / (a.u->y - y);
            double b_tg = (b.u->x - b_x) / (b.u->y - y); // delta y will be always positive
            res = std::abs(a_tg - b_tg) < 0.001 ? 0 : a_tg < b_tg ? -1 : 1;
        }

        return res;
    } else { // a and b from different points on y
        return std::abs(a_x - b_x) < 0.001 ? 0 : (a_x < b_x) ? -1 : 1;
    }
}

using Points = std::map<Point, std::unique_ptr<SegPoint>>;
using Segments = std::vector<std::unique_ptr<Segment>>;
struct SegmentComparator {
    bool operator() (const Segment* a, const Segment* b) {
        return compareSegments(*a, *b, *x, *y) < 0;
    }
    std::shared_ptr<double> x = std::make_shared<double>(0);
    std::shared_ptr<double> y = std::make_shared<double>(0);
};
using State = std::set<Segment*, SegmentComparator>;

bool
findIntersection(const Segment & sa, const Segment & sb, Point & p) {
    double cross11 = (sa.u->x - sb.u->x)*(sb.d->y - sb.u->y) - (sa.u->y - sb.u->y)*(sb.d->x - sb.u->x);
    double cross12 = (sa.d->x - sb.u->x)*(sb.d->y - sb.u->y) - (sa.d->y - sb.u->y)*(sb.d->x - sb.u->x);
    double cross22 = (sb.u->x - sa.u->x)*(sa.u->y - sa.d->y) - (sb.u->y - sa.u->y)*(sa.u->x - sa.d->x);
    double cross21 = (sb.d->x - sa.u->x)*(sa.u->y - sa.d->y) - (sb.d->y - sa.u->y)*(sa.u->x - sa.d->x);
    if (cross11*cross12 <= 0 && cross21*cross22 <= 0) {
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
    *scmp.x = points.begin()->first.x;
    *scmp.y = points.begin()->first.y;
    for(Points::const_iterator it = points.begin(); it != points.end(); ++it) {

        LOG::INFO("Process point %", it->first);
        //*scmp.y = it->first.y;
        std::set<Segment*> & L = it->second->L;
        std::set<Segment*> & U = it->second->U;
        std::set<Segment*> & C = it->second->C;
        State::iterator up_left = state.end();
        State::iterator up_right = state.end();


        if (L.empty() && C.empty()) {
            LOG::INFO("L && C are empty, don't remove");
        } else {
            LOG::SECTION sect("Remove L and C");
            for (Segment* segment : state) {
                LOG::DEBUG("state before deletion : %", *segment);
            }
            for (auto it_l = L.begin(); it_l != L.end(); ++it_l) {
                LOG::DEBUG("remove % from state (L)", **it_l);
                up_right = state.erase(state.find(*it_l));
                for (Segment* segment : state) {
                    LOG::DEBUG("state after deletion(1) : %, %", *segment, segment);
                }
            }

            for (auto it_c = C.begin(); it_c != C.end(); ++it_c) {
                LOG::DEBUG("remove % , % from state (C)", **it_c, *it_c);
                State::iterator found = state.find(*it_c);
                if (found == state.end()) {
                    LOG::INFO("% Already deleted", **it_c, *it_c);
                } else {
                    up_right = state.erase(found);
                }
                for (Segment* segment : state) {
                    LOG::DEBUG("state after deletion(2) : %, %", *segment, segment);
                }
            }
            if (!state.empty() && up_right != state.begin()) {
                up_left = std::prev(up_right);
            }
            for (Segment* segment : state) {
                LOG::DEBUG("state after deletion : %", *segment);
            }
        }


        if (up_right != state.end()) {
            LOG::INFO("Segment right to the last removed from L & C : %", **up_right);
        } else {
            LOG::INFO("Segment right to the last removed from L & C not found");
        }
        if (up_left != state.end()) {
            LOG::INFO("Segment left to the last removed from L & C : %", **up_left);
        } else {
            LOG::INFO("Segment left to the last removed from L & C not found");
        }

        *scmp.x = it->first.x;
        *scmp.y = it->first.y;

        State::iterator it_max, it_min;
        if (U.empty() && C.empty()) {
            LOG::INFO("U & C are empty");
            if (up_left != state.end() && up_right != state.end()) {
                LOG::INFO("Finding intersection between up_left & up_right ...");
                Point p;
                if (findIntersection(**up_left, **up_right, p)) {
                    LOG::INFO("Found intersection between up_left & up_right %", p);
                    auto inserted = points.insert(std::make_pair(p, std::make_unique<SegPoint>()));
                    if (inserted.second) {
                        LOG::INFO("Found new point");
                    }
                    if (p >= it->first) {
                        LOG::INFO("New intersection was not processed yet");
                        inserted.first->second->C.insert(*up_left);
                        inserted.first->second->C.insert(*up_right);
                    }
                } else {
                    LOG::INFO("No intersection between up_left & up_right found");
                }
            }
        } else {
            {
                LOG::SECTION sect("Inserting U and C");
                if (state.empty()) {
                    LOG::DEBUG("state before insertion : empty");
                }
                for (Segment* segment : state) {
                    LOG::DEBUG("state before insertion : %", *segment);
                }
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
                    LOG::DEBUG("Inserting % in state (U)", **cur_it.first);
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
                    LOG::DEBUG("Inserting % in state (C)", **cur_it.first);
                    if (scmp(*it_c, s_min)) {
                        it_min = cur_it.first;
                        s_min = *it_c;
                    }
                    if (scmp(s_max, *it_c)) {
                        it_max = cur_it.first;
                        s_max = *it_c;
                    }
                }

                for (Segment* segment : state) {
                    LOG::DEBUG("state after insertion : %", *segment);
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

            if (it_min == state.end()) {
                LOG::INFO("Leftmost element from U & C : not found");
            } else {
                LOG::INFO("Leftmost element from U & C : %", **it_min);
            }
            if (it_max == state.end()) {
                LOG::INFO("Rightmost element from U & C : not found");
            } else {
                LOG::INFO("Rightmost element from U & C : %", **it_max);
            }
            if (it_left == state.end()) {
                LOG::INFO("Left from leftmost element from U & C : not found");
            } else {
                LOG::INFO("Left from leftmost element from U & C : %", **it_left);
            }
            if (it_right == state.end()) {
                LOG::INFO("Right from rightmost element from U & C : not found");
            } else {
                LOG::INFO("Right from rightmost element from U & C : %", **it_right);
            }

            if (it_left != state.end() && it_min != state.end()) {
                Point p;
                if (findIntersection(**it_left, **it_min, p)) {
                    LOG::INFO("Found intersection between up_left & up_right %", p);
                    auto inserted = points.insert(std::make_pair(p, std::make_unique<SegPoint>()));
                    if (inserted.second) {
                        LOG::INFO("Found new point ");
                    }
                    if (p >= it->first) {
                        LOG::INFO("New intersection was not processed yet");
                        inserted.first->second->C.insert(*it_left);
                        inserted.first->second->C.insert(*it_min);
                    }
                } else {
                    LOG::INFO("No intersection between up_left & up_right found");
                }
            }

            if (it_right != state.end() && it_max != state.end()) {
                Point p;
                if (findIntersection(**it_right, **it_max, p)) {
                    LOG::INFO("Found intersection between up_left & up_right %", p);
                    auto inserted = points.insert(std::make_pair(p, std::make_unique<SegPoint>()));
                    if (inserted.second) {
                        LOG::INFO("Found new point ");
                    }
                    if (p >= it->first) {
                        LOG::INFO("New intersection was not processed yet");
                        inserted.first->second->C.insert(*it_right);
                        inserted.first->second->C.insert(*it_max);
                    }
                } else {
                    LOG::INFO("No intersection between up_left & up_right found");
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
        double x1, y1, x2, y2;
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
    if (argc != 3) {
        std::cout << "Launch : " << argv[0] << " in.txt out.txt" << std::endl;
        exit(1);
    }
    std::ifstream ifs(argv[1]);
    std::ofstream ofs(argv[2]);

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
            ofs << it->first.x << " " << it->first.y << std::endl;
        }
    }
    ifs.close();
    ofs.close();
}
