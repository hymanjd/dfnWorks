#ifndef _QUASI2DDOMIAN_H
#define _QUASI2DDOMIAN_H

#include "structures.h"
#include <vector>

void quasi2dDomain(std::vector<Poly> &acceptedPolys, std::vector<IntPoints> &intPts, std::vector<Point> triplePoints, Stats  &pstats);
bool inQuasi2DDomain(double x, double y);

#endif
