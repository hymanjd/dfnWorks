#ifndef _wiggle_h_
#define _wiggle_h_

#include <vector>
#include "structures.h"

/******************************************************************************
 * \file wiggle.h
 * \brief Small rigid shifts ("wiggles") that try to fit a fracture FRAM rejected.
 *
 * FRAM rejects a fracture when it would create a feature smaller than h. Those
 * bad configurations are thin (about h wide), so moving the fracture by a few h
 * usually clears them. wigglePoly() tries a fixed sequence of translations and
 * keeps the smallest one that FRAM accepts. Shape, size and orientation are
 * never changed, and no random numbers are used.
 ******************************************************************************/

/*!
 * \brief Try to insert a FRAM-rejected fracture by translating it a few h.
 *
 * Shift amplitudes run h, 2h, ..., wiggleMaxSteps*h (outer loop) over the six
 * directions +/-normal and +/- two in-plane axes (inner loop), so the smallest
 * working shift wins. In-plane directions are tried first for short
 * intersection and close-to-edge rejections (-2, -6); the normal directions
 * are tried first for every other reject code.
 *
 * A shift is only accepted if the fracture still intersects every fracture it
 * intersected where the user placed it, so a wiggle never drops a connection.
 *
 * Each attempt rebuilds the polygon from \p pristineVertices, then runs
 * domainTruncation(), createBoundingBox() and intersectionChecking(). On
 * success the intersections are already committed by intersectionChecking();
 * the caller still pushes \p newPoly to the accepted list, as it does for a
 * fracture accepted on the first try.
 *
 * FRAM rejection counters are restored to their value on entry, so failed
 * wiggle attempts do not inflate the rejection statistics.
 *
 * \param newPoly           Rejected polygon. Updated in place on success;
 *                          restored to its unshifted, untruncated state on failure.
 * \param pristineVertices  Vertices of \p newPoly before domain truncation.
 * \param pristineNumNodes  Number of vertices in \p pristineVertices.
 * \param firstRejectCode   Reject code from the original intersectionChecking().
 * \param id                1-based user fracture id, for reporting.
 * \param userFractureType  -1 ellipse, -2 rectangle, -3 polygon, for reporting.
 * \param acceptedPoly      Accepted polygons.
 * \param intPts            Accepted intersections.
 * \param pstats            Program statistics.
 * \param triplePoints      Accepted triple intersection points.
 * \return 0 if a shifted copy was accepted, otherwise \p firstRejectCode.
 */
int wigglePoly(Poly &newPoly, const std::vector<double> &pristineVertices, int pristineNumNodes,
               int firstRejectCode, int id, int userFractureType,
               std::vector<Poly> &acceptedPoly, std::vector<IntPoints> &intPts,
               Stats &pstats, std::vector<Point> &triplePoints);

#endif // _wiggle_h_
