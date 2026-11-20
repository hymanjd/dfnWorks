#include "wiggle.h"
#include <cmath>
#include <string>
#include "computationalGeometry.h"
#include "domain.h"
#include "input.h"
#include "layers.h"
#include "logFile.h"

/******************************************************************************
 * Rebuild newPoly from its pristine (untruncated) vertices plus a shift, and
 * clear everything a previous insertion attempt may have set.
 ******************************************************************************/
static void resetToShifted(Poly &newPoly, const std::vector<double> &pristineVertices,
                           int pristineNumNodes, const double *pristineTranslation,
                           bool pristineTruncated, const double *shift) {
    delete[] newPoly.vertices;
    newPoly.vertices = new double[3 * pristineNumNodes];
    newPoly.numberOfNodes = pristineNumNodes;

    for (int k = 0; k < pristineNumNodes; k++) {
        int idx = 3 * k;
        newPoly.vertices[idx]     = pristineVertices[idx]     + shift[0];
        newPoly.vertices[idx + 1] = pristineVertices[idx + 1] + shift[1];
        newPoly.vertices[idx + 2] = pristineVertices[idx + 2] + shift[2];
    }

    for (int k = 0; k < 3; k++) {
        newPoly.translation[k] = pristineTranslation[k] + shift[k];
    }

    for (int k = 0; k < 6; k++) {
        newPoly.faces[k] = 0;
    }

    newPoly.truncated = pristineTruncated;
    newPoly.groupNum = 0;
    newPoly.intersectionIndex.clear();
}

/******************************************************************************
 * True if the two polygons intersect, whether or not FRAM would accept the
 * intersection. Both polygons must have their bounding boxes set.
 ******************************************************************************/
static bool touches(Poly &a, Poly &b) {
    if (!checkBoundingBox(a, b)) {
        return false;
    }

    short flag;
    findIntersections(flag, a, b);
    return flag != 0;
}

/******************************************************************************
 * Build two unit vectors u, v that span the fracture plane.
 * u follows the first polygon edge that is long enough to give a direction;
 * v = normal x u.
 * Returns false if no usable edge was found (degenerate polygon).
 ******************************************************************************/
static bool inPlaneBasis(const Poly &newPoly, const std::vector<double> &pristineVertices,
                         int pristineNumNodes, double *u, double *v) {
    const double *n = newPoly.normal;

    for (int k = 0; k < pristineNumNodes; k++) {
        int a = 3 * k;
        int b = 3 * ((k + 1) % pristineNumNodes);
        double e[3] = {pristineVertices[b]     - pristineVertices[a],
                       pristineVertices[b + 1] - pristineVertices[a + 1],
                       pristineVertices[b + 2] - pristineVertices[a + 2]
                      };
        // Remove any component along the normal so u lies exactly in the plane
        double along = e[0] * n[0] + e[1] * n[1] + e[2] * n[2];
        e[0] -= along * n[0];
        e[1] -= along * n[1];
        e[2] -= along * n[2];
        double len = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);

        if (len > eps) {
            u[0] = e[0] / len;
            u[1] = e[1] / len;
            u[2] = e[2] / len;
            v[0] = n[1] * u[2] - n[2] * u[1];
            v[1] = n[2] * u[0] - n[0] * u[2];
            v[2] = n[0] * u[1] - n[1] * u[0];
            return true;
        }
    }

    return false;
}

/******************************************************************************
 * int wigglePoly()
 * See wiggle.h for full documentation.
 ******************************************************************************/
int wigglePoly(Poly &newPoly, const std::vector<double> &pristineVertices, int pristineNumNodes,
               int firstRejectCode, int id, int userFractureType,
               std::vector<Poly> &acceptedPoly, std::vector<IntPoints> &intPts,
               Stats &pstats, std::vector<Point> &triplePoints,
               int layerIndex, bool pristineTruncated) {
    // Only user-defined fractures (negative type code) are logged and recorded one by one
    const bool isUserFracture = (userFractureType < 0);
    const double zeroShift[3] = {0, 0, 0};
    const double pristineTranslation[3] = {newPoly.translation[0], newPoly.translation[1], newPoly.translation[2]};
    double u[3], v[3];

    if (pristineNumNodes < 3 || !inPlaneBasis(newPoly, pristineVertices, pristineNumNodes, u, v)) {
        return firstRejectCode;
    }

    // Fractures the polygon intersects where the user placed it. A wiggle is
    // only accepted if it still intersects every one of them, so fitting the
    // fracture in never costs a connection. On entry newPoly is the truncated
    // original with its bounding box set.
    std::vector<unsigned int> partners;

    for (unsigned int ii = 0; ii < acceptedPoly.size(); ii++) {
        if (touches(newPoly, acceptedPoly[ii])) {
            partners.push_back(ii);
        }
    }

    const double *n = newPoly.normal;
    // Six unit directions. Short intersections and close-to-edge problems are
    // cleared by sliding in the fracture plane; spacing between intersections,
    // triple points and coplanar pairs are cleared by moving along the normal.
    bool inPlaneFirst = (firstRejectCode == -2 || firstRejectCode == -6);
    double directions[6][3];
    const double *axes[3];

    if (inPlaneFirst) {
        axes[0] = u;
        axes[1] = v;
        axes[2] = n;
    } else {
        axes[0] = n;
        axes[1] = u;
        axes[2] = v;
    }

    for (int a = 0; a < 3; a++) {
        for (int k = 0; k < 3; k++) {
            directions[2 * a][k]     =  axes[a][k];
            directions[2 * a + 1][k] = -axes[a][k];
        }
    }

    // Failed attempts must not count as FRAM rejections
    const RejectionReasons savedReasons = pstats.rejectionReasons;
    int attempts = 0;

    for (int step = 1; step <= wiggleMaxSteps; step++) {
        double amplitude = step * h;

        for (int d = 0; d < 6; d++) {
            double shift[3] = {amplitude * directions[d][0],
                               amplitude * directions[d][1],
                               amplitude * directions[d][2]
                              };
            attempts++;
            pstats.wiggleAttempts++;
            resetToShifted(newPoly, pristineVertices, pristineNumNodes, pristineTranslation, pristineTruncated, shift);

            // A stochastic fracture tied to a layer is clipped to it again after the shift
            if (layerTruncation(newPoly, layerIndex)) {
                continue;
            }

            if (domainTruncation(newPoly, domainSize) == 1) {
                continue; // Shift pushed the fracture out of the domain
            }

            createBoundingBox(newPoly);
            bool keptConnections = true;
            
            for (unsigned int p = 0; p < partners.size(); p++) {
                if (!touches(newPoly, acceptedPoly[partners[p]])) {
                    keptConnections = false;
                    break;
                }
            }
            
            if (!keptConnections) {
                continue; // Shift would drop an intersection the original placement had
            }
            
            if (intersectionChecking(newPoly, acceptedPoly, intPts, pstats, triplePoints) == 0) {
                pstats.rejectionReasons = savedReasons;
                pstats.wiggledPolyCount++;

                if (!isUserFracture) {
                    return 0;
                }

                WiggledFracture record;
                record.id = id;
                record.userFractureType = userFractureType;
                record.shift[0] = shift[0];
                record.shift[1] = shift[1];
                record.shift[2] = shift[2];
                record.attempts = attempts;
                pstats.wiggledFracture.push_back(record);
                std::string logString = "Wiggled user fracture " + std::to_string(id)
                                        + " (type " + std::to_string(userFractureType) + ") by {"
                                        + std::to_string(shift[0]) + ", " + std::to_string(shift[1]) + ", "
                                        + std::to_string(shift[2]) + "} after " + std::to_string(attempts)
                                        + " attempt(s); original reject code " + std::to_string(firstRejectCode) + "\n";
                logger.writeLogFile(INFO, logString);
                return 0;
            }
        }
    }

    // Nothing worked. Leave the polygon as it was before truncation.
    pstats.rejectionReasons = savedReasons;
    resetToShifted(newPoly, pristineVertices, pristineNumNodes, pristineTranslation, pristineTruncated, zeroShift);

    if (!isUserFracture) {
        return firstRejectCode;
    }

    std::string logString = "Wiggle failed for user fracture " + std::to_string(id) + " after "
                            + std::to_string(attempts) + " attempt(s)\n";
    logger.writeLogFile(INFO, logString);
    return firstRejectCode;
}
