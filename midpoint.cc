#include "midpoint.hh"
#include "regular-domain.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

#define GAMMA_BLEND

using namespace Geometry;


// Utilities

static double hermite(double t) {
  return std::pow(1 - t, 3) + 3.0 * std::pow(1 - t, 2) * t;
}

static double gammaBlend(double d) {
#ifdef GAMMA_BLEND
  return d / (2.0 * d + 1.0);
#else
  return d;
#endif
}


// Constructor

MidPoint::MidPoint(size_t n) : n_(n)
{
  outers_.resize(n_);
  inners_.resize(n_);
  multipliers_.resize(n_, 1.0);
  corners_.resize(n_);
  domain_ = std::make_unique<RegularDomain>(n_);
}


// Constraint modifications

std::pair<MidPoint::CurvePtr, MidPoint::CurvePtr>
MidPoint::interpolant(size_t i) const {
  return { outers_[i], inners_[i] };
}

void
MidPoint::setInterpolant(size_t i, const CurvePtr &outer, const CurvePtr &inner) {
  outers_[i] = outer;
  inners_[i] = inner;
}

void
MidPoint::updateCentralControlPoint() {
  Point2D center = domain_->center();
  central_cp_ = { 0, 0, 0 };
  double def;
  auto s = eval(center, &def);
  if (std::abs(def) < epsilon)
    central_cp_ = midpoint_;  // as good as anything else
  else
    central_cp_ = (midpoint_ - s) / def;
}

Geometry::Point3D
MidPoint::midpoint() const {
  return midpoint_;
}

void
MidPoint::setMidpoint(const Point3D &p) {
  midpoint_ = p;
  updateCentralControlPoint();
}

void
MidPoint::resetMidpoint() {
  midpoint_ = Point3D(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    midpoint_ += sideInterpolant(i, 0.5, 0.5);
  midpoint_ /= n_;
  updateCentralControlPoint();
}

double
MidPoint::multiplier(size_t i) const {
  return multipliers_[i];
}

void
MidPoint::setMultiplier(size_t i, double m) {
  multipliers_[i] = m;
}

const Domain *
MidPoint::domain() const {
  return domain_.get();
}


// Side- and corner interpolant evaluation

double
MidPoint::crossScaling(size_t i) const {
  return multipliers_[i];
}

Point3D
MidPoint::sideInterpolant(size_t i, double si, double di) const {
  si = std::clamp(si, 0.0, 1.0);
  di = std::max(gammaBlend(di), 0.0);
  auto p = outers_[i]->eval(si);
  auto dir = inners_[i]->eval(si) - p;
  return p + dir * crossScaling(i) * di;
}

Point3D
MidPoint::cornerInterpolant(size_t i, const Point2DVector &sds) const {
  size_t ip = (i + 1) % n_;
  double si = sds[i][0], si1 = sds[ip][0];
  return sideInterpolant(i, si, si1) + sideInterpolant(ip, si1, 1.0 - si)
    - cornerCorrection(i, 1.0 - si, si1);
}

static Vector3D rationalTwist(double u, double v, const Vector3D &f, const Vector3D &g) {
  if (std::abs(u + v) < epsilon)
    return Vector3D(0, 0, 0);
  return (f * u + g * v) / (u + v);
}

Point3D
MidPoint::cornerCorrection(size_t i, double s1, double s2) const {
  // Assumes that both s1 and s2 are 0 at the corner,
  // s1 increases towards corner (i-1), and s2 towards corner (i+1).
  s1 = std::clamp(gammaBlend(s1), 0.0, 1.0);
  s2 = std::clamp(gammaBlend(s2), 0.0, 1.0);
  return corners_[i].point
    + corners_[i].tangent1 * s1
    + corners_[i].tangent2 * s2
    + rationalTwist(s1, s2, corners_[i].twist2, corners_[i].twist1) * s1 * s2;
}

void
MidPoint::updateCorners() {
  for (size_t i = 0; i < n_; ++i) {
    size_t ip = (i + 1) % n_;
    VectorVector der;
    corners_[i].point = outers_[i]->eval(1.0, 1, der);
    corners_[i].tangent1 = -der[1];
    outers_[ip]->eval(0.0, 1, der);
    corners_[i].tangent2 = der[1];
    inners_[i]->eval(1.0, 1, der);
    corners_[i].twist1 = (-der[1] - corners_[i].tangent1) * crossScaling(i);
    inners_[ip]->eval(0.0, 1, der);
    corners_[i].twist2 = (der[1] - corners_[i].tangent2) * crossScaling(ip);
  }
}

void
MidPoint::updateDomain() {
}


// Evaluation

static Point2DVector localParameters(const DoubleVector &bc) {
  size_t n = bc.size();
  Point2DVector sds;
  for (size_t i = 0; i < n; ++i) {
    size_t im = (i + n - 1) % n;
    double s = 0.5;
    double denom = bc[im] + bc[i];
    if (denom > epsilon)
      s = bc[i] / denom;
    sds.emplace_back(s, 1 - denom);
  }
  return sds;
}

static DoubleVector deficientBlend(const Point2DVector &sds) {
  size_t n = sds.size();
  DoubleVector blf;
  for (size_t i = 0; i < n; ++i) {
    size_t ip = (i + 1) % n;
    if (sds[i][1] < epsilon && sds[ip][1] < epsilon) {
      blf.push_back(1.0);
      continue;
    }
    blf.push_back((sds[ip][1] * hermite(1.0 - sds[i][0]) * hermite(sds[i][1] ) +
                   sds[i][1]  * hermite(   sds[ip][0]  ) * hermite(sds[ip][1])) /
                  (sds[i][1] + sds[ip][1]));
  }
  return blf;
}

Point3D
MidPoint::eval(const Point2D &uv, double *deficiency) const {
  auto bc = domain_->barycentric(uv);
  auto sds = localParameters(bc);
  auto blends = deficientBlend(sds);
  Point3D p(0, 0, 0);
  for (size_t i = 0; i < n_; ++i)
    p += cornerInterpolant(i, sds) * blends[i];
  double def = 1.0 - std::accumulate(blends.begin(), blends.end(), 0.0);
  p += central_cp_ * def;
  if (deficiency)
    *deficiency = def;
  return p;
}

TriMesh
MidPoint::eval(size_t resolution) const {
  TriMesh mesh = domain_->meshTopology(resolution);
  Point2DVector uvs = domain_->parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) { return eval(uv); });
  mesh.setPoints(points);
  std::cout << sideInterpolant(3,0.5,0) << std::endl;
  std::cout << sideInterpolant(4,0,0.5) << std::endl;
  std::cout << cornerCorrection(3,0.5,0) << std::endl;
  std::cout << cornerInterpolant(3,{{0.5,1},{0.5,1},{1,0.5},{0.5,0},{0,0.5}}) << std::endl;
  return mesh;
}
