#include "midpoint.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

#define GAMMA_BLEND

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Geometry;


// Utilities

static double inrange(double min, double x, double max) {
  if (x < min)
    return min;
  if (x > max)
    return max;
  return x;
}

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
  // Setup Domain
  if (n_ == 4)
    domain_ = { {1,1}, {-1,1}, {-1,-1}, {1,-1} };
  else {
    double alpha = 2.0 * M_PI / n_;
    for (size_t i = 0; i < n_; ++i)
      domain_.emplace_back(std::cos(alpha * i), std::sin(alpha * i));
  }
  // Create data arrays
  outers_.resize(n_);
  inners_.resize(n_);
  corners_.resize(n_);
}


// Constraint modifications

void
MidPoint::setInterpolant(size_t i, const BSCurve &outer, const BSCurve &inner) {
  outers_[i] = outer;
  inners_[i] = inner;
}

void
MidPoint::updateCentralControlPoint() {
  Point2D center(0, 0);
  central_cp_ = { 0, 0, 0 };
  double def;
  auto s = eval(center, &def);
  if (std::abs(def) < epsilon)
    central_cp_ = midpoint_;  // as good as anything else
  else
    central_cp_ = (midpoint_ - s) / def;
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


// Side- and corner interpolant evaluation

Vector3D
MidPoint::crossDerivative(size_t i, double si) const {
  si = inrange(0, si, 1);
  auto dir = inners_[i].eval(si) - outers_[i].eval(si);
  return dir * outers_[i].basis().degree();
}

Point3D
MidPoint::sideInterpolant(size_t i, double si, double di) const {
  si = inrange(0, si, 1);
  di = std::max(gammaBlend(di), 0.0);
  return outers_[i].eval(si) + crossDerivative(i, si) * di;
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
  s1 = inrange(0, gammaBlend(s1), 1);
  s2 = inrange(0, gammaBlend(s2), 1);
  return corners_[i].point
    + corners_[i].tangent1 * s1
    + corners_[i].tangent2 * s2
    + rationalTwist(s1, s2, corners_[i].twist2, corners_[i].twist1) * s1 * s2;
}

void
MidPoint::updateCorners() {
  for (size_t i = 0; i < n_; ++i) {
    double step = 1.0e-4;
    size_t ip = (i + 1) % n_;

    VectorVector der;
    Vector3D d1, d2;
    corners_[i].point = outers_[i].eval(1.0, 1, der);
    corners_[i].tangent1 = -der[1];
    outers_[ip].eval(0.0, 1, der);
    corners_[i].tangent2 = der[1];
    d1 = crossDerivative(i, 1.0);
    d2 = crossDerivative(i, 1.0 - step);
    corners_[i].twist1 = (d2 - d1) / step;
    d1 = crossDerivative(ip, 0.0);
    d2 = crossDerivative(ip, step);
    corners_[i].twist2 = (d2 - d1) / step;
  }
}


// Parameterization

static DoubleVector wachspress(const Point2DVector &domain, const Point2D &uv) {
  size_t n = domain.size();
  Vector2DVector vectors; vectors.reserve(n);
  std::transform(domain.begin(), domain.end(), std::back_inserter(vectors),
                 [uv](const Point2D &p) { return uv - p; });

  DoubleVector areas; areas.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    const Vector2D &si = vectors[i];
    const Vector2D &si1 = vectors[(i+1)%n];
    areas.push_back((si[0] * si1[1] - si[1] * si1[0]) / 2.0);
  }

  DoubleVector l; l.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    size_t i_1 = (i + n - 1) % n, i1 = (i + 1) % n;
    double Ai = 1.0, Ai_1 = 1.0, Ai_1i = 1.0;
    for (size_t j = 0; j < n; ++j) {
      if (j == i)
        Ai_1 *= areas[j];
      else if (j == i_1)
        Ai *= areas[j];
      else {
        Ai_1 *= areas[j];
        Ai *= areas[j];
        Ai_1i *= areas[j];
      }
    }
    const Vector2D &si_1 = vectors[i_1];
    const Vector2D &si1 = vectors[i1];
    double Bi = (si_1[0] * si1[1] - si_1[1] * si1[0]) / 2.0;
    l.push_back(Ai_1 + Ai - Bi * Ai_1i);
  }

  double sum = std::accumulate(l.begin(), l.end(), 0.0);
  std::transform(l.begin(), l.end(), l.begin(), [sum](double x) { return x / sum; });
  return l;
}

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


// Evaluation

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
  auto bc = wachspress(domain_, uv);
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
  TriMesh mesh = meshTopology(resolution);
  Point2DVector uvs = parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) { return eval(uv); });
  mesh.setPoints(points);
  return mesh;
}


// Domain & Mesh

Point2DVector
MidPoint::domain() const {
  return domain_;
}

static size_t meshSize(size_t n, size_t resolution) {
  if (n == 3)
    return (resolution + 1) * (resolution + 2) / 2;
  if (n == 4)
    return (resolution + 1) * (resolution + 1);
  return 1 + n * resolution * (resolution + 1) / 2;
}

Point2DVector
MidPoint::parameters(size_t resolution) const {
  size_t size = meshSize(n_, resolution);
  Point2DVector result;
  result.reserve(size);

  if (n_ == 3) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = domain_[0] * u + domain_[2] * (1 - u);
      auto q = domain_[1] * u + domain_[2] * (1 - u);
      for (size_t k = 0; k <= j; ++k) {
        double v = j == 0 ? 1.0 : (double)k / j;
        result.push_back(p * (1 - v) + q * v);
      }
    }
  } else if (n_ == 4) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = domain_[0] * (1 - u) + domain_[1] * u;
      auto q = domain_[3] * (1 - u) + domain_[2] * u;
      for (size_t k = 0; k <= resolution; ++k) {
        double v = (double)k / resolution;
        result.push_back(p * (1 - v) + q * v);
      }
    }
  } else { // n_ > 4
    Point2D center(0.0, 0.0);
    result.push_back(center);
    for (size_t j = 1; j <= resolution; ++j) {
      double u = (double)j / (double)resolution;
      for (size_t k = 0; k < n_; ++k)
        for (size_t i = 0; i < j; ++i) {
          double v = (double)i / (double)j;
          Point2D ep = domain_[(k+n_-1)%n_] * (1.0 - v) + domain_[k] * v;
          Point2D p = center * (1.0 - u) + ep * u;
          result.push_back(p);
        }
    }
  }
  return result;
}

TriMesh
MidPoint::meshTopology(size_t resolution) const {
  TriMesh mesh;
  mesh.resizePoints(meshSize(n_, resolution));

  if (n_ == 3) {
    size_t prev = 0, current = 1;
    for (size_t i = 0; i < resolution; ++i) {
      for (size_t j = 0; j < i; ++j) {
        mesh.addTriangle(current + j, current + j + 1, prev + j);
        mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
      }
      mesh.addTriangle(current + i, current + i + 1, prev + i);
      prev = current;
      current += i + 2;
    }
  } else if (n_ == 4) {
    for (size_t i = 0; i < resolution; ++i)
      for (size_t j = 0; j < resolution; ++j) {
        size_t index = i * (resolution + 1) + j;
        mesh.addTriangle(index, index + resolution + 1, index + 1);
        mesh.addTriangle(index + 1, index + resolution + 1, index + resolution + 2);
      }
  } else { // n_ > 4
    size_t inner_start = 0, outer_vert = 1;
    for (size_t layer = 1; layer <= resolution; ++layer) {
      size_t inner_vert = inner_start, outer_start = outer_vert;
      for (size_t side = 0; side < n_; ++side) {
        size_t vert = 0;
        while(true) {
          size_t next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
          mesh.addTriangle(inner_vert, outer_vert, next_vert);
          ++outer_vert;
          if (++vert == layer)
            break;
          size_t inner_next = (side == n_ - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
          mesh.addTriangle(inner_vert, next_vert, inner_next);
          inner_vert = inner_next;
        }
      }
      inner_start = outer_start;
    }
  }
  return mesh;
}

bool
MidPoint::onEdge(size_t resolution, size_t index) const {
  if (n_ == 3) {
    if (index >= meshSize(3, resolution) - resolution - 1)
      return true;
    auto issquare = [](size_t n) {
                      size_t root = std::round(std::sqrt(n));
                      return root * root == n;
                    };
    size_t n = index * 8 + 1;
    return issquare(n) || issquare(n + 8);
  }
  if (n_ == 4) {
    return index <= resolution || index >= (resolution + 1) * resolution ||
      index % (resolution + 1) == 0 || index % (resolution + 1) == resolution;
  }
  return index >= meshSize(n_, resolution) - n_ * resolution;
}
