// Example: sphere-patch 5 0.7 1.2 30
// Example: sphere-patch 8 0.7 1.8 50

#include <cassert>
#include <cmath>
#include <fstream>

#include "midpoint.hh"

using namespace Geometry;

class Arc : public MidPoint::Curve {
public:
  Arc(size_t n, size_t i, double r, double m) : m(m) {
    double t0, t1;
    t0 = 2 * M_PI * i / n;
    t1 = 2 * M_PI * (i + 1) / n;
    p0 = Point3D(r * std::cos(t0), r * std::sin(t0), 0);
    p1 = Point3D(r * std::cos(t1), r * std::sin(t1), 0);
  }
  void setPrev(MidPoint::CurvePtr c) { VectorVector der; c->eval(1, 1, der); d0 = -der[1]; }
  void setNext(MidPoint::CurvePtr c) { VectorVector der; c->eval(0, 1, der); d1 = der[1]; }
  Point3D eval(double u) const override {
    auto p = p0 + (p1 - p0) * u;
    p[2] = std::sqrt(1 - p[0] * p[0] - p[1] * p[1]);
    if (m > 0) {
      auto d = d0 + (d1 - d0) * u;
      auto n = p.normalized();
      d -= n * (n * d);
      d *= (1 - u) * (1 - u) + 2 * (1 - u) * u * (2 * m - 1) + u * u; 
      return p + d;
    }
    return p;
  }
  Point3D eval(double u, size_t n, VectorVector &der) const override {
    assert(n == 1);
    der.resize(2);
    auto q0 = eval(u), q1 = eval(u + epsilon);
    der[0] = q0;
    der[1] = (q1 - q0) / epsilon;
    return der[0];
  }
private:
  double m;
  Point3D p0, p1;
  Vector3D d0, d1;
  MidPoint::CurvePtr prev, next;
};

int main(int argc, char **argv) {
  if (argc < 4 || argc > 5) {
    std::cerr << "Usage: " << argv[0] 
      << " <# of sides> <r> <multiplier> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 15;
  size_t n = std::atoi(argv[1]);
  double r = std::strtod(argv[2], nullptr);
  double m = std::strtod(argv[3], nullptr);
  if (argc == 5)
    resolution = std::atoi(argv[4]);

  MidPoint patch(n);
  for (size_t i = 0; i < n; ++i)
    patch.setInterpolant(i, 
        std::make_shared<Arc>(n, i, r, 0),
        std::make_shared<Arc>(n, i, r, m));
  for (size_t i = 0; i < n; ++i) {
    size_t im = (i + n - 1) % n;
    size_t ip = (i + 1) % n;
    dynamic_cast<Arc *>(patch.interpolant(i).second.get())->setPrev(patch.interpolant(im).first);
    dynamic_cast<Arc *>(patch.interpolant(i).second.get())->setNext(patch.interpolant(ip).first);
  }
  patch.updateCorners();
  patch.updateDomain();
  patch.setMidpoint({0, 0, 1});

  patch.eval(resolution).writeOBJ("sphere-test.obj");
}
