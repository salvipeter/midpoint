// Example: sphere-patch 5 0.7 0.8 30

#include <cassert>
#include <cmath>
#include <fstream>

#include "midpoint.hh"

using namespace Geometry;

class Arc : public MidPoint::Curve {
public:
  Arc(size_t n, size_t i, double z, double zo = 0.0) : z(z) {
    t0 = 2 * M_PI * i / n;
    t1 = 2 * M_PI * (i + 1) / n;
    if (zo > 0) {
      r = std::sqrt(1 - zo * zo);
      r -= zo * (z - zo) / r;
    } else
      r = std::sqrt(1 - z * z);
    p0 = Point3D(r * std::cos(t0), r * std::sin(t0), z);
    p1 = Point3D(r * std::cos(t1), r * std::sin(t1), z);
    rot = Matrix3x3::rotation((p1 - p0).normalize(), -M_PI / 6);
    //rot = Matrix3x3::identity();
  }
  Point3D eval(double u) const override {
    double t = t0 + (t1 - t0) * u;
    Point3D q(r * std::cos(t), r * std::sin(t), z);
    return p0 + rot * (q - p0);
  }
  Point3D eval(double u, size_t n, VectorVector &der) const override {
    assert(n == 1);
    der.resize(2);
    double t = t0 + (t1 - t0) * u;
    Point3D q(r * std::cos(t), r * std::sin(t), z);
    Vector3D d(-q[1], q[0], 0); 
    der[0] = p0 + rot * (q - p0);
    der[1] = rot * d * (t1 - t0);
    return der[0];
  }
private:
  double t0, t1, r, z;
  Point3D p0, p1;
  Matrix3x3 rot;
};

int main(int argc, char **argv) {
  if (argc < 4 || argc > 5) {
    std::cerr << "Usage: " << argv[0] 
      << " <# of sides> <z_out> <z_in> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 15;
  size_t n = std::atoi(argv[1]);
  double zo = std::strtod(argv[2], nullptr);
  double zi = std::strtod(argv[3], nullptr);
  if (argc == 5)
    resolution = std::atoi(argv[4]);

  MidPoint patch(n);
  for (size_t i = 0; i < n; ++i)
    patch.setInterpolant(i, 
        std::make_shared<Arc>(n, i, zo),
        std::make_shared<Arc>(n, i, zi, zo));
  patch.updateCorners();
  patch.updateDomain();
  patch.setMidpoint({0, 0, 1});

  patch.eval(resolution).writeOBJ("sphere-test.obj");

  std::ofstream f("/tmp/curve.obj");
  size_t index = 0;
  for (size_t i = 0; i < 100; ++i) {
    double u = i / 99.0;
    f << "v " << patch.interpolant(index).first->eval(u) << std::endl;
  }
  for (size_t i = 0; i < 100; ++i) {
    double u = i / 99.0;
    f << "v " << patch.interpolant(index).second->eval(u) << std::endl;
  }
  f << 'l';
  for (size_t i = 0; i < 100; ++i)
    f << ' ' << i + 1;
  f << "\nl";
  for (size_t i = 0; i < 100; ++i)
    f << ' ' << i + 101;
  f << std::endl;
}
