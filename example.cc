#include "midpoint.hh"

#include <algorithm>
#include <fstream>

using namespace Geometry;

class BSpline : public MidPoint::Curve {
public:
  BSpline(const BSCurve &curve) : curve(curve) { }
  Point3D eval(double u) const override {
    return curve.eval(u);
  }
  Point3D eval(double u, size_t n, VectorVector &der) const override {
    return curve.eval(u, n, der);
  }
private:
  BSCurve curve;
};

Point3D readPoint(std::istream &is) {
  Point3D p;
  is >> p[0] >> p[1] >> p[2];
  return p;
}

BSCurve readCurve(std::istream &is) {
  size_t d, n_knots, n_points;
  is >> d >> n_points;
  n_knots = n_points + d + 1;

  DoubleVector knots(n_knots);
  for (size_t i = 0; i < n_knots; ++i)
    is >> knots[i];

  PointVector cpts(n_points);
  for (size_t i = 0; i < n_points; ++i)
    is >> cpts[i];

  return { d, knots, cpts };
}

MidPoint readPatch(std::string filename) {
  std::ifstream f(filename.c_str());
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t n;
  f >> n;
  MidPoint result(n);

  for (size_t i = 0; i < n; ++i) {
    auto outer = readCurve(f); outer.normalize();
    auto inner = readCurve(f); inner.normalize();
    double multiplier; f >> multiplier;
    result.setInterpolant(i, 
        std::make_shared<BSpline>(outer),
        std::make_shared<BSpline>(inner));
    result.setMultiplier(i, multiplier);
  }

  result.updateCorners();
  result.updateDomain();

  try {
    auto mp = readPoint(f);
    result.setMidpoint(mp);
  } catch(std::ios_base::failure &) {
    result.resetMidpoint();
  }

  return result;
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <model.mp> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 15;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  auto patch = readPatch(argv[1]);
  patch.eval(resolution).writeOBJ("test.obj");
}
