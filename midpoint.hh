#pragma once

#include <geometry.hh>

// Note that
// - the interpolants have to be set before (re)setting the midpoint
// - updateCorners() has to be called every time when the corner twists have changed
// - the midpoint has to be (re)set before evaluation

// So the correct order is:
// 1. Set up interpolants
// 2. Call updateCorners()
// 3. Set up midpoint
// 4. Evaluate
// The patch can be directly re-evaluated after modification,
// except for a call to updateCorners() when the corner tangents change.

// Boundary curves are assumed to be parameterized in [0,1],
// and that they are in the correct order and orientation.

class MidPoint {
public:
  // Constructor
  MidPoint(size_t n);

  // Getters & setters
  std::pair<Geometry::BSCurve, Geometry::BSCurve> interpolant(size_t i) const;
  void setInterpolant(size_t i, const Geometry::BSCurve &outer, const Geometry::BSCurve &inner);
  Geometry::Point3D midpoint() const;
  void setMidpoint(const Geometry::Point3D &p);
  void resetMidpoint();
  double multiplier(size_t i) const;
  void setMultiplier(size_t i, double m);

  // Evaluation
  void updateCorners();
  Geometry::Point3D eval(const Geometry::Point2D &uv, double *deficiency = nullptr) const;
  Geometry::TriMesh eval(size_t resolution) const;

  // Mesh generation utilities
  Geometry::Point2DVector domain() const;
  Geometry::Point2DVector parameters(size_t resolution) const;
  Geometry::TriMesh meshTopology(size_t resolution) const;
  bool onEdge(size_t resolution, size_t index) const;

private:
  void updateCentralControlPoint();
  double crossScaling(size_t i) const;
  Geometry::Vector3D crossDerivative(size_t i, double si) const;
  Geometry::Point3D sideInterpolant(size_t i, double si, double di) const;
  Geometry::Point3D cornerCorrection(size_t i, double s1, double s2) const;
  Geometry::Point3D cornerInterpolant(size_t i, const Geometry::Point2DVector &sds) const;

  struct CornerData {
    Geometry::Point3D point;
    Geometry::Vector3D tangent1, tangent2, twist1, twist2;
  };

  size_t n_;
  Geometry::Point3D central_cp_, midpoint_;
  std::vector<Geometry::BSCurve> outers_, inners_;
  std::vector<double> multipliers_;
  std::vector<CornerData> corners_;
  Geometry::Point2DVector domain_;
};
