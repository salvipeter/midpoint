#pragma once

#include <geometry.hh>

class Domain {
public:
  virtual Geometry::Point2D center() const = 0;
  virtual Geometry::DoubleVector barycentric(const Geometry::Point2D &uv) const = 0;
  virtual Geometry::Point2DVector parameters(size_t resolution) const = 0;
  virtual Geometry::TriMesh meshTopology(size_t resolution) const = 0;
  virtual bool onEdge(size_t resolution, size_t index) const = 0;
};
