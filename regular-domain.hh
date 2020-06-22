#pragma once

#include "domain.hh"

class RegularDomain : public Domain {
public:
  RegularDomain(size_t n);
  Geometry::Point2D center() const override;
  Geometry::DoubleVector barycentric(const Geometry::Point2D &uv) const override;
  Geometry::Point2DVector parameters(size_t resolution) const override;
  Geometry::TriMesh meshTopology(size_t resolution) const override;
  bool onEdge(size_t resolution, size_t index) const override;

private:
  size_t n_;
  Geometry::Point2DVector points_;
};
