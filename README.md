# midpoint
Minimal library for [MidPoint patches](https://arxiv.org/abs/2002.11212).

Uses my [geometry library](https://github.com/salvipeter/libgeom/).

Based on the code in the [Transfinite library](https://github.com/salvipeter/transfinite/),
but here the ribbons are defined by pairs of curves.

The example program reads `.mp` files that have the following format:
```
<# of sides>
<1st ribbon>
...
<nth ribbon>
[midpoint]
```

If the midpoint is not given, it is put at a default position.
Each ribbon is in the format
```
<outer curve>
<inner curve>
```
where a B-spline curve is given as
```
<degree> <# of control points>
<knots>
<control points>
```

See the example file `cagd86.mp`.

Other curve types can also be used; the `sphere-patch` program shows an example.
