var _CSGDEBUG=false;
var AUTOCANONICALIZE=false;

/*
Todo: lots!

CSG.Polygon.prototype.expand: we should rotate the generated cylinders and spheres
such that they join nicely with the extruded plane even if a low precision is used.

Add a invertedClipTo() function so we can get rid of the inverts() here:
    b.invert();
    b.clipTo(a);
    b.invert();

-joostn
*/



// Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
// operations like union and intersection to combine 3D solids. This library
// implements CSG operations on meshes elegantly and concisely using BSP trees,
// and is meant to serve as an easily understandable implementation of the
// algorithm. All edge cases involving overlapping coplanar polygons in both
// solids are correctly handled.
// 
// Example usage:
// 
//     var cube = CSG.cube();
//     var sphere = CSG.sphere({ radius: 1.3 });
//     var polygons = cube.subtract(sphere).toPolygons();
// 
// ## Implementation Details
// 
// All CSG operations are implemented in terms of two functions, `clipTo()` and
// `invert()`, which remove parts of a BSP tree inside another BSP tree and swap
// solid and empty space, respectively. To find the union of `a` and `b`, we
// want to remove everything in `a` inside `b` and everything in `b` inside `a`,
// then combine polygons from `a` and `b` into one solid:
// 
//     a.clipTo(b);
//     b.clipTo(a);
//     a.build(b.allPolygons());
// 
// The only tricky part is handling overlapping coplanar polygons in both trees.
// The code above keeps both copies, but we need to keep them in one tree and
// remove them in the other tree. To remove them from `b` we can clip the
// inverse of `b` against `a`. The code for union now looks like this:
// 
//     a.clipTo(b);
//     b.clipTo(a);
//     b.invert();
//     b.clipTo(a);
//     b.invert();
//     a.build(b.allPolygons());
// 
// Subtraction and intersection naturally follow from set operations. If
// union is `A | B`, subtraction is `A - B = ~(~A | B)` and intersection is
// `A & B = ~(~A | ~B)` where `~` is the complement operator.
// 
// ## License
// 
// Copyright (c) 2011 Evan Wallace (http://madebyevan.com/), under the MIT license.
// Parts Copyright (c) 2012 Joost Nieuwenhuijse (joost@newhouse.nl) under the MIT license.

// # class CSG

// Holds a binary space partition tree representing a 3D solid. Two solids can
// be combined using the `union()`, `subtract()`, and `intersect()` methods.

CSG = function() {
  this.polygons = [];
};

// Construct a CSG solid from a list of `CSG.Polygon` instances.
CSG.fromPolygons = function(polygons) {
  var csg = new CSG();
  csg.polygons = polygons;
  return csg;
};

CSG.prototype = {
  toPolygons: function() {
    return this.polygons;
  },

  // Return a new CSG solid representing space in either this solid or in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.union(B)
  // 
  //     +-------+            +-------+
  //     |       |            |       |
  //     |   A   |            |       |
  //     |    +--+----+   =   |       +----+
  //     +----+--+    |       +----+       |
  //          |   B   |            |       |
  //          |       |            |       |
  //          +-------+            +-------+
  // 
  union: function(csg) {
    var a = new CSG.Tree(this.polygons);
    var b = new CSG.Tree(csg.polygons);
    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    var newpolygons = a.allPolygons().concat(b.allPolygons());
    var csg = CSG.fromPolygons(newpolygons);
    if(AUTOCANONICALIZE) csg = csg.reTesselate();
    return csg;
  },

  // Return a new CSG solid representing space in this solid but not in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.subtract(B)
  // 
  //     +-------+            +-------+
  //     |       |            |       |
  //     |   A   |            |       |
  //     |    +--+----+   =   |    +--+
  //     +----+--+    |       +----+
  //          |   B   |
  //          |       |
  //          +-------+
  // 
  subtract: function(csg) {
    var a = new CSG.Tree(this.polygons);
    var b = new CSG.Tree(csg.polygons);
    a.invert();
    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    a.addPolygons(b.allPolygons());
    a.invert();
    var csg = CSG.fromPolygons(a.allPolygons());
    if(AUTOCANONICALIZE) csg = csg.reTesselate();
    return csg;
  },

  // Return a new CSG solid representing space both this solid and in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.intersect(B)
  // 
  //     +-------+
  //     |       |
  //     |   A   |
  //     |    +--+----+   =   +--+
  //     +----+--+    |       +--+
  //          |   B   |
  //          |       |
  //          +-------+
  // 
  intersect: function(csg) {
    var a = new CSG.Tree(this.polygons);
    var b = new CSG.Tree(csg.polygons);
    a.invert();
    b.clipTo(a);
    b.invert();
    a.clipTo(b);
    b.clipTo(a);
    a.addPolygons(b.allPolygons());
    a.invert();
    var csg = CSG.fromPolygons(a.allPolygons());
    if(AUTOCANONICALIZE) csg = csg.reTesselate();
    return csg;
  },

  // Return a new CSG solid with solid and empty space switched. This solid is
  // not modified.
  inverse: function() {
    var flippedpolygons = this.polygons.map(function(p) { p.flipped(); });
    return CSG.fromPolygons(flippedpolygons);
  },
  
  // Affine transformation of CSG object. Returns a new CSG object
  transform: function(matrix4x4) {
    var newpolygons = this.polygons.map(function(p) { return p.transform(matrix4x4); } );
    return CSG.fromPolygons(newpolygons);  
  },
  
  translate: function(v) {
    return this.transform(CSG.Matrix4x4.translation(v));
  },
  
  scale: function(f) {
    return this.transform(CSG.Matrix4x4.scaling(f));
  },
  
  rotateX: function(deg) {
    return this.transform(CSG.Matrix4x4.rotationX(deg));
  },
  
  rotateY: function(deg) {
    return this.transform(CSG.Matrix4x4.rotationY(deg));
  },
  
  rotateZ: function(deg) {
    return this.transform(CSG.Matrix4x4.rotationZ(deg));
  },
  
  toStlString: function() {
    var result="solid csg.js\n";
    this.polygons.map(function(p){ result += p.toStlString(); });
    result += "endsolid csg.js\n";
    return result;
  },
  
  toString: function() {
    var result = "";
    this.polygons.map(function(p){ result += p.toString(); });
    return result;
  },

};

// Construct an axis-aligned solid cuboid. Optional parameters are `center` and
// `radius`, which default to `[0, 0, 0]` and `[1, 1, 1]`. The radius can be
// specified using a single number or a list of three numbers, one for each axis.
// 
// Example code:
// 
//     var cube = CSG.cube({
//       center: [0, 0, 0],
//       radius: 1
//     });
CSG.cube = function(options) {
  options = options || {};
  var c = new CSG.Vector(options.center || [0, 0, 0]);
  var r = !options.radius ? [1, 1, 1] : options.radius.length ?
           options.radius : [options.radius, options.radius, options.radius];
  return CSG.fromPolygons([
    [[0, 4, 6, 2], [-1, 0, 0]],
    [[1, 3, 7, 5], [+1, 0, 0]],
    [[0, 1, 5, 4], [0, -1, 0]],
    [[2, 6, 7, 3], [0, +1, 0]],
    [[0, 2, 3, 1], [0, 0, -1]],
    [[4, 5, 7, 6], [0, 0, +1]]
  ].map(function(info) {
    var normal = new CSG.Vector(info[1]);
    //var plane = new CSG.Plane(normal, 1);
    var vertices = info[0].map(function(i) {
      var pos = new CSG.Vector(
        c.x + r[0] * (2 * !!(i & 1) - 1),
        c.y + r[1] * (2 * !!(i & 2) - 1),
        c.z + r[2] * (2 * !!(i & 4) - 1)
      );
      return new CSG.Vertex(pos);
    });
    return new CSG.Polygon(vertices, null /* , plane */);
  }));
};

// Construct a solid sphere. Optional parameters are `center`, `radius`,
// `slices`, and `stacks`, which default to `[0, 0, 0]`, `1`, `16`, and `8`.
// The `slices` and `stacks` parameters control the tessellation along the
// longitude and latitude directions.
// 
// Example usage:
// 
//     var sphere = CSG.sphere({
//       center: [0, 0, 0],
//       radius: 1,
//       slices: 16,
//       stacks: 8
//     });
CSG.sphere = function(options) {
  options = options || {};
  var c = new CSG.Vector(options.center || [0, 0, 0]);
  var r = options.radius || 1;
  var slices = options.slices || 16;
  var stacks = options.stacks || 8;
  var polygons = [], vertices;
  function vertex(theta, phi) {
    theta *= Math.PI * 2;
    phi *= Math.PI;
    var dir = new CSG.Vector(
      Math.cos(theta) * Math.sin(phi),
      Math.cos(phi),
      Math.sin(theta) * Math.sin(phi)
    );
    vertices.push(new CSG.Vertex(c.plus(dir.times(r))));
  }
  for (var i = 0; i < slices; i++) {
    for (var j = 0; j < stacks; j++) {
      vertices = [];
      vertex(i / slices, j / stacks);
      if (j > 0) vertex((i + 1) / slices, j / stacks);
      if (j < stacks - 1) vertex((i + 1) / slices, (j + 1) / stacks);
      vertex(i / slices, (j + 1) / stacks);
      polygons.push(new CSG.Polygon(vertices));
    }
  }
  return CSG.fromPolygons(polygons);
};

// Construct a solid cylinder. Optional parameters are `start`, `end`,
// `radius`, and `slices`, which default to `[0, -1, 0]`, `[0, 1, 0]`, `1`, and
// `16`. The `slices` parameter controls the tessellation.
// 
// Example usage:
// 
//     var cylinder = CSG.cylinder({
//       start: [0, -1, 0],
//       end: [0, 1, 0],
//       radius: 1,
//       slices: 16
//     });
CSG.cylinder = function(options) {
  options = options || {};
  var s = new CSG.Vector(options.start || [0, -1, 0]);
  var e = new CSG.Vector(options.end || [0, 1, 0]);
  var ray = e.minus(s);
  var r = options.radius || 1;
  var slices = options.slices || 16;
  var axisZ = ray.unit(), isY = (Math.abs(axisZ.y) > 0.5);
  var axisX = new CSG.Vector(isY, !isY, 0).cross(axisZ).unit();
  var axisY = axisX.cross(axisZ).unit();
  var start = new CSG.Vertex(s);
  var end = new CSG.Vertex(e);
  var polygons = [];
  function point(stack, slice, normalBlend) {
    var angle = slice * Math.PI * 2;
    var out = axisX.times(Math.cos(angle)).plus(axisY.times(Math.sin(angle)));
    var pos = s.plus(ray.times(stack)).plus(out.times(r));
    var normal = out.times(1 - Math.abs(normalBlend)).plus(axisZ.times(normalBlend));
    return new CSG.Vertex(pos);
  }
  for (var i = 0; i < slices; i++) {
    var t0 = i / slices, t1 = (i + 1) / slices;
    polygons.push(new CSG.Polygon([start, point(0, t0, -1), point(0, t1, -1)]));
    polygons.push(new CSG.Polygon([point(0, t1, 0), point(0, t0, 0), point(1, t0, 0), point(1, t1, 0)]));
    polygons.push(new CSG.Polygon([end, point(1, t1, 1), point(1, t0, 1)]));
  }
  return CSG.fromPolygons(polygons);
};

// # class Vector

// Represents a 3D vector.
// 
// Example usage:
// 
//     new CSG.Vector(1, 2, 3);
//     new CSG.Vector([1, 2, 3]);
//     new CSG.Vector({ x: 1, y: 2, z: 3 });

CSG.Vector = function(x, y, z) {
  if (arguments.length == 3) {
    this.x = x;
    this.y = y;
    this.z = z;
  } else if ('x' in x) {
    this.x = x.x;
    this.y = x.y;
    this.z = x.z;
  } else {
    this.x = x[0];
    this.y = x[1];
    this.z = x[2];
  }
};

CSG.Vector.prototype = {
  clone: function() {
    return new CSG.Vector(this.x, this.y, this.z);
  },

  negated: function() {
    return new CSG.Vector(-this.x, -this.y, -this.z);
  },

  plus: function(a) {
    return new CSG.Vector(this.x + a.x, this.y + a.y, this.z + a.z);
  },

  minus: function(a) {
    return new CSG.Vector(this.x - a.x, this.y - a.y, this.z - a.z);
  },

  times: function(a) {
    return new CSG.Vector(this.x * a, this.y * a, this.z * a);
  },

  dividedBy: function(a) {
    return new CSG.Vector(this.x / a, this.y / a, this.z / a);
  },

  dot: function(a) {
    return this.x * a.x + this.y * a.y + this.z * a.z;
  },

  lerp: function(a, t) {
    return this.plus(a.minus(this).times(t));
  },

  lengthSquared: function() {
    return this.dot(this);
  },

  length: function() {
    return Math.sqrt(this.lengthSquared());
  },

  unit: function() {
    return this.dividedBy(this.length());
  },

  cross: function(a) {
    return new CSG.Vector(
      this.y * a.z - this.z * a.y,
      this.z * a.x - this.x * a.z,
      this.x * a.y - this.y * a.x
    );
  },
  
  distanceTo: function(a) {
    return this.minus(a).length();
  },

  distanceToSquared: function(a) {
    return this.minus(a).lengthSquared();
  },

  equals: function(a) {
    return (this.x == a.x) && (this.y == a.y) && (this.z == a.z);
  },
  
  // Right multiply by a 4x4 matrix (the vector is interpreted as a row vector)
  // Returns a new CSG.Vector
  multiply4x4: function(matrix4x4) {
    return matrix4x4.rightMultiply1x3Vector(this);
  },
  
  toStlString: function() {
    return this.x+" "+this.y+" "+this.z;
  },
  
  toString: function() {
    return "("+this.x+", "+this.y+", "+this.z+")";
  },
  
};

// # class Vertex

// Represents a vertex of a polygon. Use your own vertex class instead of this
// one to provide additional features like texture coordinates and vertex
// colors. Custom vertex classes need to provide a `pos` property
// `flipped()`, and `interpolate()` methods that behave analogous to the ones
// defined by `CSG.Vertex`.

CSG.Vertex = function(pos) {
  this.pos = pos;
};

CSG.Vertex.prototype = {
  // Return a vertex with all orientation-specific data (e.g. vertex normal) flipped. Called when the
  // orientation of a polygon is flipped.
  flipped: function() {
    return this;
  },

  getTag: function() {
    var result = this.tag;
    if(!result)
    {
      result = CSG.getTag();
      this.tag = result;
    }
    return result;
  },

  // Create a new vertex between this vertex and `other` by linearly
  // interpolating all properties using a parameter of `t`. Subclasses should
  // override this to interpolate additional properties.
  interpolate: function(other, t) {
    var newpos = this.pos.lerp(other.pos, t);
    return new CSG.Vertex(newpos);
  },
  
  // Affine transformation of vertex. Returns a new CSG.Vertex
  transform: function(matrix4x4) {
    var newpos = this.pos.multiply4x4(matrix4x4);
    return new CSG.Vertex(newpos);  
  },
  
  toStlString: function() {
    return "vertex "+this.pos.toStlString()+"\n";
  },
  
  toString: function() {
    return this.pos.toString();
  },
};

// # class Plane

// Represents a plane in 3D space.

CSG.Plane = function(normal, w) {
  this.normal = normal;
  this.w = w;
};

// `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
// point is on the plane.
CSG.Plane.EPSILON = 1e-5;

CSG.Plane.fromPoints = function(a, b, c) {
  var n = b.minus(a).cross(c.minus(a)).unit();
  return new CSG.Plane(n, n.dot(a));
};

CSG.Plane.prototype = {
  flipped: function() {
    return new CSG.Plane(this.normal.negated(), -this.w);
  },
  
  getTag: function() {
    var result = this.tag;
    if(!result)
    {
      result = CSG.getTag();
      this.tag = result;
    }
    return result;
  },

  equals: function(n) {
    return this.normal.equals(n.normal) && this.w == n.w;
  },
  
  transform: function(matrix4x4) {
    var origin = new CSG.Vector(0,0,0);
    var pointOnPlane = this.normal.times(this.w);
    var neworigin = origin.multiply4x4(matrix4x4);
    var neworiginPlusNormal = this.normal.multiply4x4(matrix4x4);
    var newnormal = neworiginPlusNormal.minus(neworigin);
    var newpointOnPlane = pointOnPlane.multiply4x4(matrix4x4);
    var neww = newnormal.dot(newpointOnPlane);
    return new CSG.Plane(newnormal, neww);
  },

  // Returns object:
  // .type:
  //   0: coplanar-front
  //   1: coplanar-back
  //   2: front
  //   3: back
  //   4: spanning
  // In case the polygon is spanning, returns:
  // .front: a CSG.Polygon of the front part 
  // .back: a CSG.Polygon of the back part 
  splitPolygonNew: function(polygon) {
    var result = {
      type: null,
      front: null,
      back: null,
    };    
    // cache in local vars (speedup):
    var planenormal = this.normal;      
    var vertices = polygon.vertices;
    var numvertices = vertices.length;
    if(polygon.plane.equals(this))
    {
      result.type = 0;
    }
    else
    {   
      var hasfront = false;
      var hasback = false;
      var vertexIsBack = [];
      var thisw = this.w;
      var EPS = CSG.Plane.EPSILON;
      var MINEPS = -EPS;
      for (var i = 0; i < numvertices; i++) {
        var t = planenormal.dot(vertices[i].pos) - thisw;
        var isback = (t < 0);
        vertexIsBack.push(isback);
        if(t > EPS) hasfront = true; 
        if(t < MINEPS) hasback = true; 
      }
      if( (!hasfront) && (!hasback) )
      {
        // all points coplanar
        var t = planenormal.dot(polygon.plane.normal);
        result.type = (t >= 0)? 0:1;
      }
      else if(!hasback)
      {
        result.type = 2;
      }
      else if(!hasfront)
      {
        result.type = 3;
      }
      else
      {
        // spanning
        result.type = 4;
        var frontvertices = [], backvertices = [];
        var isback = vertexIsBack[0];
        for(var vertexindex = 0; vertexindex < numvertices; vertexindex++)
        {
          var vertex = vertices[vertexindex];
          var nextvertexindex = vertexindex + 1;
          if(nextvertexindex >= numvertices) nextvertexindex = 0;
          var nextisback = vertexIsBack[nextvertexindex];
          if(isback == nextisback)
          {
            // line segment is on one side of the plane:
            if(isback)
            {
              backvertices.push(vertex);
            }
            else
            {
              frontvertices.push(vertex);
            }          
          }
          else
          {
            // line segment intersects plane:
            var point = vertex.pos;
            var nextpoint = vertices[nextvertexindex].pos;
            var line = CSG.Line3D.fromPoints(point, nextpoint);
            var intersectionpoint =  this.intersectWithLine(line);
            var intersectionvertex = new CSG.Vertex(intersectionpoint);
            if(isback)
            {
              backvertices.push(vertex);
              backvertices.push(intersectionvertex);
              frontvertices.push(intersectionvertex);
            }
            else
            {
              frontvertices.push(vertex);
              frontvertices.push(intersectionvertex);
              backvertices.push(intersectionvertex);
            }          
          }
          isback = nextisback;
        }  // for vertexindex

        // remove duplicate vertices:
        var EPS_SQUARED = CSG.Plane.EPSILON * CSG.Plane.EPSILON;  
        if(backvertices.length >= 3)
        {
          var prevvertex = backvertices[backvertices.length - 1];
          for(var vertexindex = 0; vertexindex < backvertices.length; vertexindex++)
          {
            var vertex = backvertices[vertexindex];
            if(vertex.pos.distanceToSquared(prevvertex.pos) < EPS_SQUARED)
            {
              backvertices.splice(vertexindex,1);
              vertexindex--;
            }
            prevvertex = vertex;
          }        
        }
        if(frontvertices.length >= 3)
        {
          var prevvertex = frontvertices[frontvertices.length - 1];
          for(var vertexindex = 0; vertexindex < frontvertices.length; vertexindex++)
          {
            var vertex = frontvertices[vertexindex];
            if(vertex.pos.distanceToSquared(prevvertex.pos) < EPS_SQUARED)
            {
              frontvertices.splice(vertexindex,1);
              vertexindex--;
            }
            prevvertex = vertex;
          }        
        }
  
        if (frontvertices.length >= 3)
        {
          result.front = new CSG.Polygon(frontvertices, polygon.shared, polygon.plane); 
        }
        if (backvertices.length >= 3)
        {
          result.back = new CSG.Polygon(backvertices, polygon.shared, polygon.plane); 
        }
      }
    }
    return result;
  },

  // returns CSG.Point3D
  intersectWithLine: function(line3d) {
    return line3d.intersectWithPlane(this);
  },

  // intersection of two planes
  intersectWithPlane: function(plane) {
    return CSG.Line3D.fromPlanes(this, plane);
  },

  signedDistanceToPoint: function(point) {
    var t = this.normal.dot(point) - this.w;
    return t;
  },

  toString: function() {
    return "[normal: "+this.normal.toString()+", w: "+this.w+"]";
  },
};


// # class Polygon

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. They do not have to be `CSG.Vertex`
// instances but they must behave similarly (duck typing can be used for
// customization).
// 
// Each convex polygon has a `shared` property, which is shared between all
// polygons that are clones of each other or were split from the same polygon.
// This can be used to define per-polygon properties (such as surface color).
// 
// The plane of the polygon is calculated from the vertex coordinates
// To avoid unnecessary recalculation, the plane can alternatively be
// passed as the third argument 
CSG.Polygon = function(vertices, shared, plane) {
  this.vertices = vertices;
  this.shared = shared;
  var numvertices = vertices.length;

  if(arguments.length >= 3)
  {
    this.plane = plane;
  }
  else
  {
    this.plane = CSG.Plane.fromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos);
  }

  if(_CSGDEBUG)
  {
    this.checkIfConvex();
  }
};

CSG.Polygon.prototype = {
  // check whether the polygon is convex (it should be, otherwise we will get unexpected results)
  checkIfConvex: function() {
    if(! CSG.Polygon.verticesConvex(this.vertices, this.plane.normal))
    {
      throw new Error("Not convex!");
    }
  },

  flipped: function() {
    var newvertices = this.vertices.map(function(v) { return v.flipped(); });
    newvertices.reverse();
    var newplane = this.plane.flipped();
    return new CSG.Polygon(newvertices, this.shared, newplane);
  },
  
  // Affine transformation of polygon. Returns a new CSG.Polygon
  transform: function(matrix4x4) {
    var newvertices = this.vertices.map(function(v) { return v.transform(matrix4x4); } );
    var newplane = this.plane.transform(matrix4x4);
    return new CSG.Polygon(newvertices, this.shared, newplane);
  },
  
  toStlString: function() {
    var result="";
    if(this.vertices.length >= 3) // should be!
    {
      // STL requires triangular polygons. If our polygon has more vertices, create
      // multiple triangles:
      var firstVertexStl = this.vertices[0].toStlString();
      for(var i=0; i < this.vertices.length-2; i++)
      {
        result += "facet normal "+this.plane.normal.toStlString()+"\nouter loop\n";
        result += firstVertexStl;
        result += this.vertices[i+1].toStlString();
        result += this.vertices[i+2].toStlString();
        result += "endloop\nendfacet\n";    
      } 
    }
    return result;
  },
  
  toString: function() {
    var result = "Polygon plane: "+this.plane.toString()+"\n";
    this.vertices.map(function(vertex) {
      result += "  "+vertex.toString()+"\n";
    });
    return result;
  },  
};

CSG.Polygon.verticesConvex = function(vertices, planenormal) {
  var numvertices = vertices.length;
  if(numvertices > 2)
  {
    var prevprevpos=vertices[numvertices-2].pos;
    var prevpos=vertices[numvertices-1].pos;
    for(var i=0; i < numvertices; i++)
    {
      var pos=vertices[i].pos;
      if(!CSG.Polygon.isConvexPoint(prevprevpos, prevpos, pos, planenormal))
      {
        return false;
      }
      prevprevpos=prevpos;
      prevpos=pos;
    }
  }
  return true;
};

CSG.Polygon.verticesStrictlyConvex = function(vertices, planenormal) {
  var numvertices = vertices.length;
  if(numvertices > 2)
  {
    var prevprevpos=vertices[numvertices-2].pos;
    var prevpos=vertices[numvertices-1].pos;
    for(var i=0; i < numvertices; i++)
    {
      var pos=vertices[i].pos;
      if(!CSG.Polygon.isStrictlyConvexPoint(prevprevpos, prevpos, pos, planenormal))
      {
        return false;
      }
      prevprevpos=prevpos;
      prevpos=pos;
    }
  }
  return true;
};

// Create a polygon from the given points
CSG.Polygon.createFromPoints = function(points, shared, plane) {
  var normal;
  if(arguments.length < 3)
  {
    // initially set a dummy vertex normal:
    normal = new CSG.Vector(0, 0, 0);
  }
  else
  {
    normal = plane.normal;
  }
  var vertices = [];
  points.map( function(p) {
    var vec = new CSG.Vector(p);
    var vertex = new CSG.Vertex(vec);
    vertices.push(vertex); 
  });            
  var polygon;
  if(arguments.length < 3)
  {
    polygon = new CSG.Polygon(vertices, shared);
  }
  else
  {
    polygon = new CSG.Polygon(vertices, shared, plane);
  }
  return polygon;
};

// calculate whether three points form a convex corner 
//  prevpoint, point, nextpoint: the 3 coordinates (CSG.Vector instances)
//  normal: the normal vector of the plane
CSG.Polygon.isConvexPoint = function(prevpoint, point, nextpoint, normal) {
  var crossproduct=point.minus(prevpoint).cross(nextpoint.minus(point));
  var crossdotnormal=crossproduct.dot(normal);
  return (crossdotnormal >= 0);
};

CSG.Polygon.isStrictlyConvexPoint = function(prevpoint, point, nextpoint, normal) {
  var crossproduct=point.minus(prevpoint).cross(nextpoint.minus(point));
  var crossdotnormal=crossproduct.dot(normal);
  return (crossdotnormal >= 1e-5);
};

// # class PolygonTreeNode

// This class manages hierarchical splits of polygons
// At the top is a root node which doesn hold a polygon, only child PolygonTreeNodes
// Below that are zero or more 'top' nodes; each holds a polygon. The polygons can be in different planes 
// splitByPlane() splits a node by a plane. If the plane intersects the polygon, two new child nodes
// are created holding the splitted polygon.
// getPolygons() retrieves the polygon from the tree. If for PolygonTreeNode the polygon is split but 
// the two split parts (child nodes) are still intact, then the unsplit polygon is returned.
// This ensures that we can safely split a polygon into many fragments. If the fragments are untouched,
//  getPolygons() will return the original unsplit polygon instead of the fragments.
// remove() removes a polygon from the tree. Once a polygon is removed, the parent polygons are invalidated 
// since they are no longer intact. 

// constructor creates the root node:
CSG.PolygonTreeNode = function() {
  this.parent = null;
  this.children = [];
  this.polygon = null;
};

CSG.PolygonTreeNode.prototype = {
  // fill the tree with polygons. Should be called on the root node only; child nodes must
  // always be a derivate (split) of the parent node.
  addPolygons: function(polygons) {
    if(!this.isRootNode()) throw new Error("Assertion failed");  // new polygons can only be added to root node; children can only be splitted polygons
    var _this = this;
    polygons.map(function(polygon) {
      _this.addChild(polygon);
    });
  },
  
  // remove a node
  // - the siblings become toplevel nodes
  // - the parent is removed recursively
  remove: function() {
    if(this.isRootNode()) throw new Error("Assertion failed");  // can't remove root node
    if(this.children.length) throw new Error("Assertion failed"); // we shouldn't remove nodes with children
    
    // remove ourselves from the parent's children list:
    var parentschildren = this.parent.children;
    var i = parentschildren.indexOf(this);
    if(i < 0) throw new Error("Assertion failed");
    parentschildren.splice(i,1);
    
    // invalidate the parent's polygon, and of all parents above it:
    this.parent.recursivelyInvalidatePolygon();
  },

  isRootNode: function() {
    return !this.parent;
  },  

  // invert all polygons in the tree. Call on the root node
  invert: function() {
    if(!this.isRootNode()) throw new Error("Assertion failed");  // can only call this on the root node
    this.invertSub();
  },

  getPolygon: function () {
    if(!this.polygon) throw new Error("Assertion failed");  // doesn't have a polygon, which means that it has been broken down
    return this.polygon;
  },

  getPolygons: function (result) {
    if(this.polygon)
    {
      // the polygon hasn't been broken yet. We can ignore the children and return our polygon:
      result.push(this.polygon);
    }
    else
    {
      // our polygon has been split up and broken, so gather all subpolygons from the children:
      var childpolygons = [];
      this.children.map(function(child) {
        child.getPolygons(childpolygons);
      });
      childpolygons.map(function(p) {
        result.push(p);
      });
    }
  },

  // split the node by a plane; add the resulting nodes to the frontnodes and backnodes array  
  // If the plane doesn't intersect the polygon, the 'this' object is added to one of the arrays
  // If the plane does intersect the polygon, two new child nodes are created for the front and back fragments,
  //  and added to both arrays. 
  splitByPlane: function(plane, coplanarfrontnodes, coplanarbacknodes, frontnodes, backnodes) {
    var children = this.children;
    var numchildren = children.length; 
    if(numchildren > 0)
    {
      // if we have children, split the children
      for(var i = 0; i < numchildren; i++)
      {
        children[i].splitByPlane(plane, coplanarfrontnodes, coplanarbacknodes, frontnodes, backnodes);
      }
    }
    else
    {
      // no children. Split the polygon:
      if(this.polygon)
      {
        var splitresult = plane.splitPolygonNew(this.polygon);
        switch(splitresult.type)
        {
          case 0:   // coplanar front:
            coplanarfrontnodes.push(this);
            break;
            
          case 1:   // coplanar back:
            coplanarbacknodes.push(this);
            break;
            
          case 2:   // front:
            frontnodes.push(this);
            break;
            
          case 3:   // back:
            backnodes.push(this);
            break;
            
          case 4:  // spanning:
            if(splitresult.front)
            {
              var frontnode = this.addChild(splitresult.front);
              frontnodes.push(frontnode);
            }
            if(splitresult.back)
            {
              var backnode = this.addChild(splitresult.back);
              backnodes.push(backnode);
            }
            break;
        }
        

/*
        var coplanarfrontpolygons = [], coplanarbackpolygons = [];
        var frontpolygons = [], backpolygons = [];
        plane.splitPolygon(this.polygon, coplanarfrontpolygons, coplanarbackpolygons, frontpolygons, backpolygons);
        var numpolygons=coplanarfrontpolygons.length + coplanarbackpolygons.length + frontpolygons.length + backpolygons.length; 
        if(numpolygons > 2) throw new Error("Assertion failed");
        if(numpolygons == 1)
        {
          // the polygon was not split:
          if(frontpolygons.length == 1)
          {
            frontnodes.push(this);
          }
          else if(backpolygons.length == 1)
          {
            backnodes.push(this);
          }
          else if(coplanarfrontpolygons.length == 1)
          {
            coplanarfrontnodes.push(this);
          }
          else // if(coplanarbackpolygons.length == 1)
          {
            coplanarbacknodes.push(this);
          }
        }
        else if(numpolygons == 2)
        {
          // the polygon was split. Split our node by creating two children:
          if(frontpolygons.length != 1) throw new Error("Assertion failed");
          if(backpolygons.length != 1) throw new Error("Assertion failed");
          var frontnode = this.addChild(frontpolygons[0]);
          var backnode = this.addChild(backpolygons[0]);
          frontnodes.push(frontnode);
          backnodes.push(backnode);          
        }
*/        
      }      
    }
  },
  
 
  // PRIVATE methods from here:

  // add child to a node
  // this should be called whenever the polygon is split
  // a child should be created for every fragment of the split polygon 
  // returns the newly created child
  addChild: function(polygon) {
    var newchild = new CSG.PolygonTreeNode();
    newchild.parent = this;
    newchild.polygon = polygon;
    this.children.push(newchild);
    return newchild;
  },

  invertSub: function() {
    if(this.polygon)
    {
      this.polygon = this.polygon.flipped();
    }
    this.children.map(function(child) {
      child.invertSub();
    });
  },
  
  recursivelyInvalidatePolygon: function() {
    this.polygon = null;
    if(this.parent)
    {
      this.parent.recursivelyInvalidatePolygon();
    }
  },
  
};



// # class Tree
// This is the root of a BSP tree
// We are using this separate class for the root of the tree, to hold the PolygonTreeNode root
// The actual tree is kept in this.rootnode
CSG.Tree = function(polygons) {
  this.polygonTree = new CSG.PolygonTreeNode();
  this.rootnode = new CSG.Node();
  if (polygons) this.addPolygons(polygons);
};

CSG.Tree.prototype = {
  invert: function() {
    this.polygonTree.invert();
    this.rootnode.invert();
  },
  
  // Remove all polygons in this BSP tree that are inside the other BSP tree
  // `tree`.
  clipTo: function(tree) {
    this.rootnode.clipTo(tree);
  },

  allPolygons: function() {
    var result = [];
    this.polygonTree.getPolygons(result);
    return result;
  },

  addPolygons: function(polygons) {
    var _this = this;
    polygons.map(function(p) {
      _this.addPolygon(p);
    });
  },

  addPolygon: function(polygon) {
    var polygontreenode=this.polygonTree.addChild(polygon);
    this.rootnode.addPolygonTreeNode(polygontreenode);
  },  
};

// # class Node

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along.
// Polygons are not stored directly in the tree, but in PolygonTreeNodes, stored in
// this.polygontreenodes. Those PolygonTreeNodes are children of the owning
// CSG.Tree.polygonTree
// This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.

CSG.Node = function() {
  this.plane = null;
  this.front = null;
  this.back = null;
  this.polygontreenodes = [];
};

CSG.Node.prototype = {
  // Convert solid space to empty space and empty space to solid space.
  invert: function() {
    this.plane = this.plane.flipped();
    if (this.front) this.front.invert();
    if (this.back) this.back.invert();
    var temp = this.front;
    this.front = this.back;
    this.back = temp;
  },

  // clip polygontreenodes to our plane
  // returns a new array of PolygonTreeNodes with the clipped polygons
  // also calls remove() for all clipped PolygonTreeNodes, so the polygon tree is modified!
  clipPolygons: function(polygontreenodes) {
    var clippednodes;
    if(this.plane)
    {
      var backnodes = [];
      var frontnodes = [];
      var plane = this.plane;
      polygontreenodes.map(function(node) {
        node.splitByPlane(plane, frontnodes, backnodes, frontnodes, backnodes);
      });
      var clippedfrontnodes = [];
      var clippedbacknodes = [];
      if(this.front)
      {
        clippedfrontnodes = this.front.clipPolygons(frontnodes);
      }
      else
      {
        clippedfrontnodes = frontnodes;
      }
      if(this.back)
      {
        clippedbacknodes = this.back.clipPolygons(backnodes);
      }
      else
      {
        // there's nothing behind this plane. Delete the nodes behind this plane:
        backnodes.map( function(node) {
          node.remove();
        });
      }
      clippednodes = clippedfrontnodes.concat(clippedbacknodes);
    }
    else
    {
      clippednodes=polygontreenodes.slice();
    }
    return clippednodes;
  },

  // Remove all polygons in this BSP tree that are inside the other BSP tree
  // `tree`.
  clipTo: function(tree) {
    var origpolygontreenodes = this.polygontreenodes;
    if(origpolygontreenodes.length > 0)
    {
      var clippedtreenodes = tree.rootnode.clipPolygons(origpolygontreenodes);
      this.polygontreenodes = clippedtreenodes;
    }
    if (this.front) this.front.clipTo(tree);
    if (this.back) this.back.clipTo(tree);
  },
  
  addPolygonTreeNode: function(polygontreenode) {
    if(!this.plane)
    {
      this.plane = polygontreenode.getPolygon().plane;
    }
    var frontnodes = [];
    var backnodes = [];
    polygontreenode.splitByPlane(this.plane, this.polygontreenodes, this.polygontreenodes, frontnodes, backnodes);
    if(frontnodes.length > 0)
    {
      if (!this.front) this.front = new CSG.Node();
      this.front.addPolygonTreeNode(frontnodes[0]);
    }
    if(backnodes.length > 0)
    {
      if (!this.back) this.back = new CSG.Node();
      this.back.addPolygonTreeNode(backnodes[0]);
    }
  },
};

//////////

// # class Matrix4x4:
// Represents a 4x4 matrix. Elements are specified in row order
CSG.Matrix4x4 = function(elements) {
  if (arguments.length >= 1) {
    this.elements=elements;
  }
  else
  {
    // if no arguments passed: create unity matrix  
    this.elements=[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
  }
}

CSG.Matrix4x4.prototype = {
  plus: function(m) {
    var r=[];
    for(var i=0; i < 16; i++)
    {
      r[i]=this.elements[i]+m.elements[i];
    }
    return new CSG.Matrix4x4(r);
  },
  
  minus: function(m) {
    var r=[];
    for(var i=0; i < 16; i++)
    {
      r[i]=this.elements[i]-m.elements[i];
    }
    return new CSG.Matrix4x4(r);
  },

  // right multiply by another 4x4 matrix:
  multiply: function(m) {
    // cache elements in local variables, for speedup:
    var this0=this.elements[0];
    var this1=this.elements[1];
    var this2=this.elements[2];
    var this3=this.elements[3];
    var this4=this.elements[4];
    var this5=this.elements[5];
    var this6=this.elements[6];
    var this7=this.elements[7];
    var this8=this.elements[8];
    var this9=this.elements[9];
    var this10=this.elements[10];
    var this11=this.elements[11];
    var this12=this.elements[12];
    var this13=this.elements[13];
    var this14=this.elements[14];
    var this15=this.elements[15];
    var m0=m.elements[0];
    var m1=m.elements[1];
    var m2=m.elements[2];
    var m3=m.elements[3];
    var m4=m.elements[4];
    var m5=m.elements[5];
    var m6=m.elements[6];
    var m7=m.elements[7];
    var m8=m.elements[8];
    var m9=m.elements[9];
    var m10=m.elements[10];
    var m11=m.elements[11];
    var m12=m.elements[12];
    var m13=m.elements[13];
    var m14=m.elements[14];
    var m15=m.elements[15];
    
    var result=[];
    result[0] = this0*m0 + this1*m4 + this2*m8 + this3*m12;
    result[1] = this0*m1 + this1*m5 + this2*m9 + this3*m13;
    result[2] = this0*m2 + this1*m6 + this2*m10 + this3*m14;
    result[3] = this0*m3 + this1*m7 + this2*m11 + this3*m15;
    result[4] = this4*m0 + this5*m4 + this6*m8 + this7*m12;
    result[5] = this4*m1 + this5*m5 + this6*m9 + this7*m13;
    result[6] = this4*m2 + this5*m6 + this6*m10 + this7*m14;
    result[7] = this4*m3 + this5*m7 + this6*m11 + this7*m15;
    result[8] = this8*m0 + this9*m4 + this10*m8 + this11*m12;
    result[9] = this8*m1 + this9*m5 + this10*m9 + this11*m13;
    result[10] = this8*m2 + this9*m6 + this10*m10 + this11*m14;
    result[11] = this8*m3 + this9*m7 + this10*m11 + this11*m15;
    result[12] = this12*m0 + this13*m4 + this14*m8 + this15*m12;
    result[13] = this12*m1 + this13*m5 + this14*m9 + this15*m13;
    result[14] = this12*m2 + this13*m6 + this14*m10 + this15*m14;
    result[15] = this12*m3 + this13*m7 + this14*m11 + this15*m15;
    return new CSG.Matrix4x4(result);
  },
  
  clone: function() {
    var elements = this.elements.map(function(p) { return p; }); 
    return new CSG.Matrix4x4(elements);
  },
  
  // Multiply a CSG.Vector (interpreted as 1 row, 3 column) by this matrix 
  // Fourth element is taken as 1
  rightMultiply1x3Vector: function(v) {
    var v0 = v.x;
    var v1 = v.y;
    var v2 = v.z;
    var v3 = 1;    
    var x = v0*this.elements[0] + v1*this.elements[1] + v2*this.elements[2] + v3*this.elements[3];    
    var y = v0*this.elements[4] + v1*this.elements[5] + v2*this.elements[6] + v3*this.elements[7];    
    var z = v0*this.elements[8] + v1*this.elements[9] + v2*this.elements[10] + v3*this.elements[11];    
    var w = v0*this.elements[12] + v1*this.elements[13] + v2*this.elements[14] + v3*this.elements[15];
    // scale such that fourth element becomes 1:
    if(w != 1)
    {
      var invw=1.0/w;
      x *= invw;
      y *= invw;
      z *= invw;
    }
    return new CSG.Vector(x,y,z);       
  },
  
  // Multiply a CSG.Vector2D (interpreted as 1 row, 2 column) by this matrix 
  // Fourth element is taken as 1
  rightMultiply1x2Vector: function(v) {
    var v0 = v.x;
    var v1 = v.y;
    var v2 = 0;
    var v3 = 1;    
    var x = v0*this.elements[0] + v1*this.elements[1] + v2*this.elements[2] + v3*this.elements[3];    
    var y = v0*this.elements[4] + v1*this.elements[5] + v2*this.elements[6] + v3*this.elements[7];    
    var z = v0*this.elements[8] + v1*this.elements[9] + v2*this.elements[10] + v3*this.elements[11];    
    var w = v0*this.elements[12] + v1*this.elements[13] + v2*this.elements[14] + v3*this.elements[15];
    // scale such that fourth element becomes 1:
    if(w != 1)
    {
      var invw=1.0/w;
      x *= invw;
      y *= invw;
      z *= invw;
    }
    return new CSG.Vector2D(x,y);       
  },
};

// return the unity matrix
CSG.Matrix4x4.unity = function() {
  return new CSG.Matrix4x4(); 
};

// Create a rotation matrix for rotating around the x axis
CSG.Matrix4x4.rotationX = function(degrees) {
  var radians = degrees * Math.PI * (1.0/180.0);
  var cos = Math.cos(radians);
  var sin = Math.sin(radians);
  var els = [
    1, 0, 0, 0,
    0, cos, -sin, 0,
    0, sin, cos, 0,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

// Create a rotation matrix for rotating around the y axis
CSG.Matrix4x4.rotationY = function(degrees) {
  var radians = degrees * Math.PI * (1.0/180.0);
  var cos = Math.cos(radians);
  var sin = Math.sin(radians);
  var els = [
    cos, 0, sin, 0,
    0, 1, 0, 0,
    -sin, 0, cos, 0,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

// Create a rotation matrix for rotating around the z axis
CSG.Matrix4x4.rotationZ = function(degrees) {
  var radians = degrees * Math.PI * (1.0/180.0);
  var cos = Math.cos(radians);
  var sin = Math.sin(radians);
  var els = [
    cos, -sin, 0, 0,
    sin, cos, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

// Create an affine matrix for translation:
CSG.Matrix4x4.translation = function(v) {
  // parse as CSG.Vector, so we can pass an array or a CSG.Vector
  var vec = new CSG.Vector(v);
  var els = [
    1, 0, 0, vec.x,
    0, 1, 0, vec.y,
    0, 0, 1, vec.z,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

// Create an affine matrix for scaling:
CSG.Matrix4x4.scaling = function(v) {
  // parse as CSG.Vector, so we can pass an array or a CSG.Vector
  var vec = new CSG.Vector(v);
  var els = [
    vec.x, 0, 0, 0,
    0, vec.y, 0, 0,
    0, 0, vec.z, 0,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

///////////////////////////////////////////////////

// # class Vector2D:
// Represents a 2 element vector
CSG.Vector2D = function(x, y) {
  if (arguments.length == 2) {
    this.x = x;
    this.y = y;
  } else if ('x' in x) {
    this.x = x.x;
    this.y = x.y;
  } else {
    this.x = x[0];
    this.y = x[1];
  }
};

CSG.Vector2D.prototype = {
  // extend to a 3D vector by adding a z coordinate:
  toVector3D: function(z) {
    return new CSG.Vector2D(this.x, this.y, z);
  },
  
  equals: function(a) {
    return (this.x == a.x) && (this.y == a.y);
  },
  
  clone: function() {
    return new CSG.Vector2D(this.x, this.y);
  },

  negated: function() {
    return new CSG.Vector2D(-this.x, -this.y);
  },

  plus: function(a) {
    return new CSG.Vector2D(this.x + a.x, this.y + a.y);
  },

  minus: function(a) {
    return new CSG.Vector2D(this.x - a.x, this.y - a.y);
  },

  times: function(a) {
    return new CSG.Vector2D(this.x * a, this.y * a);
  },

  dividedBy: function(a) {
    return new CSG.Vector2D(this.x / a, this.y / a);
  },

  dot: function(a) {
    return this.x * a.x + this.y * a.y;
  },

  lerp: function(a, t) {
    return this.plus(a.minus(this).times(t));
  },

  length: function() {
    return Math.sqrt(this.dot(this));
  },

  distanceTo: function(a) {
    return this.minus(a).length();
  },

  unit: function() {
    return this.dividedBy(this.length());
  },

  // returns the vector rotated by 90 degrees clockwise
  normal: function() {
    return new CSG.Vector2D(this.y, -this.x);
  },

  // Right multiply by a 4x4 matrix (the vector is interpreted as a row vector)
  // Returns a new CSG.Vector2D
  multiply4x4: function(matrix4x4) {
    return matrix4x4.rightMultiply1x2Vector(this);
  },
};

// A polygon in 2D space:
CSG.Polygon2D = function(points, shared) {
  var vectors = [];
  if(arguments.length >= 1) {
    points.map( function(p) {
      vectors.push(new CSG.Vector2D(p) );
    });    
  }
  this.points = vectors;
  this.shared = shared;
};

CSG.Polygon2D.prototype = {
  // Matrix transformation of polygon. Returns a new CSG.Polygon2D
  transform: function(matrix4x4) {
    var newpoints = this.points.map(function(p) { return p.multiply4x4(matrix4x4); } );
    return new CSG.Polygon2D(newpoints, this.shared);
  },
  
  translate: function(v) {
    v=new CSG.Vector2D(v);
    return this.transform(CSG.Matrix4x4.translation(v.toVector3D(0)));
  },
  
  scale: function(f) {
    f=new CSG.Vector2D(f);
    return this.transform(CSG.Matrix4x4.scaling(f.toVector3D(1)));
  },
  
  rotate: function(deg) {
    return this.transform(CSG.Matrix4x4.rotationZ(deg));
  },    
  
  // convert into a CSG.Polygon; set z coordinate to the given value
  toPolygon3D: function(z) {
    var points3d=[];
    this.points.map( function(p) {
      var vec3d = p.toVector3D(z);      
      points3d.push(vec3d);
    });
    var polygon = CSG.Polygon.createFromPoints(points3d, this.shared);
    polygon.checkIfConvex();
    return polygon;
  },
  
  // extruded=shape2d.extrude({offset: [0,0,10], twistangle: 360, twiststeps: 100});
  // linear extrusion of 2D polygon, with optional twist
  // The 2d polygon is placed in in z=0 plane and extruded into direction <offset> (a CSG.Vector)
  // The final face is rotated <twistangle> degrees. Rotation is done around the origin of the 2d shape (i.e. x=0, y=0)
  // twiststeps determines the resolution of the twist (should be >= 1)  
  // returns a CSG object
  extrude: function(params) {
    // parse parameters:
    if(!params) params={};
    var offsetvector;
    if("offset" in params)
    {
      offsetvector = new CSG.Vector(params.offset); // reparse as a CSG.Vector
    }
    else
    {
      offsetvector = new CSG.Vector(0,0,1);
    }
    
    var twistangle=0;
    if("twistangle" in params)
    {
      twistangle = params.twistangle;
    }
    
    var twiststeps=10;
    if("twiststeps" in params)
    {
      twiststeps = params.twiststeps;
    }
    if(twistangle == 0) twiststeps=1;
    if(twiststeps < 1) twiststeps=1;

    // create the polygons:        
    var newpolygons = [];
    
    // bottom face polygon:
    var bottomfacepolygon = this.toPolygon3D(0);
    var direction = bottomfacepolygon.plane.normal.dot(offsetvector);
    if(direction > 0)
    {
      bottomfacepolygon = bottomfacepolygon.flipped();
    }
    newpolygons.push(bottomfacepolygon);
    
    var getTwistedPolygon = function(twiststep) {
      var fraction = (twiststep + 1) / twiststeps;
      var rotation = twistangle * fraction;
      var offset = offsetvector.times(fraction);
      var transformmatrix = CSG.Matrix4x4.rotationZ(rotation).multiply( CSG.Matrix4x4.translation(offset) );
      var polygon = bottomfacepolygon.transform(transformmatrix);      
      return polygon;
    };

    // create the side face polygons:
    var numvertices = bottomfacepolygon.vertices.length;
    var prevlevelpolygon = bottomfacepolygon;
    for(var twiststep=0; twiststep < twiststeps; ++twiststep)
    {
      var levelpolygon = getTwistedPolygon(twiststep);
      for(var i=0; i < numvertices; i++)
      {
        var sidefacepoints = [];
        var nexti = (i < (numvertices-1))? i+1:0;
        sidefacepoints.push(prevlevelpolygon.vertices[i].pos);
        sidefacepoints.push(levelpolygon.vertices[i].pos);
        sidefacepoints.push(levelpolygon.vertices[nexti].pos);
        sidefacepoints.push(prevlevelpolygon.vertices[nexti].pos);
        var sidefacepolygon=CSG.Polygon.createFromPoints(sidefacepoints, this.shared);
        newpolygons.push(sidefacepolygon);
      }
      if(twiststep == (twiststeps -1) )
      {
        // last level; add the top face polygon:
        levelpolygon = levelpolygon.flipped(); // flip so that the normal points outwards
        newpolygons.push(levelpolygon);
      }
      prevlevelpolygon = levelpolygon;
    }

    return CSG.fromPolygons(newpolygons);
  }
};

CSG.Polygon.prototype.extrude = function(offsetvector) {
  var newpolygons = [];

  var polygon1=this;
  var direction = polygon1.plane.normal.dot(offsetvector);
  if(direction > 0)
  {
    polygon1 = polygon1.flipped();
  }
  newpolygons.push(polygon1);
  var polygon2=polygon1.translate(offsetvector);
  var numvertices=this.vertices.length;
  for(var i=0; i < numvertices; i++)
  {
    var sidefacepoints = [];
    var nexti = (i < (numvertices-1))? i+1:0;
    sidefacepoints.push(polygon1.vertices[i].pos);
    sidefacepoints.push(polygon2.vertices[i].pos);
    sidefacepoints.push(polygon2.vertices[nexti].pos);
    sidefacepoints.push(polygon1.vertices[nexti].pos);
    var sidefacepolygon=CSG.Polygon.createFromPoints(sidefacepoints);
    newpolygons.push(sidefacepolygon);
  }
  polygon2 = polygon2.flipped();
  newpolygons.push(polygon2);
  return CSG.fromPolygons(newpolygons);
};

CSG.Polygon.prototype.translate = function(offset) {
  return this.transform(CSG.Matrix4x4.translation(offset));
};

// Expand the polygon with a certain radius
// This extrudes the face of the polygon and adds rounded corners 
// Returns a CSG object (not a polygon anymore!)
// resolution: number of points per 360 degree for the rounded corners
CSG.Polygon.prototype.expand = function(radius, resolution) {
  resolution=resolution || 8;
  var result=new CSG();
  
  // expand each side of the polygon. The expansion of a line is a cylinder with
  // two spheres at the end:
  var numvertices=this.vertices.length;
  for(var i=0; i < numvertices; i++)
  {
    var previ = (i == 0) ? (numvertices-1):i-1;
    var p1 = this.vertices[previ].pos;
    var p2 = this.vertices[i].pos;
    var cylinder=CSG.cylinder({start: p1, end: p2, radius: radius, slices: resolution});
    result = result.union(cylinder);
    var sphere = CSG.sphere({center: p1, radius: radius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
  }
  var extrudevector=this.plane.normal.unit().times(2*radius);
  var translatedpolygon = this.translate(extrudevector.times(-0.5));
  var extrudedface = translatedpolygon.extrude(extrudevector);  
  result=result.union(extrudedface);
  return result;
};

// Expand the solid
// resolution: number of points per 360 degree for the rounded corners
CSG.prototype.expand = function(radius, resolution) {
  var result=this;
  this.polygons.map(function(p) {
    var expanded=p.expand(radius, resolution);
    result=result.union(expanded);
  });
  return result;
};

// Contract the solid
// resolution: number of points per 360 degree for the rounded corners
CSG.prototype.contract = function(radius, resolution) {
  var result=this;
  this.polygons.map(function(p) {
    var expanded=p.expand(radius, resolution);
    result=result.subtract(expanded);
  });
  return result;
};

CSG.roundedCube = function(cuberadius, roundradius, resolution) {
  resolution = resolution || 8;
  cuberadius=new CSG.Vector(cuberadius);
  var innercuberadius=cuberadius.clone();
  innercuberadius.x -= roundradius;
  innercuberadius.y -= roundradius;
  innercuberadius.z -= roundradius;
  var result = CSG.cube({radius: [cuberadius.x, innercuberadius.y, innercuberadius.z]});
  result = result.union( CSG.cube({radius: [innercuberadius.x, cuberadius.y, innercuberadius.z]}));
  result = result.union( CSG.cube({radius: [innercuberadius.x, innercuberadius.y, cuberadius.z]}));
  for(var level=0; level < 2; level++)
  {
    var z = innercuberadius.z;
    if(level == 1) z = -z;
    var p1 = new CSG.Vector(innercuberadius.x, innercuberadius.y, z);
    var p2 = new CSG.Vector(innercuberadius.x, -innercuberadius.y, z);
    var p3 = new CSG.Vector(-innercuberadius.x, -innercuberadius.y, z);
    var p4 = new CSG.Vector(-innercuberadius.x, innercuberadius.y, z);
    var sphere = CSG.sphere({center: p1, radius: roundradius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
    sphere = CSG.sphere({center: p2, radius: roundradius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
    sphere = CSG.sphere({center: p3, radius: roundradius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
    sphere = CSG.sphere({center: p4, radius: roundradius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
    var cylinder = CSG.cylinder({start:p1, end: p2, radius: roundradius, slices: resolution});
    result = result.union(cylinder);
    cylinder = CSG.cylinder({start:p2, end: p3, radius: roundradius, slices: resolution});
    result = result.union(cylinder);
    cylinder = CSG.cylinder({start:p3, end: p4, radius: roundradius, slices: resolution});
    result = result.union(cylinder);
    cylinder = CSG.cylinder({start:p4, end: p1, radius: roundradius, slices: resolution});
    result = result.union(cylinder);
    if(level == 0) {
      var d = new CSG.Vector(0, 0, -2*z);
      cylinder = CSG.cylinder({start:p1, end: p1.plus(d), radius: roundradius, slices: resolution});
      result = result.union(cylinder);
      cylinder = CSG.cylinder({start:p2, end: p2.plus(d), radius: roundradius, slices: resolution});
      result = result.union(cylinder);
      cylinder = CSG.cylinder({start:p3, end: p3.plus(d), radius: roundradius, slices: resolution});
      result = result.union(cylinder);
      cylinder = CSG.cylinder({start:p4, end: p4.plus(d), radius: roundradius, slices: resolution});
      result = result.union(cylinder);
    }
  }
  return result;
}


// # class Line2D

// Represents a directional line in 2D space
// A line is parametrized by its normal vector (perpendicular to the line, rotated 90 degrees counter clockwise)
// and w. The line passes through the point <normal>.times(w).
// normal must be a unit vector!
// Equation: p is on line if normal.dot(p)==w
CSG.Line2D = function(normal, w) {
  this.normal = normal;
  this.w = w;
};

CSG.Line2D.fromPoints = function(p1, p2) {
  var direction = p2.minus(p1);
  var normal = direction.normal().negated().unit();
  var w = p1.dot(normal);
  return new CSG.Line2D(normal, w); 
};

CSG.Line2D.prototype = {
  // same line but opposite direction:
  inverse: function() {
    return new CSG.Line2D(this.normal.negated(), -this.w);
  },
  
  equals: function(l) {
    return (l.normal.equals(this.normal) && (l.w == this.w));
  },
  
  origin: function() {
    return this.normal.times(this.w);
  },

  direction: function() {
    return this.normal.normal(); 
  },
  
  xAtY: function(y) {
    // (py == y) && (normal * p == w)
    // -> px = (w - normal.y * y) / normal.x
    var x = (this.w - this.normal.y * y) / this.normal.x;
    return x; 
  },
  
  absDistanceToPoint: function(point) {
    var point_projected = point.dot(this.normal);
    var distance = Math.abs(point_projected - this.w);
    return distance;
  },
  
  closestPoint: function(point) {
    var vector = point.dot(this.direction());
    return origin.plus(vector);  
  },
};

// # class Line3D

// Represents a line in 3D space
// direction must be a unit vector 
// point is a random point on the line
CSG.Line3D = function(point, direction) {
  this.point = point;
  this.direction = direction;
};

CSG.Line3D.fromPoints = function(p1, p2) {
  // make the line (somewhat) canonical by always using the point closest to the origin
  // as the reference point
  var direction = p2.minus(p1).unit();
  var l1 = p1.dot(p1);
  var l2 = p2.dot(p2);
  if(l1 > l2)
  {
    return new CSG.Line3D(p2, direction);
  }
  else
  {
    return new CSG.Line3D(p1, direction);
  }
};

CSG.Line3D.fromPlanes = function(p1, p2) {
  var direction = p1.normal.cross(p2.normal);
  var l=direction.length();
  if(l < 1e-10)
  {
    throw new Error("Parallel planes");
  }
  direction = direction.times(1.0/l);

  var mabsx = Math.abs(direction.x);
  var mabsy = Math.abs(direction.y);
  var mabsz = Math.abs(direction.z);
  var origin;
  if( (mabsx >= mabsy) && (mabsx >= mabsz) )
  {
    // direction vector is mostly pointing towards x
    // find a point p for which x is zero:
    var r = CSG.Line3D.Solve2Linear(p1.normal.y, p1.normal.z, p2.normal.y, p2.normal.z, p1.w, p2.w);    
    origin = new CSG.Vector(0, r[0], r[1]);
  }
  else if( (mabsy >= mabsx) && (mabsy >= mabsz) )
  {
    // find a point p for which y is zero:
    var r = CSG.Line3D.Solve2Linear(p1.normal.x, p1.normal.z, p2.normal.x, p2.normal.z, p1.w, p2.w);    
    origin = new CSG.Vector(r[0], 0, r[1]);
  }
  else
  {
    // find a point p for which z is zero:
    var r = CSG.Line3D.Solve2Linear(p1.normal.x, p1.normal.y, p2.normal.x, p2.normal.y, p1.w, p2.w);    
    origin = new CSG.Vector(r[0], r[1], 0);
  }
  return new CSG.Line3D(origin, direction);
};

// solve
// [ab][x] = [u]
// [cd][y]   [v]
CSG.Line3D.Solve2Linear = function(a,b,c,d,u,v) {
  var det = a*d - b*c;
  var invdet = 1.0/det;
  var x = u*d - b*v;
  var y = -u*c + a*v;
  x *= invdet;
  y *= invdet;
  return [x,y];
};

CSG.Line3D.prototype = {
  intersectWithPlane: function(plane) {
    // plane: plane.normal * p = plane.w
    // line: p=line.point + labda * line.direction
    var labda = (plane.w - plane.normal.dot(this.point)) / plane.normal.dot(this.direction);
    var point = this.point.plus(this.direction.times(labda));
    return point;
  },
  
  clone: function(line) {
    return new CSG.Line3D(this.point.clone(), this.direction.clone());
  },
  
  reverse: function() {
    return new CSG.Line3D(this.point.clone(), this.direction.negated());
  },
  
  transform: function(matrix4x4) {
    var newpoint = this.point.multiply4x4(matrix4x4);
    var pointPlusDirection = this.point.plus(this.direction);
    var newPointPlusDirection = pointPlusDirection.multiply4x4(matrix4x4);
    var newdirection = newPointPlusDirection.minus(newpoint);
    return new CSG.Line3D(newpoint, newdirection);  
  },  
  
  closestPointOnLine: function(point) {
    var t = point.minus(this.point).dot(this.direction) / this.direction.dot(this.direction);
    var closestpoint = this.point.plus(this.direction.times(t));
    return closestpoint;
  },
  
  distanceToPoint: function(point) {
    var closestpoint = this.closestPointOnLine(point);
    var distancevector = point.minus(closestpoint);
    var distance = distancevector.length();
    return distance;
  },
  
  equals: function(line3d) {
    if(!this.direction.equals(line3d.direction)) return false;
    var distance = this.distanceToPoint(line3d.point);
    if(distance > 1e-8) return false;
    return true;    
  },
};

// # class OrthoNormalBasis

CSG.OrthoNormalBasis = function (plane, rot90) {
  var rightvector;
  if(Math.abs(plane.normal.x) > Math.abs(plane.normal.y))
  {
    rightvector = new CSG.Vector(0, 1, 0);
  }
  else
  {
    rightvector = new CSG.Vector(1, 0, 0);
  }
  this.v = rightvector.cross(plane.normal).unit();
  this.u = plane.normal.cross(this.v);
  if(rot90)
  {
    var v = this.u.negated();
    this.u = this.v;
    this.v = v;    
  }
  this.planeorigin = plane.normal.times(plane.w);
};

CSG.OrthoNormalBasis.prototype = {
  to2D: function(vec3) {
    return new CSG.Vector2D(vec3.dot(this.u), vec3.dot(this.v));
  },
  
  to3D: function(vec2) {
    return this.planeorigin.plus(this.u.times(vec2.x)).plus(this.v.times(vec2.y));
  },
  
  line3Dto2D: function(line3d) {
    var a = line3d.point;
    var b = line3d.direction.plus(a);
    var a2d = this.to2D(a);
    var b2d = this.to2D(b);
    return CSG.Line2D.fromPoints(a2d, b2d);
  },

  line2Dto3D: function(line2d) {
    var a = line2d.origin();
    var b = line2d.direction().plus(a);
    var a3d = this.to3D(a);
    var b3d = this.to3D(b);
    return CSG.Line3D.fromPoints(a3d, b3d);
  },
};

function insertSorted(array, element, comparefunc) {
  var leftbound = 0;
  var rightbound = array.length;
  while(rightbound > leftbound)
  {
    var testindex = Math.floor( (leftbound + rightbound) / 2);
    var testelement = array[testindex];
    var compareresult = comparefunc(element, testelement);
    if(compareresult > 0)   // element > testelement
    {
      leftbound = testindex + 1;
    }
    else
    {
      rightbound = testindex;
    }    
  }
  array.splice(leftbound,0,element);
}

// Get the x coordinate of a point with a certain y coordinate, interpolated between two
// points (CSG.Vector2D).
// Interpolation is robust even if the points have the same y coordinate
CSG.interpolateBetween2DPointsForY = function(point1, point2, y) {
  var f1 = y - point1.y;
  var f2 = point2.y - point1.y;
  if(f2 < 0)
  {
    f1 = -f1;
    f2 = -f2;
  }
  var t;
  if(f1 <= 0)
  {
    t = 0.0;
  }
  else if(f1 >= f2)
  {
    t = 1.0;
  }
  else if(f2 < 1e-10)
  {
    t = 0.5;
  }
  else
  {
    t = f1 / f2;
  }
  var result = point1.x + t * (point2.x - point1.x);
  return result;
};


CSG.tesselate = function(sourcepolygons, destpolygons, rot90)
{
  var EPS = 1e-5;
  var numpolygons = sourcepolygons.length;
  if(numpolygons > 0)
  {
    var plane = sourcepolygons[0].plane;
    var orthobasis = new CSG.OrthoNormalBasis(plane, rot90);
    var polygonvertices2d = [];    // array of array of CSG.Vector2D
    var polygontopvertexindexes = []; // array of indexes of topmost vertex per polygon
    var topy2polygonindexes = {};
    var ycoordinatetopolygonindexes = {};
    
    var xcoordinatebins = {};    
    var ycoordinatebins = {};    
    
    // convert all polygon vertices to 2D
    // Make a list of all encountered y coordinates
    // And build a map of all polygons that have a vertex at a certain y coordinate:    
    for(var polygonindex=0; polygonindex < numpolygons; polygonindex++)
    {
      var poly3d = sourcepolygons[polygonindex];
      var vertices2d = [];
      var numvertices = poly3d.vertices.length;
      var minindex = -1;
      if(numvertices > 0)
      {
        var miny, maxy, maxindex;
        for(var i=0; i < numvertices; i++)
        {
          var pos2d = orthobasis.to2D(poly3d.vertices[i].pos);
          // perform binning of coordinates: If we have multiple vertexes very
          // close to each other, give them the same coordinates:
          var xcoordinatebin = Math.floor(pos2d.x * 1e6);
          var ycoordinatebin = Math.floor(pos2d.y * 1e6);
/*          
          if(xcoordinatebin in xcoordinatebins)
          {
            pos2d.x = xcoordinatebins[xcoordinatebin];
          }
          else if(xcoordinatebin+1 in xcoordinatebins)
          {
            pos2d.x = xcoordinatebins[xcoordinatebin+1];
          }
          else if(xcoordinatebin-1 in xcoordinatebins)
          {
            pos2d.x = xcoordinatebins[xcoordinatebin-1];
          }
          else
          {
            xcoordinatebins[xcoordinatebin] = pos2d.x;
          }
*/          
          if(ycoordinatebin in ycoordinatebins)
          {
            pos2d.y = ycoordinatebins[ycoordinatebin];
          }
          else if(ycoordinatebin+1 in ycoordinatebins)
          {
            pos2d.y = ycoordinatebins[ycoordinatebin+1];
          }
          else if(ycoordinatebin-1 in ycoordinatebins)
          {
            pos2d.y = ycoordinatebins[ycoordinatebin-1];
          }
          else
          {
            ycoordinatebins[ycoordinatebin] = pos2d.y;
          }
          vertices2d.push(pos2d);
          var y = pos2d.y;
          if( (i == 0) || (y < miny) )
          {
            miny = y;
            minindex = i;
          }
          if( (i == 0) || (y > maxy) )
          {
            maxy = y;
            maxindex = i;
          }
          if(! (y in ycoordinatetopolygonindexes))
          {
            ycoordinatetopolygonindexes[y] = {};
          }
          ycoordinatetopolygonindexes[y][polygonindex]=true;
        }
        if(miny >= maxy)
        {
          // degenerate polygon, all vertices have same y coordinate. Just ignore it from now:
          vertices2d = [];
        }
        else
        {
          if(! (miny in topy2polygonindexes))
          {
            topy2polygonindexes[miny] = [];
          }
          topy2polygonindexes[miny].push(polygonindex);          
        }
      }  // if(numvertices > 0)
      polygonvertices2d.push(vertices2d); 
      polygontopvertexindexes.push(minindex); 
    }
    var ycoordinates = [];
    for(var ycoordinate in ycoordinatetopolygonindexes) ycoordinates.push(ycoordinate);
    ycoordinates.sort(function(a,b) {return a-b});

    // Now we will iterate over all y coordinates, from lowest to highest y coordinate
    // activepolygons: source polygons that are 'active', i.e. intersect with our y coordinate
    var activepolygons = [];
    var prevoutpolygonrow = [];
    for(var yindex = 0; yindex < ycoordinates.length; yindex++)
    {
      var newoutpolygonrow = [];
      var ycoordinate_as_string = ycoordinates[yindex];
      var ycoordinate = Number(ycoordinate_as_string);
      
      // update activepolygons for this y coordinate:
      // - Remove any polygons that end at this y coordinate
      // - update leftvertexindex and rightvertexindex (which point to the current vertex index 
      //   at the the left and right side of the polygon
      // Iterate over all polygons that have a corner at this y coordinate:
      var polygonindexeswithcorner = ycoordinatetopolygonindexes[ycoordinate_as_string];
      for(var activepolygonindex = 0; activepolygonindex < activepolygons.length; ++activepolygonindex)  
      {
        var activepolygon = activepolygons[activepolygonindex];
        var polygonindex = activepolygon.polygonindex;
        if(polygonindexeswithcorner[polygonindex])
        {
          // this active polygon has a corner at this y coordinate:
          var vertices2d = polygonvertices2d[polygonindex];
          var numvertices = vertices2d.length;
          var newleftvertexindex = activepolygon.leftvertexindex;
          var newrightvertexindex = activepolygon.rightvertexindex;
          // See if we need to increase leftvertexindex or decrease rightvertexindex:
          while(true)
          {
            var nextleftvertexindex = newleftvertexindex+1;
            if(nextleftvertexindex >= numvertices) nextleftvertexindex = 0;
            if(vertices2d[nextleftvertexindex].y != ycoordinate) break;
            newleftvertexindex = nextleftvertexindex;
          }
          var nextrightvertexindex = newrightvertexindex-1;
          if(nextrightvertexindex < 0) nextrightvertexindex = numvertices-1;
          if(vertices2d[nextrightvertexindex].y == ycoordinate)
          {
            newrightvertexindex = nextrightvertexindex;
          }
          if( (newleftvertexindex != activepolygon.leftvertexindex) && (newleftvertexindex == newrightvertexindex) )
          {
            // We have increased leftvertexindex or decreased rightvertexindex, and now they point to the same vertex
            // This means that this is the bottom point of the polygon. We'll remove it:
            activepolygons.splice(activepolygonindex, 1);
            --activepolygonindex;            
          }
          else 
          {
            activepolygon.leftvertexindex = newleftvertexindex;
            activepolygon.rightvertexindex = newrightvertexindex;
            activepolygon.topleft = vertices2d[newleftvertexindex];
            activepolygon.topright = vertices2d[newrightvertexindex];
            var nextleftvertexindex = newleftvertexindex+1;
            if(nextleftvertexindex >= numvertices) nextleftvertexindex = 0;
            activepolygon.bottomleft = vertices2d[nextleftvertexindex];
            var nextrightvertexindex = newrightvertexindex-1;
            if(nextrightvertexindex < 0) nextrightvertexindex = numvertices-1;
            activepolygon.bottomright = vertices2d[nextrightvertexindex];            
          } 
        } // if polygon has corner here
      }  // for activepolygonindex

      var nextycoordinate;      
      if(yindex >= ycoordinates.length-1)
      {
        // last row, all polygons must be finished here:
        activepolygons = [];
        nextycoordinate = null;
      }
      else // yindex < ycoordinates.length-1
      {
        nextycoordinate = Number(ycoordinates[yindex+1]);
        var middleycoordinate = 0.5 * (ycoordinate + nextycoordinate);
        // update activepolygons by adding any polygons that start here: 
        var startingpolygonindexes = topy2polygonindexes[ycoordinate_as_string];      
        for(var polygonindex_key in startingpolygonindexes)
        {
          var polygonindex = startingpolygonindexes[polygonindex_key];
          var vertices2d = polygonvertices2d[polygonindex];
          var numvertices = vertices2d.length;
          var topvertexindex = polygontopvertexindexes[polygonindex];
          // the top of the polygon may be a horizontal line. In that case topvertexindex can point to any point on this line.
          // Find the left and right topmost vertices which have the current y coordinate:
          var topleftvertexindex = topvertexindex;
          while(true)
          {
            var i = topleftvertexindex + 1;
            if(i >= numvertices) i = 0;
            if(vertices2d[i].y != ycoordinate) break;
            if(i == topvertexindex) break; // should not happen, but just to prevent endless loops
            topleftvertexindex = i;          
          }
          var toprightvertexindex = topvertexindex;
          while(true)
          {
            var i = toprightvertexindex - 1;
            if(i < 0) i = numvertices - 1;
            if(vertices2d[i].y != ycoordinate) break;
            if(i == topleftvertexindex) break; // should not happen, but just to prevent endless loops
            toprightvertexindex = i;          
          }
          var nextleftvertexindex = topleftvertexindex+1;
          if(nextleftvertexindex >= numvertices) nextleftvertexindex = 0;
          var nextrightvertexindex = toprightvertexindex-1;
          if(nextrightvertexindex < 0) nextrightvertexindex = numvertices-1;
          var newactivepolygon = {
            polygonindex: polygonindex,
            leftvertexindex: topleftvertexindex,
            rightvertexindex: toprightvertexindex,
            topleft: vertices2d[topleftvertexindex],
            topright: vertices2d[toprightvertexindex],
            bottomleft: vertices2d[nextleftvertexindex],
            bottomright: vertices2d[nextrightvertexindex],            
          };
          insertSorted(activepolygons, newactivepolygon, function(el1, el2) {
            var x1 = CSG.interpolateBetween2DPointsForY(el1.topleft, el1.bottomleft, middleycoordinate); 
            var x2 = CSG.interpolateBetween2DPointsForY(el2.topleft, el2.bottomleft, middleycoordinate); 
            if(x1 > x2) return 1;
            if(x1 < x2) return -1;
            return 0;
          });
        } // for(var polygonindex in startingpolygonindexes)
      } //  yindex < ycoordinates.length-1
      //if( (yindex == ycoordinates.length-1) || (nextycoordinate - ycoordinate > EPS) )
      if(true)
      {
        // Now activepolygons is up to date
        // Build the output polygons for the next row in newoutpolygonrow:
        for(var activepolygon_key in activepolygons)
        {
          var activepolygon = activepolygons[activepolygon_key];
          var polygonindex = activepolygon.polygonindex;
          var vertices2d = polygonvertices2d[polygonindex];
          var numvertices = vertices2d.length;

          var x = CSG.interpolateBetween2DPointsForY(activepolygon.topleft, activepolygon.bottomleft, ycoordinate);          
          var topleft=new CSG.Vector2D(x, ycoordinate); 
          x = CSG.interpolateBetween2DPointsForY(activepolygon.topright, activepolygon.bottomright, ycoordinate);          
          var topright=new CSG.Vector2D(x, ycoordinate); 
          x = CSG.interpolateBetween2DPointsForY(activepolygon.topleft, activepolygon.bottomleft, nextycoordinate);          
          var bottomleft=new CSG.Vector2D(x, nextycoordinate); 
          x = CSG.interpolateBetween2DPointsForY(activepolygon.topright, activepolygon.bottomright, nextycoordinate);          
          var bottomright=new CSG.Vector2D(x, nextycoordinate);                      
          var outpolygon = {
            topleft: topleft, 
            topright: topright,
            bottomleft: bottomleft, 
            bottomright: bottomright,
            leftline: CSG.Line2D.fromPoints(topleft, bottomleft),
            rightline: CSG.Line2D.fromPoints(bottomright, topright),
          };
          if(newoutpolygonrow.length > 0)
          {
            var prevoutpolygon = newoutpolygonrow[newoutpolygonrow.length - 1];
            var d1 = outpolygon.topleft.distanceTo(prevoutpolygon.topright);
            var d2 = outpolygon.bottomleft.distanceTo(prevoutpolygon.bottomright);
            if( (d1 < EPS) && (d2 < EPS) )
            {          
              // we can join this polygon with the one to the left:
              outpolygon.topleft = prevoutpolygon.topleft;
              outpolygon.leftline = prevoutpolygon.leftline;            
              outpolygon.bottomleft = prevoutpolygon.bottomleft;
              newoutpolygonrow.splice(newoutpolygonrow.length - 1, 1);
            }          
          }
          newoutpolygonrow.push(outpolygon);
        } // for(activepolygon in activepolygons)
        if(yindex > 0)
        {
          // try to match the new polygons against the previous row:
          var prevcontinuedindexes = {};
          var matchedindexes = {};
          for(var i = 0; i < newoutpolygonrow.length; i++)
          {
            var thispolygon = newoutpolygonrow[i];
            for(var ii = 0; ii < prevoutpolygonrow.length; ii++)
            {
              if(!matchedindexes[ii])   // not already processed?
              {
                // We have a match if the sidelines are equal or if the top coordinates
                // are on the sidelines of the previous polygon
                var prevpolygon = prevoutpolygonrow[ii];
                if(prevpolygon.bottomleft.distanceTo(thispolygon.topleft) < EPS)
                {
                  if(prevpolygon.bottomright.distanceTo(thispolygon.topright) < EPS)
                  {
                    // Yes, the top of this polygon matches the bottom of the previous:
                    matchedindexes[ii] = true;
                    // Now check if the joined polygon would remain convex:
                    var d1 = thispolygon.leftline.direction().x - prevpolygon.leftline.direction().x;
                    var d2 = thispolygon.rightline.direction().x - prevpolygon.rightline.direction().x;                    
                    var leftlinecontinues = Math.abs(d1) < EPS;
                    var rightlinecontinues = Math.abs(d2) < EPS;
                    var leftlineisconvex = leftlinecontinues || (d1 >= 0);
                    var rightlineisconvex = rightlinecontinues || (d2 >= 0);
                    if(leftlineisconvex && rightlineisconvex)
                    {
                      // yes, both sides have convex corners:
                      // This polygon will continue the previous polygon
                      thispolygon.outpolygon = prevpolygon.outpolygon;
                      thispolygon.leftlinecontinues = leftlinecontinues;
                      thispolygon.rightlinecontinues = rightlinecontinues;
                      prevcontinuedindexes[ii] = true;
                    }
                    break;                  
                  }
                }
              } // if(!prevcontinuedindexes[ii])
            } // for ii
          } // for i
          for(var ii = 0; ii < prevoutpolygonrow.length; ii++)
          {
            if(!prevcontinuedindexes[ii])
            {
              // polygon ends here
              // Finish the polygon with the last point(s):
              var prevpolygon = prevoutpolygonrow[ii];
              prevpolygon.outpolygon.rightpoints.push(prevpolygon.bottomright);
              if(prevpolygon.bottomright.distanceTo(prevpolygon.bottomleft) > EPS)
              {
                // polygon ends with a horizontal line:
                prevpolygon.outpolygon.leftpoints.push(prevpolygon.bottomleft);
              }
              // reverse the right half so we get a counterclockwise circle:
              prevpolygon.outpolygon.rightpoints.reverse();
              var points2d = prevpolygon.outpolygon.leftpoints.concat(prevpolygon.outpolygon.rightpoints); 
              var vertices3d = [];
              points2d.map(function(point2d) {
                var point3d = orthobasis.to3D(point2d);
                var vertex3d = new CSG.Vertex(point3d);
                vertices3d.push(vertex3d);              
              });
              var shared = null;
              var polygon = new CSG.Polygon(vertices3d, shared, plane);
              destpolygons.push(polygon);
            }
          }                
        } // if(yindex > 0)
        for(var i = 0; i < newoutpolygonrow.length; i++)
        {
          var thispolygon = newoutpolygonrow[i];
          if(!thispolygon.outpolygon)
          {
            // polygon starts here:
            thispolygon.outpolygon = {
              leftpoints: [],
              rightpoints: [],
            };
            thispolygon.outpolygon.leftpoints.push(thispolygon.topleft);
            if(thispolygon.topleft.distanceTo(thispolygon.topright) > EPS)
            {
              // we have a horizontal line at the top:
              thispolygon.outpolygon.rightpoints.push(thispolygon.topright);
            }
          }
          else
          {
            // continuation of a previous row
            if(! thispolygon.leftlinecontinues )
            {
              thispolygon.outpolygon.leftpoints.push(thispolygon.topleft);
            }
            if(! thispolygon.rightlinecontinues )
            {
              thispolygon.outpolygon.rightpoints.push(thispolygon.topright);
            }
          }
        }
        prevoutpolygonrow = newoutpolygonrow;
      }
    } // for yindex 
  } // if(numpolygons > 0)
}

/*
CSG.testTesselate = function()
{
  var plane = new CSG.Plane(new CSG.Vector(0,0,1), 0);  // the z=0 plane
  var shared = null;

  var inpolygons = [];
  inpolygons.push(CSG.Polygon.createFromPoints([[1, 1, 0], [3, 1, 0], [2, 1.5, 0]], shared, plane));
  inpolygons.push(CSG.Polygon.createFromPoints([[1, 1, 0], [4, 0, 0], [3, 1, 0]], shared, plane));
  inpolygons.push(CSG.Polygon.createFromPoints([[0, 0, 0], [4, 0, 0], [1, 1, 0]], shared, plane));

  var outpolygons = [];
  CSG.tesselate(inpolygons, outpolygons);
};
*/

CSG.prototype.canonicalized = function () {
  if(this.isCanonicalized)
  {
    return this;
  }
  else
  {
    var factory = new CSG.fuzzyCSGFactory();
    var result = factory.getCSG(this);
    result.isCanonicalized = true;
    return result;
  }
};

CSG.prototype.reTesselate = function () {
  var csg = this.canonicalized();
  var polygonsPerPlane = {};
  csg.polygons.map(function(polygon) {
    var planetag = polygon.plane.getTag();
    if(! (planetag in polygonsPerPlane) )
    {
      polygonsPerPlane[planetag] = [];
    }
    polygonsPerPlane[planetag].push(polygon);
  });
  var destpolygons = [];
  for(planetag in polygonsPerPlane)
  {
    var sourcepolygons = polygonsPerPlane[planetag];
    if(sourcepolygons.length < 2)
    {
      destpolygons = destpolygons.concat(sourcepolygons);
    }
    else
    {
      var retesselayedpolygons = [];
      CSG.tesselate(sourcepolygons, retesselayedpolygons, false);
      destpolygons = destpolygons.concat(retesselayedpolygons);
    }
  }
  var result = CSG.fromPolygons(destpolygons);
  result = result.canonicalized(); 
  return result;
};

////////////////////////////////

CSG.fuzzyFactory = function(numdimensions, tolerance) {
  var lookuptable = [];
  for(var i=0; i < numdimensions; i++)
  {
    lookuptable.push({});
  }
  this.lookuptable = lookuptable;
  this.nextElementId = 1;
  this.multiplier = 1.0 / tolerance;
  this.objectTable = {};
};

CSG.fuzzyFactory.insertKey = function(key, lookuptable, quantizedvalue) {
  if(quantizedvalue in lookuptable)
  {
    lookuptable[quantizedvalue][key] = true;
  }
  else
  {
    var newset = {};
    newset[key] = true;
    lookuptable[quantizedvalue] = newset;
  }
};

CSG.fuzzyFactory.prototype = {
  // var obj = f.lookupOrCreate([el1, el2, el3], function(elements) {/* create the new object */});
  // Performs a fuzzy lookup of the object with the specified elements.
  // If found, returns the existing object
  // If not found, calls the supplied callback function which should create a new object with
  // the specified properties. This object is inserted in the lookup database.
  lookupOrCreate: function(els, creatorCallback) {
    var object;
    var key = this.lookupKey(els);
    if(key === null)
    {
      object = creatorCallback(els);
      key = this.nextElementId++;
      this.objectTable[key] = object;
      for(var dimension = 0; dimension < els.length; dimension++)
      {
        var elementLookupTable = this.lookuptable[dimension];
        var value = els[dimension];
        var valueMultiplied = value * this.multiplier;
        var valueQuantized1 = Math.floor(valueMultiplied);
        var valueQuantized2 = Math.ceil(valueMultiplied);
        CSG.fuzzyFactory.insertKey(key, elementLookupTable, valueQuantized1);
        CSG.fuzzyFactory.insertKey(key, elementLookupTable, valueQuantized2);
      }      
    }
    else
    {
      object = this.objectTable[key];
    }
    return object;
  },

  // ----------- PRIVATE METHODS:
  lookupKey: function(els) {
    var keyset = {};
    for(var i=0; i < els.length; i++)
    {
      var r = this.lookupKeySetForDimension(i, els[i]);
      if(i == 0)
      {
        keyset = r;
      }
      else
      {
        keyset = CSG.fuzzyFactory.intersectSets(keyset, r);
      }
      if(CSG.fuzzyFactory.isEmptySet(keyset)) return null;
    }
    // return first matching key:
    for(var key in keyset) return key;
    return null;
  },

  lookupKeySetForDimension: function(dimension, value) {
    var result;
    var elementLookupTable = this.lookuptable[dimension];
    var valueMultiplied = value * this.multiplier;
    var valueQuantized1 = Math.floor(valueMultiplied);
    var valueQuantized2 = Math.ceil(valueMultiplied);
    if(valueQuantized1 in elementLookupTable)
    {
      result = elementLookupTable[valueQuantized1]; 
      if(valueQuantized2 in elementLookupTable)
      {
        result = CSG.fuzzyFactory.joinSets(result, elementLookupTable[valueQuantized2]);
      }
    }
    else if(valueQuantized2 in elementLookupTable)
    {
      result = elementLookupTable[valueQuantized2]; 
    }
    else
    {
      result = {};
    }
    return result;
  },
};

CSG.fuzzyFactory.isEmptySet = function(obj) {
  for(var key in obj) return false;
  return true;
};

CSG.fuzzyFactory.intersectSets = function(set1, set2) {
  var result = {};
  for(var key in set1)
  {
    if(key in set2)
    {
      result[key] = true;
    }
  }
  return result;
};

CSG.fuzzyFactory.joinSets = function(set1, set2) {
  var result = {};
  for(var key in set1)
  {
    result[key] = true;
  }
  for(var key in set2)
  {
    result[key] = true;
  }
  return result;
};

//////////////////////////////////////

CSG.fuzzyCSGFactory = function() {
  this.vertexfactory = new CSG.fuzzyFactory(3, 1e-5);
  this.planefactory = new CSG.fuzzyFactory(4, 1e-5);
};

CSG.fuzzyCSGFactory.prototype = {
  getVertex: function(sourcevertex) {
    var elements = [sourcevertex.pos.x, sourcevertex.pos.y, sourcevertex.pos.z]; 
    var result = this.vertexfactory.lookupOrCreate(elements, function(els) {
      return sourcevertex;
    });
    return result;
  },
  
  getPlane: function(sourceplane) {
    var elements = [sourceplane.normal.x, sourceplane.normal.y, sourceplane.normal.z, sourceplane.w]; 
    var result = this.planefactory.lookupOrCreate(elements, function(els) {
      return sourceplane;
    });
    return result;
  },

  getPolygon: function(sourcepolygon) {
    var newplane = this.getPlane(sourcepolygon.plane);
    var _this = this;
    var newvertices = sourcepolygon.vertices.map(function(vertex) {
      return _this.getVertex(vertex);
    });
    return new CSG.Polygon(newvertices, sourcepolygon.shared, newplane);    
  },

  getCSG: function(sourcecsg) {
    var _this = this;
    var newpolygons = sourcecsg.polygons.map(function(polygon) {
      return _this.getPolygon(polygon);
    });
    return CSG.fromPolygons(newpolygons);
  },
};

//////////////////////////////////////

CSG.staticTag = 1;

CSG.getTag = function () {
  return CSG.staticTag++;
};