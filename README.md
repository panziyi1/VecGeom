Overview
========

What is VecGeom ?
-----------------
VecGeom is a library of that provides geometrical shape primitives, 
ways to describe a model geometry and to navigate within it, 
tailored for use in cutting edge particle transport simulation

Shape primitives have the full set of methods for intersection, 
distance to nearest boundary of volume required for navigation of tracks 
for particle transport simulations.  
A distinguishing feature of VecGeom is that it provides methods with SIMD
signatures to cope with multiple tracks in one call, for all appropriate 
methods of the shape/solid primitives.

List of topics from Google document
- How we create geometry objects via factories, 
- The roles of ‘unplaced’, ‘placed’ and ‘specialised’ classes
- The struct describing each shape
- Navigation techniques/features

Geometry volume primitives 
---------------------------
### Unplaced volumes (solids) 
An unplaced volume represents a geometry shape (primitive) and offers
interfaces to query distance, location, containment, etc. 

The volume is typically placed in its "natural" system of coordinates, e.g. a  the center of coordinates for a sphere and for a 
rectangular parallelilepiped ('box') are at their centers.

This is the same concept as Geant4 'solids' (G4VSolid) and TGeo 'shape' (TGeoShape).

### Placed volumes (solids) and specialisations
A placed volume represents a geometry shape which has been located 
with either at a different position or with a rotation or both.
It is the primary object from which a geometry model of a setup is created. 

For reasons of efficiency different versions exist for each combination, 
depending on whether or not a translation or a rotation exists.

### How we create geometry objects via factories
These different specialisations of placed volume can be created using a factor 

[ Creating a shape/solid optimally in VecGeom ]

### The struct describing each shape

How to create a geometry setup in VecGeom
-----------------------------------------
### Logical volumes

Navigation techniques/features
------------------------------




