/**************************************
 * meso  - MESH module                *
 *  contains:                         *
 *         mesh object, geom function *
 *                Nov 6, 2023  by MYC *
 **************************************/

#ifndef HEADER_MESH
#define HEADER_MESH

#ifndef HEADER_CORE
#include "core.h"
#endif

/// MESH boundary
#include "mesh/mesh_boundary.h"
/// MESH object
#include "mesh/mesh_object.h"
/// MESH geom
#include "mesh/mesh_geom.h"
/// neu reader
#include "mesh/neu_reader.h"
/// MESH generator
#include "mesh/mesh_generator.h"
/// MESH writer
#include "mesh/mesh_writer.h"

/// Discrete Phase - Particle Tracker
#include "mesh/dpm_object.h"

#endif  // HEADER_MESH
