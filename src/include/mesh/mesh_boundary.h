/*
#define MESH_BC_INTERFACE           0
#define MESH_BC_INLET               1
#define MESH_BC_OUTLET              2
#define MESH_BC_ISOTHERMAL_WALL     3
#define MESH_BC_ABDIABAT_WALL       4
#define MESH_BC_SYMMETRY            5
#define MESH_BC_PERIODIC            6
*/

enum {
    MeshBC_interface,
    MeshBC_inlet,
    MeshBC_outlet,
    MeshBC_isothermal_wall,
    MeshBC_abdiabat_wall,
    MeshBC_symmetry,
    MeshBC_periodic,
};
