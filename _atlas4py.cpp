/*cppimport
<%
cfg['compiler_args'] = ['-std=c++17', '-fopenmp', '-O0']
cfg['linker_args'] = [
    '-L/home/lukas/documents/work/eckit/install/lib/',
    '-L/home/lukas/documents/work/atlas/install/lib/',
    '-Wl,-rpath,/home/lukas/documents/work/atlas/install/lib:/home/lukas/documents/work/eckit/install/lib',
    '-fopenmp'
    ]
cfg['include_dirs'] = [
    '/home/lukas/documents/work/atlas4py/build/_deps/pybind11_fetch-src/include',
    '/home/lukas/documents/work/atlas/install/include',
    '/home/lukas/documents/work/eckit/install/include'
    ]
cfg['libraries'] = ['eckit', 'atlas']

setup_pybind11(cfg)
%>
*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"

namespace py = ::pybind11;
using namespace atlas;
using namespace pybind11::literals;

namespace pybind11 {
namespace detail {
template <>
struct type_caster<array::ArrayStrides>
    : public type_caster<std::vector<array::ArrayStrides::value_type>> {};
template <>
struct type_caster<array::ArrayShape>
    : public type_caster<std::vector<array::ArrayShape::value_type>> {};

}  // namespace detail
}  // namespace pybind11

PYBIND11_MODULE(_atlas4py, m) {
    py::class_<Grid>(m, "Grid")
        .def_property_readonly("name", &StructuredGrid::name)
        .def_property_readonly("uid", &StructuredGrid::uid)
        .def_property_readonly("size", &StructuredGrid::size);

    py::class_<PointLonLat>(m, "PointLonLat")
        .def(py::init([](double lon, double lat) {
                 return PointLonLat({lon, lat});
             }),
             "lon"_a, "lat"_a)
        .def_property_readonly(
            "lon", py::overload_cast<>(&PointLonLat::lon, py::const_))
        .def_property_readonly(
            "lat", py::overload_cast<>(&PointLonLat::lat, py::const_))
        .def("__repr__", [](PointLonLat const& p) {
            return "_atlas4py.PointLonLat(lon=" + std::to_string(p.lon()) +
                   ", lat=" + std::to_string(p.lat()) + ")";
        });
    py::class_<PointXY>(m, "PointXY")
        .def(py::init([](double x, double y) {
                 return PointXY({x, y});
             }),
             "x"_a, "y"_a)
        .def_property_readonly("x",
                               py::overload_cast<>(&PointXY::x, py::const_))
        .def_property_readonly("y",
                               py::overload_cast<>(&PointXY::y, py::const_))
        .def("__repr__", [](PointXY const& p) {
            return "_atlas4py.PointXY(x=" + std::to_string(p.x()) +
                   ", y=" + std::to_string(p.y()) + ")";
        });

    py::class_<StructuredGrid, Grid>(m, "StructuredGrid")
        .def(py::init([](std::string const& s) { return StructuredGrid{s}; }),
             "gridname"_a)
        .def_property_readonly("valid", &StructuredGrid::valid)
        .def_property_readonly("ny", &StructuredGrid::ny)
        .def_property_readonly(
            "nx", py::overload_cast<>(&StructuredGrid::nx, py::const_))
        .def_property_readonly("nxmax", &StructuredGrid::nxmax)
        .def_property_readonly(
            "y", py::overload_cast<>(&StructuredGrid::y, py::const_))
        .def_property_readonly("x", &StructuredGrid::x)
        .def("xy",
             py::overload_cast<idx_t, idx_t>(&StructuredGrid::xy, py::const_),
             "i"_a, "j"_a)
        .def("lonlat",
             py::overload_cast<idx_t, idx_t>(&StructuredGrid::lonlat,
                                             py::const_),
             "i"_a, "j"_a)
        .def_property_readonly("reduced", &StructuredGrid::reduced)
        .def_property_readonly("regular", &StructuredGrid::regular)
        .def_property_readonly("periodic", &StructuredGrid::periodic);

    py::class_<StructuredMeshGenerator>(m, "StructuredMeshGenerator")
        .def(py::init())
        .def("generate", py::overload_cast<Grid const&>(
                             &StructuredMeshGenerator::generate, py::const_));

    py::class_<Mesh>(m, "Mesh")
        .def_property_readonly("grid", &Mesh::grid)
        .def_property("nodes", py::overload_cast<>(&Mesh::nodes, py::const_),
                      py::overload_cast<>(&Mesh::nodes))
        .def_property("edges", py::overload_cast<>(&Mesh::edges, py::const_),
                      py::overload_cast<>(&Mesh::edges))
        .def_property("cells", py::overload_cast<>(&Mesh::cells, py::const_),
                      py::overload_cast<>(&Mesh::cells))
        .def_property_readonly("generated", &Mesh::generated);
    m.def("build_edges", py::overload_cast<Mesh&>(&mesh::actions::build_edges));
    m.def("build_node_to_edge_connectivity",
          py::overload_cast<Mesh&>(
              &mesh::actions::build_node_to_edge_connectivity));

    py::class_<mesh::MultiBlockConnectivity>(m, "MultiBlockConnectivity");

    py::class_<mesh::HybridElements>(m, "HybridElements")
        .def_property_readonly("size", &mesh::HybridElements::size)
        .def("nb_nodes", &mesh::HybridElements::nb_nodes)
        .def("nb_edges", &mesh::HybridElements::nb_edges)
        .def_property_readonly(
            "node_connectivity",
            py::overload_cast<>(&mesh::HybridElements::node_connectivity,
                                py::const_))
        .def_property_readonly(
            "edge_connectivity",
            py::overload_cast<>(&mesh::HybridElements::edge_connectivity,
                                py::const_))
        .def_property_readonly(
            "cell_connectivity",
            py::overload_cast<>(&mesh::HybridElements::cell_connectivity,
                                py::const_));

    auto m_fs = m.def_submodule("functionspace");
    py::class_<FunctionSpace>(m_fs, "FunctionSpace")
        .def_property_readonly("size", &FunctionSpace::size)
        .def_property_readonly("type", &FunctionSpace::type)
        .def("create_field",
             [](FunctionSpace const& fs, std::optional<std::string> const& name,
                std::optional<int> levels, std::string type) {
                 util::Config config;
                 if (name) config = config | option::name(*name);
                 if (levels) config = config | option::levels(*levels);
                 config = config | option::datatype(type);
                 return fs.createField(config);
             },
             "name"_a = std::nullopt, "levels"_a = std::nullopt, "type"_a);
    py::class_<functionspace::EdgeColumns, FunctionSpace>(m_fs, "EdgeColumns")
        .def(py::init(
            [](Mesh const& m) { return functionspace::EdgeColumns(m); }))
        .def_property_readonly("nb_edges",
                               &functionspace::EdgeColumns::nb_edges)
        .def_property_readonly("mesh", &functionspace::EdgeColumns::mesh)
        .def_property_readonly("edges", &functionspace::EdgeColumns::edges)
        .def_property_readonly("valid", &functionspace::EdgeColumns::valid);
    py::class_<functionspace::NodeColumns, FunctionSpace>(m_fs, "NodeColumns")
        .def(py::init(
            [](Mesh const& m) { return functionspace::NodeColumns(m); }))
        .def_property_readonly("nb_nodes",
                               &functionspace::NodeColumns::nb_nodes)
        .def_property_readonly("mesh", &functionspace::NodeColumns::mesh)
        .def_property_readonly("nodes", &functionspace::NodeColumns::nodes)
        .def_property_readonly("valid", &functionspace::NodeColumns::valid);
    py::class_<functionspace::CellColumns, FunctionSpace>(m_fs, "CellColumns")
        .def(py::init(
            [](Mesh const& m) { return functionspace::CellColumns(m); }))
        .def_property_readonly("nb_cells",
                               &functionspace::CellColumns::nb_cells)
        .def_property_readonly("mesh", &functionspace::CellColumns::mesh)
        .def_property_readonly("cells", &functionspace::CellColumns::cells)
        .def_property_readonly("valid", &functionspace::CellColumns::valid);

    py::class_<Field>(m, "Field", py::buffer_protocol())
        .def_property_readonly("name", &Field::name)
        .def_property_readonly("datatype", &Field::datatype)
        .def_property_readonly("strides", &Field::strides)
        .def_property_readonly("shape",
                               py::overload_cast<>(&Field::shape, py::const_))
        .def_property_readonly("size", &Field::size)
        .def_property_readonly("rank", &Field::rank)
        .def_buffer([](Field& f) {
            auto strides = f.strides();
            std::transform(strides.begin(), strides.end(), strides.begin(),
                           [&](auto const& stride) {
                               return stride * f.datatype().size();
                           });
            return py::buffer_info(f.storage(), f.datatype().size(),
                                   py::format_descriptor<double>::format(),
                                   f.rank(), f.shape(), strides);
        });

    py::class_<output::Gmsh>(m, "Gmsh")
        .def(py::init(
                 [](std::string const& path) { return output::Gmsh{path}; }),
             "path"_a)
        .def("write",
             [](output::Gmsh& gmsh, Mesh const& mesh) { gmsh.write(mesh); },
             "mesh"_a)
        .def("write",
             [](output::Gmsh& gmsh, Field const& field) { gmsh.write(field); },
             "field"_a)
        .def("write",
             [](output::Gmsh& gmsh, Field const& field,
                FunctionSpace const& fs) { gmsh.write(field, fs); },
             "field"_a, "functionspace"_a);
}
