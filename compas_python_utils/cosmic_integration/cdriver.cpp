#include <integrator/integrator.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;



py::array_t<double> sample_from_imf(int n)
{

  double *samples= integrator::kroupa_imf::sample_n_masses(n);

    // Create a NumPy array from the double* array
    py::array_t<double> array(
        { n },  // Shape
        { sizeof(double) },  // Strides
        samples,  // Data pointer
        py::capsule(samples, [](void* ptr) { delete[] static_cast<double*>(ptr); })  // Capsule to free memory
    );

    return array;
}



// Binding code for python
PYBIND11_MODULE(cdriver, m) {
  m.doc() = R"doc(
    The computation engine for cosmic integrator
)doc";
  m.def("sample_from_imf", &sample_from_imf, py::arg("n"));
  m.def("compute_star_forming_mass_per_binary", &integrator::kroupa_imf::compute_star_forming_mass_per_binary,
  py::arg("binaryFraction"), py::arg("Mlower"), py::arg("Mupper"), py::arg("m2_min"), py::arg("n"));
}