#include <integrator/integrator.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <iostream>


namespace py = pybind11;


namespace cdriver {


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

} // namespace cdriver

PYBIND11_MODULE(cdriver, m) {
  m.doc() = R"doc(
    The computation engine for cosmic integrator
)doc";
  m.def("sample_from_imf", &cdriver::sample_from_imf, py::arg("n"));
}