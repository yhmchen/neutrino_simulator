#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"
#include "../include/simulator.hpp"
#include "../include/constants.hpp"

namespace py = pybind11;

PYBIND11_MODULE(neutrino_simulator, m) {
    m.doc() = "C++ 實作的中微子振盪模擬器";
    /*
    // 保留原始：自動從檔案讀密度並決定距離
    m.def("simulate_oscillation", &neutrino::simulate_oscillation,
        py::arg("E_nu"),
        py::arg("solar_density_profile_path"),
        py::arg("r_start"),
        py::arg("r_end"),
        py::arg("num_points"),
        "模擬中微子振盪（自動從密度檔案讀取）");
*/
    // 公開新版模擬函數（純粹用 input_r_vals 控制模擬距離點）
    m.def("simulate_custom_rvals", &neutrino::simulate_custom_rvals, 
        py::arg("E_nu"), 
        py::arg("input_r_vals"),
        py::arg("density_r_vals"),
        py::arg("Ne_vals"),
        "模擬中微子振盪（距離與密度由 Python 傳入）");

    
    // 公開一些常數
    m.attr("R_sun") = neutrino::R_sun;
    m.attr("R_earth_orbit") = neutrino::R_earth_orbit;
    
    // 將自定義結構轉換為Python物件
    py::class_<neutrino::SimulationResult>(m, "SimulationResult")
        .def_readonly("r_vals", &neutrino::SimulationResult::r_vals)
        .def_readonly("probs", &neutrino::SimulationResult::probs);
}