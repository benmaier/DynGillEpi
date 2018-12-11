/* 
 * The MIT License (MIT)
 * Copyright (c) 2018, Benjamin Maier
 *
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without 
 * restriction, including without limitation the rights to use, 
 * copy, modify, merge, publish, distribute, sublicense, and/or 
 * sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall 
 * be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-
 * INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 */

#include "Utilities.h"
#include "SIS_Poisson_homogeneous.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace std;
namespace py = pybind11;

PYBIND11_PLUGIN(DynGillEpi) {
    py::module m("DynGillEpi", "Module to perform fast flockwork simulations");
    
    m.def("SIS_Poisson_homogeneous", &SIS_Poisson_homogeneous, "Simulate an SIS process on a time-dependent contact list.",
            py::arg("N"),
            py::arg("list_of_contact_lists"),
            py::arg("infection_rate_per_dt"),
            py::arg("recovery_rate_per_dt"),
            py::arg("T_simulation") = 0,
            py::arg("output_time_resolution_in_dt") = 1,
            py::arg("number_of_simulations") = 1,
            py::arg("initial_number_of_infected") = 1,
            py::arg("seed") = 0,
            py::arg("t_infection_start") = 0,
            py::arg("verbose") = false
            );

    py::class_<SI_result>(m,"SI_result")
        .def(py::init<>())
        .def_readwrite("true_I", &SI_result::true_I)
        .def_readwrite("true_SI", &SI_result::true_SI)
        .def_readwrite("true_t", &SI_result::true_t)
        .def_readwrite("I", &SI_result::I)
        .def_readwrite("SI", &SI_result::SI)
        .def_readwrite("hist", &SI_result::hist)
        ;

    return m.ptr();

}
