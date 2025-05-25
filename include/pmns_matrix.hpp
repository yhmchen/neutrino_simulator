#pragma once

#include <complex>
#include <array>
#include <string>
#include <iostream>

namespace neutrino {

// 使用 array 代替 Eigen 矩陣，以簡化範例
using ComplexMatrix = std::array<std::array<std::complex<double>, 3>, 3>;

void printMatrix(const ComplexMatrix& M, const std::string& name);

/**
 * 計算 PMNS 矩陣
 * @return PMNS 矩陣 (3x3 複數矩陣)
 */
ComplexMatrix get_pmns_matrix();


} // namespace neutrino