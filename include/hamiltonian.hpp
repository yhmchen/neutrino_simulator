#pragma once

#include <complex>
#include <vector>
#include <array>

namespace neutrino {

using ComplexMatrix = std::array<std::array<std::complex<double>, 3>, 3>;

/**
 * 計算有效哈密頓量
 * @param r 從太陽中心的距離 (以太陽半徑為單位)
 * @param E_nu 中微子能量 (MeV)
 * @param r_vals 從檔案讀取的距離值
 * @param Ne_vals 從檔案讀取的電子密度值
 * @return 有效哈密頓量 (eV)
 */
ComplexMatrix hamiltonian(double r, 
                         double E_nu, 
                         const std::vector<double>& r_vals, 
                         const std::vector<double>& Ne_vals);

/**
 * 計算中微子態矢的導數
 * @param t 傳播距離 (以太陽半徑為單位)
 * @param psi 中微子態矢 (3x3 複數矩陣)
 * @param E_nu 中微子能量 (MeV)
 * @param r_vals 從檔案讀取的距離值
 * @param Ne_vals 從檔案讀取的電子密度值
 * @param dpsi_dt 輸出: 中微子態矢的導數
 */
void neutrino_evolution(double t, 
                      const std::array<std::complex<double>, 9>& psi,
                      double E_nu,
                      const std::vector<double>& r_vals,
                      const std::vector<double>& Ne_vals,
                      std::array<std::complex<double>, 9>& dpsi_dt);

} // namespace neutrino