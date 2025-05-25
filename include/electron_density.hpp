#pragma once

#include <vector>
#include <string>

namespace neutrino {

/**
 * 從檔案讀取太陽電子密度數據
 * @param filename 密度數據檔案路徑
 * @param r_vals 輸出: 從太陽中心的距離 (以太陽半徑為單位)
 * @param Ne_vals 輸出: 對應的電子數密度值
 * @return 是否成功讀取數據
 */
bool solar_electron_density(const std::string& filename, 
                           std::vector<double>& r_vals, 
                           std::vector<double>& Ne_vals);

/**
 * 計算給定位置的電子勢能
 * @param r 從太陽中心的距離 (以太陽半徑為單位)
 * @param r_vals 從檔案讀取的距離值
 * @param Ne_vals 從檔案讀取的電子密度值
 * @return 電子勢能 (eV)
 */
double electron_potential(double r, 
                         const std::vector<double>& r_vals, 
                         const std::vector<double>& Ne_vals);

} // namespace neutrino