#include "../include/electron_density.hpp"
#include "../include/constants.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>


namespace neutrino {

// 新增自定義電子密度函數
double custom_density_function(double r) {
    // 添加平滑因子，避免太陽表面附近急劇變化
    const double surface = 1.0;
    const double smoothing = 0.01;
    if (r > surface - smoothing) {
        double factor = (surface - r) / smoothing;
        return std::pow(10.0, 26.0 - 10* (surface - smoothing)) * factor;
    }
    // n(r) = 10^(-13-4.3r) GeV^3 = 10^(14-4.3r) eV^3
    return std::pow(10.0, 26.0 -  r);
}


bool solar_electron_density(const std::string& filename, 
    std::vector<double>& r_vals, 
    std::vector<double>& Ne_vals,
    const std::string& model = "BP2000") {

    if (model == "custom") {
        // 使用自定義模型
        r_vals.clear();
        Ne_vals.clear();

        // 生成適當數量的點，範圍根據需要調整
        const int num_points = 1000;  // 可以根據需要調整點數
        const double r_min = 0.0;    // 太陽中心
        const double r_max = 1.00135;    // 太陽表面

        for (int i = 0; i < num_points; i++) {
        double r = r_min + i * (r_max - r_min) / (num_points - 1);
        r_vals.push_back(r);

        // 使用自定義密度函數
        Ne_vals.push_back(custom_density_function(r));
    }

    std::cout << "Using custom electron density model: n(r) = 10^(14-4.3r) eV^3" << std::endl;
    std::cout << "Number of r_vals points: " << r_vals.size() << std::endl;

    // 檢查 Ne_vals 是否計算正確
    std::cout << "Sample Ne_vals: ";
    for (size_t i = 0; i < std::min<size_t>(5, Ne_vals.size()); ++i) {
        std::cout << Ne_vals[i] << " ";
    }
    std::cout << "..." << std::endl;

    return true;
  }
  else {
        // 原始的BP2000模型實現
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }

        // 跳過前兩行
        std::string line;
        for (int i = 0; i < 2; i++) {
            if (!std::getline(file, line)) {
                return false;
            }
        }

        // 讀取數據
        r_vals.clear();
        Ne_vals.clear();

        while (std::getline(file, line)) {
        // 排除最後一行及空白行
        if (line.empty() || line.find("BP2000") != std::string::npos) {
            continue;
        }

        std::istringstream iss(line);
        double r, log_Ne_NA;

        if (iss >> r >> log_Ne_NA) {
            r_vals.push_back(r);  //將讀取到的距離r儲存到r_vals
            // 將 log(N_e/N_A) 轉換為 N_e/N_A=n_e
            Ne_vals.push_back(std::pow(10, log_Ne_NA)*N_A); //eV* std::pow(10, 13), std::exp(log_Ne_NA) 2*std::pow(10,-12)
        }
    }
        std::cout << "Using BP2000 electron density model" << std::endl;
        std::cout << "Number of r_vals points: " << r_vals.size() << std::endl;
        // 檢查 Ne_vals 是否讀入正確
        std::cout << "Ne_vals: ";
        for (size_t i = 0; i < std::min<size_t>(5, Ne_vals.size()); ++i) {
            std::cout << Ne_vals[i] << " ";
        }
        std::cout << "..." << std::endl;

        return !r_vals.empty();
    }
}

double electron_potential(double r, 
                         const std::vector<double>& r_vals, 
                         const std::vector<double>& Ne_vals) {
    // 簡單線性插值
    if (r < r_vals.front()) {
        // 外推: 假設在太陽中心之前的密度等於太陽中心的密度
        return std::sqrt(2) * G_F * Ne_vals.front(); // 轉換為MeV
    }
    else if (r > r_vals.back()) {
        // 使用指數衰減外插
        double last_r = r_vals.back();
        double last_Ne = Ne_vals.back();
        double decay_length = 0.1 * r_vals.back(); // 衰減長度可以調整
        double Ne_extrap = last_Ne * std::exp(-(r - last_r) / decay_length);
        return std::sqrt(2) * G_F * Ne_extrap;
    }
    
    
    else {
        // 找到最接近的較小值索引
        auto it = std::lower_bound(r_vals.begin(), r_vals.end(), r);
        size_t idx;
        
        if (it == r_vals.end()) {
            idx = r_vals.size() - 2;
        }
        else if (it == r_vals.begin()) {
            idx = 0;
        }
        else {
            idx = std::distance(r_vals.begin(), it) - 1;
        }
        
        // 線性插值
        double r1 = r_vals[idx];
        double r2 = r_vals[idx + 1];
        double Ne1 = Ne_vals[idx];
        double Ne2 = Ne_vals[idx + 1];
        
        double Ne_interp = Ne1 + (r - r1) * (Ne2 - Ne1) / (r2 - r1);
        return std::sqrt(2) * G_F * Ne_interp ; // 轉換為 MeV, eV:1e6
    }
}

} // namespace neutrino