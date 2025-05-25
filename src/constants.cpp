#include "../include/constants.hpp"
#include <cmath>

namespace neutrino {

// 物理常數
const double hbar = 6.582119569e-16;  // 約化普朗克常數 (eV·s)
const double c = 299792458;           // 光速 (m/s)
const double eV = 1.602176634e-19;    // 電子伏特 (J)
const double MeV = eV * 1e6;          // 兆電子伏特 (J)
const double GeV = eV * 1e9;          // 吉電子伏特 (J)
const double km = 1000.0;             // 千米 (m)
const double cm = 0.01;               // 厘米 (m)
const double N_A = 6.02214076e23;     // 阿伏加德羅常數

// 太陽和地球參數
const double R_sun = 6.957e8;        // 太陽半徑 (m)
const double R_earth_orbit = 1.496e11; // 地球軌道半徑 (m)
const double G_F = 1.166e-23 ; // 費米耦合常數 (eV^-2)

// 中微子振盪參數 (使用當前最佳擬合值)
const double theta12 = 33.45* M_PI / 180.0 ;  // 太陽角 (rad)
const double theta13 = 8.62 * M_PI / 180.0;   // 反應堆角 (rad)
const double theta23 = 42.1* M_PI / 180.0 ;   // 大氣角 (rad)
const double delta_cp = 230.0 * M_PI / 180.0; // CP違反相位 (rad)

// 質量平方差 (eV^2)
const double delta_m21_squared = 7.42e-5;  // 太陽中微子質量平方差 (eV^2)
const double delta_m31_squared = 2.51e-3;  // 大氣中微子質量平方差 (eV^2)

//python setup.py build_ext --inplace

 
} // namespace neutrino