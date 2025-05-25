#pragma once

#include <complex>
#include <cmath>

namespace neutrino {

// 物理常數
extern const double hbar;        // 約化普朗克常數 (J·s)
extern const double c;           // 光速 (m/s)
extern const double eV;          // 電子伏特 (J)
extern const double MeV;         // 兆電子伏特 (J)
extern const double GeV;         // 吉電子伏特 (J)
extern const double km;          // 千米 (m)
extern const double cm;          // 厘米 (m)
extern const double N_A;         // 阿伏加德羅常數

// 太陽和地球參數
extern const double R_sun;       // 太陽半徑 (m)
extern const double R_earth_orbit; // 地球軌道半徑 (m)
extern const double G_F;         // 費米耦合常數 (GeV^-2)

// 中微子振盪參數
extern const double theta12;     // 太陽角 (rad)
extern const double theta13;     // 反應堆角 (rad)
extern const double theta23;     // 大氣角 (rad)
extern const double delta_cp;    // CP違反相位 (rad)

// 質量平方差 (eV^2)
extern const double delta_m21_squared; // 太陽中微子質量平方差 (eV^2)
extern const double delta_m31_squared; // 大氣中微子質量平方差 (eV^2)

} // namespace neutrino