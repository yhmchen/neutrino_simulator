#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <complex>
#include <iomanip>

// 定义复数类型
using Complex = std::complex<double>;

// 定义3x3复数矩阵类型
class Matrix3c {
private:
    std::vector<Complex> data;
public:
    Matrix3c() : data(9, 0.0) {}
    
    Complex& operator()(int i, int j) {
        return data[i * 3 + j];
    }
    
    Complex operator()(int i, int j) const {
        return data[i * 3 + j];
    }
    
    // 矩阵乘法
    Matrix3c operator*(const Matrix3c& other) const {
        Matrix3c result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result(i, j) = 0.0;
                for (int k = 0; k < 3; ++k) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }
        return result;
    }
    
    // 矩阵乘以向量
    std::vector<Complex> operator*(const std::vector<Complex>& vec) const {
        std::vector<Complex> result(3, 0.0);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i] += (*this)(i, j) * vec[j];
            }
        }
        return result;
    }
    
    // 矩阵指数函数 (简化版本，实际应用中可能需要更复杂的实现)
    Matrix3c exp() const {
        Matrix3c result;
        // 对角矩阵的指数很简单
        for (int i = 0; i < 3; ++i) {
            result(i, i) = std::exp((*this)(i, i));
        }
        return result;
    }
};

// 太阳电子密度数据结构
struct SolarDensityData {
    std::vector<double> radius;     // 太阳半径 (以太阳半径为单位)
    std::vector<double> ne;         // 电子数密度 (cm^-3)
};

// 混合参数
struct OscillationParameters {
    double theta12, theta13, theta23;  // 混合角 (弧度)
    double deltaCP;                   // CP违反相位 (弧度)
    double dm21_sq, dm31_sq;          // 质量平方差 (eV^2)
    double energy;                     // 中微子能量 (MeV)
};

// 读取BP2000电子密度数据文件
SolarDensityData readDensityData(const std::string& filename) {
    SolarDensityData data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return data;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // 跳过注释行或非数据行
        if (line.empty() || line[0] == '#' || line[0] == '/') {
            continue;
        }
        
        std::istringstream iss(line);
        double r, ne_value;
        
        // 尝试读取两个数值
        if (iss >> r >> ne_value) {
            data.radius.push_back(r);
            data.ne.push_back(ne_value);
        }
    }
    
    file.close();
    return data;
}

// 计算PMNS矩阵
Matrix3c calculatePMNSMatrix(const OscillationParameters& params) {
    Matrix3c U;
    
    // 简化计算，使用 c_ij = cos(theta_ij), s_ij = sin(theta_ij)
    double c12 = std::cos(params.theta12);
    double s12 = std::sin(params.theta12);
    double c13 = std::cos(params.theta13);
    double s13 = std::sin(params.theta13);
    double c23 = std::cos(params.theta23);
    double s23 = std::sin(params.theta23);
    
    Complex phase = std::exp(Complex(0, -params.deltaCP));
    
    // PMNS矩阵元素
    U(0, 0) = c12 * c13;
    U(0, 1) = s12 * c13;
    U(0, 2) = s13 * std::conj(phase);
    
    U(1, 0) = -s12 * c23 - c12 * s23 * s13 * phase;
    U(1, 1) = c12 * c23 - s12 * s23 * s13 * phase;
    U(1, 2) = s23 * c13;
    
    U(2, 0) = s12 * s23 - c12 * c23 * s13 * phase;
    U(2, 1) = -c12 * s23 - s12 * c23 * s13 * phase;
    U(2, 2) = c23 * c13;
    
    return U;
}

// 计算哈密顿量
Matrix3c calculateHamiltonian(double ne, const OscillationParameters& params) {
    Matrix3c H;
    
    // 真空项
    const double GF = 1.166378e-23; // Fermi常数 (eV^-2)
    double dm21 = params.dm21_sq;
    double dm31 = params.dm31_sq;
    double E = params.energy;
    
    // 物质效应项
    double V = std::sqrt(2) * GF * ne;
    
    // 质量矩阵在物质中 (假设对角基)
    H(0, 0) = 0.0 + V;
    H(1, 1) = dm21 / (2 * E);
    H(2, 2) = dm31 / (2 * E);
    
    return H;
}

// 计算演化矩阵 (在给定位置和步长的情况下)
Matrix3c calculateEvolutionMatrix(double r, double dr, const SolarDensityData& data, const OscillationParameters& params) {
    // 找到最接近的电子密度数据点
    int idx = 0;
    double min_diff = std::abs(data.radius[0] - r);
    
    for (size_t i = 1; i < data.radius.size(); ++i) {
        double diff = std::abs(data.radius[i] - r);
        if (diff < min_diff) {
            min_diff = diff;
            idx = i;
        }
    }
    
    double ne = data.ne[idx];
    
    // 计算该位置的哈密顿量
    Matrix3c H = calculateHamiltonian(ne, params);
    
    // 计算演化矩阵 (简化版本，实际应用可能需要更复杂的数值积分)
    // U(r+dr, r) = exp(-i * H * dr)
    Matrix3c evolution;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            evolution(i, j) = (i == j) ? std::exp(Complex(0, -1) * H(i, i) * dr) : 0.0;
        }
    }
    
    return evolution;
}

// 模拟中微子振荡
std::vector<std::vector<Complex>> simulateOscillation(const SolarDensityData& data, const OscillationParameters& params) {
    // 初始中微子状态 (假设是纯电子中微子)
    std::vector<Complex> neutrino_state = {1.0, 0.0, 0.0};
    
    // 定义步长
    double r_start = data.radius.front(); // 太阳中心
    double r_end = data.radius.back();    // 太阳表面
    int steps = 1000;  // 模拟步数
    double dr = (r_end - r_start) / steps;
    
    // 存储每一步的中微子状态
    std::vector<std::vector<Complex>> states;
    states.push_back(neutrino_state);
    
    // PMNS矩阵
    Matrix3c U = calculatePMNSMatrix(params);
    Matrix3c U_dag;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            U_dag(i, j) = std::conj(U(j, i));
        }
    }
    
    // 从太阳中心到表面的演化
    for (int i = 0; i < steps; ++i) {
        double r_current = r_start + i * dr;
        
        // 计算该步长的演化矩阵
        Matrix3c evolution = calculateEvolutionMatrix(r_current, dr, data, params);
        
        // 更新中微子状态
        neutrino_state = evolution * neutrino_state;
        
        // 保存状态
        states.push_back(neutrino_state);
    }
    
    // 从太阳表面到地球的演化 (真空中)
    // 在真空中，只需考虑标准的振荡公式
    // 近似为自由传播
    
    // 假设太阳到地球的距离 (AU)
    double sunToEarth = 1.0; // 1 AU, 约1.5e11 米
    double R_sun = 6.957e8;  // 太阳半径 (米)
    
    // 太阳到地球的振荡
    // 简化处理，实际可能需要更复杂的计算
    // 在真空中传播
    Matrix3c vacuum_evolution;
    for (int i = 0; i < 3; ++i) {
        vacuum_evolution(i, i) = std::exp(Complex(0, -params.dm21_sq * sunToEarth * R_sun / (2 * params.energy)));
    }
    
    neutrino_state = vacuum_evolution * neutrino_state;
    states.push_back(neutrino_state);
    
    return states;
}

// 计算生存概率
std::vector<double> calculateSurvivalProbabilities(const std::vector<std::vector<Complex>>& states) {
    std::vector<double> probabilities;
    
    for (const auto& state : states) {
        // 计算 |<νe|ν(r)>|^2
        double prob = std::norm(state[0]);
        probabilities.push_back(prob);
    }
    
    return probabilities;
}

int main() {
    // 读取电子密度数据
    std::string filename = "/Users/hcjhuang/Documents/各種資料/neutrino_simulator/BP2000 electron density.txt";
    SolarDensityData data = readDensityData(filename);
    
    if (data.radius.empty()) {
        std::cerr << "无法读取有效的电子密度数据" << std::endl;
        return 1;
    }
    
    std::cout << "成功读取了 " << data.radius.size() << " 个数据点" << std::endl;
    
    // 设置振荡参数 (使用当前最佳拟合值)
    OscillationParameters params;
    params.theta12 = 33.45 * M_PI / 180.0;  // 太阳角
    params.theta13 = 8.62 * M_PI / 180.0;   // 反应堆角
    params.theta23 = 42.1 * M_PI / 180.0;   // 大气角
    params.deltaCP = 230.0 * M_PI / 180.0;  // CP违反相位
    params.dm21_sq = 7.42e-5;              // 太阳质量平方差 (eV^2)
    params.dm31_sq = 2.515e-3;             // 大气质量平方差 (eV^2)
    params.energy = 10.0;                  // 中微子能量 (MeV)
    
    // 模拟中微子振荡
    std::vector<std::vector<Complex>> states = simulateOscillation(data, params);
    
    // 计算生存概率
    std::vector<double> probabilities = calculateSurvivalProbabilities(states);
    
    // 输出结果
    std::cout << "中微子振荡模拟结果：" << std::endl;
    std::cout << "太阳中心电子中微子生存概率: " << probabilities.front() << std::endl;
    std::cout << "太阳表面电子中微子生存概率: " << probabilities[probabilities.size() - 2] << std::endl;
    std::cout << "到达地球的电子中微子生存概率: " << probabilities.back() << std::endl;
    
    // 输出详细结果到文件
    std::ofstream outfile("neutrino_oscillation_results.txt");
    if (outfile.is_open()) {
        outfile << "# 太阳半径 (R_sun)\t电子中微子生存概率" << std::endl;
        double r_start = data.radius.front();
        double r_end = data.radius.back();
        double dr = (r_end - r_start) / (probabilities.size() - 2);
        
        for (size_t i = 0; i < probabilities.size() - 1; ++i) {
            double r = r_start + i * dr;
            outfile << r << "\t" << probabilities[i] << std::endl;
        }
        
        outfile << "# 地球" << "\t" << probabilities.back() << std::endl;
        outfile.close();
        
        std::cout << "结果已保存到 neutrino_oscillation_results.txt" << std::endl;
    } else {
        std::cerr << "无法创建输出文件" << std::endl;
    }
    
    return 0;
}