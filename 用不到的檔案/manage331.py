import numpy as np
import matplotlib.pyplot as plt
import time
import multiprocessing as mp
import os
import neutrino_simulator  # 導入 C++ 擴展模組
import subprocess

def run_simulation():
    # 執行 C++ 程式（假設編譯後的可執行檔名為 "neutrino_simulator"）
    subprocess.run(["./neutrino_simulator"])
# 從文件讀取 C++ 輸出的數據

def plot_results():
    # 從文件讀取 C++ 輸出的數據（太陽內部）
    data_sun = np.loadtxt("output_data.txt")
    r_vals_sun = data_sun[:, 0]          # 距離 (R/R_sun)
    probs_e_to_e_sun = data_sun[:, 1]    # P(νe -> νe)
    probs_e_to_mu_sun = data_sun[:, 2]   # P(νe -> νμ)
    probs_e_to_tau_sun = data_sun[:, 3]  # P(νe -> ντ)

    # 太陽外部距離點
    r_outer = np.linspace(1.00135, 1.496e11, 1000)  # 距離 (m)
    r_outer_Rsun = r_outer / 6.96e8  # 將距離轉換為 R/R_sun

    # 假設太陽外部的機率保持不變（這裡需要根據你的實際需求進行修改）
    probs_e_to_e_outer = np.ones_like(r_outer_Rsun)  # P(νe -> νe)
    probs_e_to_mu_outer = np.zeros_like(r_outer_Rsun)  # P(νe -> νμ)
    probs_e_to_tau_outer = np.zeros_like(r_outer_Rsun)  # P(νe -> ντ)
    

    # 繪製太陽內部的生存機率
    plt.figure(figsize=(10, 6))
    plt.plot(r_vals_sun, probs_e_to_e_sun, label='P(νe -> νe)')
    plt.plot(r_vals_sun, probs_e_to_mu_sun, label='P(νe -> νμ)')
    plt.plot(r_vals_sun, probs_e_to_tau_sun, label='P(νe -> ντ)')

    plt.xlabel('Distance (R/R_sun)')
    plt.ylabel('Survival Probability')
    plt.title('Neutrino Oscillation Probabilities Inside sun')
    plt.legend()
    plt.grid(True)
    plt.show()

    # 繪製太陽外部的生存機率
    plt.figure(figsize=(10, 6))
    plt.plot(r_outer_Rsun, probs_e_to_e_outer, label='P(νe -> νe)')
    plt.plot(r_outer_Rsun, probs_e_to_mu_outer, label='P(νe -> νμ)')
    plt.plot(r_outer_Rsun, probs_e_to_tau_outer, label='P(νe -> ντ)')

    plt.xlabel('Distance (R/R_sun)')
    plt.ylabel('Survival Probability')
    plt.title('Neutrino Oscillation Probabilities Outside Sun')
    plt.legend()
    plt.grid(True)
    plt.show()


#try:
#    import neutrino_simulator
#    print("C++ 模組成功導入")
#except ImportError as e:
#    print(f"C++ 模組導入失敗：{e}")
#    exit()  # 如果導入失敗，直接退出程式

import sys
sys.path.append('/Users/hcjhuang/Documents/各種資料/neutrino_simulator')

def plot_neutrino_oscillations(valid_results):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # 使用更多的绘图样式和颜色
    colors = plt.cm.viridis(np.linspace(0, 1, len(valid_results)))
    #line_styles = ['-', '--', '-.', ':']  # 提供一些线型
    #markers = ['o', 's', '^', 'D']      # 提供一些标记
    
    # 全距离图
    for (E_nu, r_vals, probs), color in zip(valid_results, colors):
        # P(νe -> νe)
        survival_probs = [prob[0] for prob in probs]
        ax1.plot(r_vals, survival_probs,label=f'{E_nu:.2f} MeV (νe -> νe)',color='red',markersize=4,alpha=0.7)
        
        # P(νe -> νμ)
        oscillation_probs_mu = [prob[1] for prob in probs]
        ax1.plot(r_vals, oscillation_probs_mu,label=f'{E_nu:.2f} MeV (νe -> νμ)',color='green',markersize=4,alpha=0.7)
        
        # P(νe -> ντ)
        oscillation_probs_tau = [prob[2] for prob in probs]
        ax1.plot(r_vals, oscillation_probs_tau,  label=f'{E_nu:.2f} MeV (νe -> ντ)', color='blue',markersize=4,alpha=0.7)
    
    ax1.set_xlabel('Distance from the center of the sun (solar radius)', fontsize=10)
    ax1.set_ylabel('Neutrino Oscillation Probability', fontsize=10)
    ax1.set_title('Electron neutrino survival probability at different energies', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()
        
    # 太阳内部图
    for (E_nu, r_vals, probs), color in zip(valid_results, colors):
        sun_indices = [i for i, r in enumerate(r_vals) if r <= 1.0]
        if sun_indices:
                r_sun_only = [r_vals[i] for i in sun_indices]
                
                # P(νe -> νe)
                probs_sun_only_ee = [probs[i][0] for i in sun_indices]
                ax2.plot(r_sun_only, probs_sun_only_ee,label=f'{E_nu:.2f} MeV (νe -> νe)',color='red',markersize=4,alpha=0.7)

                # P(νe -> νμ)
                probs_sun_only_mu = [probs[i][1] for i in sun_indices]
                ax2.plot(r_sun_only, probs_sun_only_mu,label=f'{E_nu:.2f} MeV (νe -> νμ)', color='green',markersize=4, alpha=0.7)

                # P(νe -> ντ)
                probs_sun_only_tau = [probs[i][2] for i in sun_indices]
                ax2.plot(r_sun_only, probs_sun_only_tau,label=f'{E_nu:.2f} MeV (νe -> ντ)', color='blue',markersize=4, alpha=0.7)
        
        
        ax2.set_xlabel('Distance from the center of the sun (solar radius)', fontsize=10)
        ax2.set_ylabel('Neutrino Oscillation Probability', fontsize=10)
        ax2.set_title('Neutrino Oscillation Probabilities Inside the Sun', fontsize=12)
        ax2.grid(True, linestyle='--', alpha=0.7)
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig("neutrino_oscillation.png", dpi=300, bbox_inches='tight')
        plt.show()
        print("圖形已保存為 neutrino_oscillation.png")

def simulate_single_energy(args):
    print(f"simulate_single_energy()函數被呼叫，參數為：{args}")
    """為單一能量模擬中微子振盪"""
    E_nu, density_file, r_vals_space = args
    print(f"simulate_single_energy()函數被呼叫，參數為：{args}")
    # 將 E_nu 轉換為 Python 的 float 類型
    E_nu = float(E_nu)

    # 使用 C++ 實作的模擬器
    r_start = r_vals_space[0]
    r_end = r_vals_space[-1]
    num_points = len(r_vals_space)
    
    start_time = time.time()
    try:
        result = neutrino_simulator.simulate_oscillation(
            E_nu, 
            density_file, 
            r_start, 
            r_end, 
            num_points
        )
        # 新增：印出詳細的模擬結果
        print("Simulation Result:")
        print(f"  Energy: {E_nu:.2f} MeV")
        print(f"  Distance range: {r_start:.2f} to {r_end:.2f}")
        print(f"  Number of points: {num_points}")
        print(f"  R values: {result.r_vals}")
        print(f"  Probabilities: {result.probs}")  # 印出生存機率
        
        duration = time.time() - start_time
        print(f"Energy {E_nu:.2f} MeV simulation completed in {duration:.2f} seconds")
        return E_nu, result.r_vals, result.probs
    except Exception as e:
        print(f"Energy {E_nu:.2f} MeV simulation failed: {str(e)}")
        return E_nu, None, None


def main():
    # 讀取太陽電子密度數據文件路徑
    base_dir = os.path.dirname(os.path.abspath(__file__))
    density_file = input(f"請輸入太陽電子密度數據文件的路徑 (默認: {base_dir}/BP2000 electron density.txt): ")
    if not density_file:
        density_file = os.path.join(base_dir, "BP2000 electron density.txt")
    
    # 確認文件存在
    if not os.path.exists(density_file):
        print(f"找不到文件：{density_file}")
        return
    
    # 設定中微子能量範圍 (MeV)
    energies_input = input("請輸入中微子能量，以逗號分隔 (默認: 1.0,5.0,10.0 MeV): ")
    if energies_input:
        try:
            energies = np.array([float(e.strip()) for e in energies_input.split(",")])
        except ValueError:
            print("能量輸入格式錯誤，使用默認值")
            energies = np.array([1.0, 5.0, 10.0])
    else:
        energies = np.array([1.0, 5.0, 10.0])
    
    # 設定距離範圍
    r_sun = np.linspace(0, 1.00135, 500)  # 太陽內部
    
    # 太陽表面到地球軌道的距離
    r_earth = neutrino_simulator.R_earth_orbit / neutrino_simulator.R_sun
    r_outer = np.linspace(1.00135, r_earth, 1000)
    
    # 合併距離點，並確保數值唯一
    r_vals_space = np.concatenate([r_sun, r_outer])
    r_vals_space = np.sort(np.unique(r_vals_space))
    
    print(f"模擬 {len(energies)} 個不同能量值的中微子振盪")
    print(f"距離範圍: 0 到 {r_earth:.2f} 太陽半徑")
    
    # 準備多處理參數
    num_cores = max(1, mp.cpu_count() - 1)
    print(f"使用 {num_cores} 個 CPU 核心進行並行計算")
    
    # 打包參數
    simulation_args = [(E_nu, density_file, r_vals_space) for E_nu in energies]
    
    # 使用多處理加速計算
    with mp.Pool(num_cores) as pool:
        results = pool.map(simulate_single_energy, simulation_args)
    print(f"simulation_args 列表：{simulation_args}")

    with mp.Pool(num_cores) as pool:
        results = pool.map(simulate_single_energy, simulation_args)
    print(f"results 列表：{results}")
    # 繪圖
    valid_results = [(E_nu, r_vals, probs) for E_nu, r_vals, probs in results if r_vals is not None]
    
    if valid_results:
        plot_neutrino_oscillations(valid_results) 
    else:
        print("所有模擬都失敗了，請檢查錯誤訊息。")

if __name__ == "__main__":
    main()