import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import time
import os
import math
import neutrino_simulator  # 導入 C++ 擴展模組

# === 常數 ===
R_earth = neutrino_simulator.R_earth_orbit / neutrino_simulator.R_sun
density_file = "/Users/hcjhuang/Documents/各種資料/neutrino_simulator/BP2000 electron density.txt"
energies = [11.0e6]
#energies = [1e6 * x for x in range(1, 21)]  # 1~20 MeV
#energies = [1.0e6, 5.0e6, 10.0e6]  # MeV，修改這裡可以模擬多個能量
#energies = np.logspace(5.0, 6.6, 100)  # 300k eV ~ 100 MeV


# === 輔助函數：讀取密度檔案 ===
def load_density_profile(filename):
    r_vals, Ne_vals = [], []
    with open(filename, 'r') as f:
        for line in f.readlines()[2:]:  # 跳過前兩行
            if not line.strip() or "BP2000" in line:
                continue
            parts = line.strip().split()
            try:
                r = float(parts[0])
                log_Ne = float(parts[1])
                r_vals.append(r)
                Ne_vals.append(10**log_Ne*6.02214076e23)
            except:
                continue
    return np.array(r_vals), np.array(Ne_vals)


# === 新增：自定義密度模型函數 ===
def custom_density_model(r_vals):
    """
    計算自定義密度模型: n(r) = 10^(14-4.3r) eV^3
    """
    return np.power(10.0, 14.0 - 4.3 * r_vals)

# === 新增：獲取不同模型的密度值 ===
def get_density_values(r_vals, model="BP2000", density_file=density_file):
    """
    根據選擇的模型返回對應的密度值
    
    參數:
        r_vals: 距離點數組
        model: 模型名稱 ("BP2000" 或 "custom")
        density_file: BP2000模型的文件路徑
        
    返回:
        (密度點的r值, 對應的Ne值)
    """
    if model == "custom":
        # 使用自定義模型，直接計算每個r點的密度值
        density_r_vals = np.linspace(0.0, 1.00135, 1000)  # 產生均勻分布的點來代表密度分佈
        Ne_vals = custom_density_model(density_r_vals)
        print(f"[INFO] 使用自定義密度模型: n(r) = 10^(14-4.3r) eV^3")
        print(f"[INFO] 太陽中心 Ne = {Ne_vals[0]:.2e}, 太陽表面 Ne = {Ne_vals[-1]:.2e}")
        return density_r_vals, Ne_vals
    else:
        # 使用BP2000模型，從文件讀取
        density_r_vals, Ne_vals = load_density_profile(density_file)
        print(f"[INFO] 使用BP2000密度模型，共{len(density_r_vals)}個點")
        print(f"[INFO] 太陽中心 Ne = {Ne_vals[0]:.2e}, 太陽表面 Ne = {Ne_vals[-1]:.2e}")
        return density_r_vals, Ne_vals


import traceback
from concurrent.futures import ProcessPoolExecutor, TimeoutError, as_completed
# === 輔助函數：模擬單一能量 ===
def simulate_single_energy(args):
    E_nu, r_vals, density_r_vals, Ne_vals = args
    try:
        print(f"[DEBUG] simulate E = {E_nu:.2e}, total r = {len(r_vals)}")

        # 只跑前幾個點先測試是否初始化就爆
        #test_r_vals = r_vals[:20]

        print("[STEP 1] 呼叫 simulate_custom_rvals() ...")
        result = neutrino_simulator.simulate_custom_rvals(
            float(E_nu),
            r_vals.tolist(),
            density_r_vals.tolist(),
            Ne_vals.tolist()
        )
        print("[STEP 2] 成功回傳 simulate_custom_rvals() 結果")

         # ✅ 立即確認是否回傳有效
        if not result.r_vals or not result.probs:
            print(f"❌ 回傳空值：r_vals 或 probs 無效，跳過 E = {E_nu}")
            return E_nu, None, None

        if len(result.r_vals) != len(result.probs):
            print(f"❌ r_vals 和 probs 長度不符，跳過 E = {E_nu}")
            return E_nu, None, None

        print(f"✅ r[0] = {result.r_vals[0]:.3f}, r[-1] = {result.r_vals[-1]:.3f}")
        print(f"✅ P_ee[0] = {result.probs[0][0]:.3f}, P_ee[-1] = {result.probs[-1][0]:.3f}")
        
        # Step 2 後面加上自動比對
        #cpp_pees = [float(line.split("P_ee =")[1].strip())
            #for line in open("cpp_debug.log") if "r =" in line and "P_ee =" in line]
        import re

        # Step 2 後面加上自動比對
        print("\n🔍 自動比對 C++ 和 Python 的 P_ee 差異：")

        # ✅ 改成用正規表達式解析 "r = ..., P_ee = ..."
        cpp_pees = []
        pattern = re.compile(r"r\s*=\s*([\d\.eE\+\-]+),\s*P_ee\s*=\s*([\d\.eE\+\-]+)")
        try:    
            # 檢查文件是否存在
            with open("cpp_debug.log", "r") as f:
                # 讀取格式: r值 P_ee值
                for line in f:
                    parts = line.strip().split()
                    if len(parts) == 2:
                        try:
                            r_val = float(parts[0])
                            pee_val = float(parts[1])
                            cpp_pees.append((r_val, pee_val))
                        except ValueError:
                            pass  # 忽略不能轉換為浮點數的行
        except Exception as e:
            print(f"❌ 讀取 cpp_debug.log 時發生錯誤: {e}")


        py_pees = []
        if result.r_vals and result.probs and len(result.r_vals) == len(result.probs):
            py_pees = [(r, row[0]) for r, row in zip(result.r_vals, result.probs)]
        
        if not py_pees:
            print("❌ 無法從結果中獲取有效的 Python P_ee 數據！")
            return E_nu, result.r_vals, result.probs  # 仍然返回可用數據
        
        r_py = np.array([x[0] for x in py_pees])
        pee_py = np.array([x[1] for x in py_pees])

        if not cpp_pees:
            print("❌ 沒有從 cpp_debug.log 讀到任何 P_ee 資料！")
        else:
            r_cpp = np.array([x[0] for x in cpp_pees])
            pee_cpp = np.array([x[1] for x in cpp_pees])

            min_len = min(len(pee_cpp), len(pee_py))
            print(f"C++ 數據點數: {len(pee_cpp)}, Python 數據點數: {len(pee_py)}")
            
            diff_count = 0
            for i in range(min_len):
                diff = abs(pee_cpp[i] - pee_py[i])
                if diff > 1e-6:
                    diff_count += 1
                    if diff_count <= 5:  # 只顯示前5個差異
                        print(f"⚠️ 差異超過：i = {i}, C++ = {pee_cpp[i]:.6f}, Python = {pee_py[i]:.6f}, 差異 = {diff:.2e}")
            
            if diff_count > 5:
                print(f"... 還有 {diff_count-5} 個差異點未顯示")

            if len(r_cpp) > 0:
                print(f"✅ r[0] = {r_cpp[0]:.3f}, P_ee[0] = {pee_cpp[0]:.3f}")
                print(f"✅ r[-1] = {r_cpp[-1]:.3f}, P_ee[-1] = {pee_cpp[-1]:.3f}")

        '''for i in range(min(len(cpp_pees), len(py_pees))):
            r_cpp, pee_cpp = cpp_pees[i]
            r_py, pee_py = py_pees[i]
            diff = abs(pee_cpp - pee_py)
            if diff > 1e-4:
                print(f"⚠️ 差異超過：i={i}, r={r_py:.4f}, C++={pee_cpp:.4f}, Python={pee_py:.4f}, 差={diff:.2e}")

        for i in range(min(len(cpp_pees), len(py_pees))):
            diff = abs(cpp_pees[i] - py_pees[i])
            if diff > 1e-6:
                print(f"⚠️ 差異超過：i = {i}, C++ = {cpp_pees[i]:.6f}, Python = {py_pees[i]:.6f}, 差異 = {diff:.2e}")'''


        # ✅ 插入 debug 輸出
        print(f"\n🌟 E_nu = {E_nu} eV")
        sample_size = min(10, len(result.r_vals))
        sample_indices = [i * len(result.r_vals) // sample_size for i in range(sample_size)]
        
        for i in sample_indices:  # 只印出部分樣本點
            r = result.r_vals[i]
            p = result.probs[i]
            print(f"  r = {r:.3f}, P_ee = {p[0]:.3f}, P_emu = {p[1]:.3f}, P_etau = {p[2]:.3f}")

        return E_nu, result.r_vals, result.probs
    except Exception as e:
        print(f"\n❌ simulate_single_energy() 發生錯誤！E = {E_nu:.2e}")
        traceback.print_exc()
        return E_nu, None, None


# === Step 1.5: 自動共振點細化 ===
def generate_auto_refined_rvals(energy, delta_m21_squared, density_r_vals, Ne_vals):
    """
    根據 energy 和 delta_m2 自動找共振點，並加密共振區附近的 r 值分布
    """
    G_F = 1.166e-23
    delta_m21_squared = 7.42e-5
    energy = energies[0]  # 拿第一個能量
    delta = delta_m21_squared / (2 * energy) # 共振條件 Ve = Δ

    # 把 Ve 對應出來
    Ve_vals = np.sqrt(2) * G_F * Ne_vals

    # 找出 Ve 最接近 delta 的位置
    idx_res = np.argmin(np.abs(Ve_vals - delta))
    r_res = density_r_vals[idx_res]
    print(f"[INFO] E = {energy:.2e} eV 時，共振點大約在 r = {r_res:.4f} R_sun")

    # 給定 margin，細化這段區域
    margin = 0.05
    r_lo = max(0.0, r_res - margin)
    r_hi = min(1.0, r_res + margin)
    r_refined = np.linspace(r_lo, r_hi, 2000)  # 這段加密

    return r_refined


'''# === 額外圖：E vs P_ee（在地球觀測點）===
def plot_energy_vs_pee(results):
    print("📊 [DEBUG] 開始畫 E vs P_ee 圖...")
    
    Es = []
    P_ees = []

    for (E_nu, r_vals, probs) in results:
        r_vals = np.array(r_vals)
        probs = np.array(probs)

        if len(r_vals) == 0 or probs is None:
            print(f"[⚠️] 忽略 E = {E_nu}, 資料為空")
            continue

        # 找最遠點（代表地球）
        P_ee_at_surface = probs[-1][0]  # 最後一個 r 的電子微中子機率
        Es.append(E_nu / 1e6)  # 轉成 MeV
        P_ees.append(P_ee_at_surface)
    if not Es:
        print("❌ 沒有有效的能量資料可畫 E vs P_ee 圖！")
        return

    print(f"✅ 有 {len(Es)} 筆資料要畫 E vs P_ee")

    plt.figure(figsize=(8, 5))
    plt.plot(Es, P_ees, marker='o', linestyle='-', color='crimson')
    plt.xscale("log")
    plt.xlabel("Neutrino Energy (MeV)")
    plt.ylabel("P(νe → νe) at Earth")
    plt.title("Electron Neutrino Survival Probability vs Energy")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("Pee_vs_Energy.png", dpi=300)
    path = os.path.abspath("neutrino_oscillation.png")
    print(f"📁 圖檔儲存於：{path}")
    #plt.show()
    print("✅ 圖片已儲存為 Pee_vs_Energy.png")'''


# === 主程式入口 ===
def main(density_model="BP2000"):
    open("cpp_debug.log", "w").close()  # 清空之前的 log

    # Step 1: 讀密度檔案
    density_r_vals, Ne_vals = get_density_values(None, model=density_model)
    print(f"使用 {density_model} 模型，讀入 {len(density_r_vals)} 個密度點")

    # ✅ 加入 sanity check
    if len(Ne_vals) != len(density_r_vals):
        print("❌ 錯誤：Ne_vals 和 density_r_vals 長度不符！")
        return
    '''
    # ✅ 插入 r = 0.0 點，對應的 Ne 就用 Ne_vals[0]（中心密度）
    if density_r_vals[0] > 0.0:
        density_r_vals = np.insert(density_r_vals, 0, 0.0)
        Ne_vals = np.insert(Ne_vals, 0, Ne_vals[0])  # 等於中心密度

    # sanity check
    print(f"最小密度點 r = {density_r_vals[0]}，對應 Ne = {Ne_vals[0]:.2e}")'''
    print(f"[DEBUG] 太陽中心 Ne = {Ne_vals[0]:.2e} cm⁻³")
    G_F = 1.166e-23
    Ve = np.sqrt(2) *G_F * Ne_vals[0]  # 單位：eV
    print(f"[DEBUG] 對應 Ve = {Ve:.2e} eV")
    delta_m21_squared = 7.42e-5

    # Step 2: 自訂距離點：增加太陽內部的密集取樣點（避免錯過共振），外部均勻取樣再 subsample
    r_core_all = []

    for E in energies:
        r_res_part = generate_auto_refined_rvals(E, delta_m21_squared, density_r_vals, Ne_vals)
        r_core_all.append(r_res_part)

    r_core_combined = np.concatenate(r_core_all)
    r_core_combined = np.concatenate([r_core_combined, density_r_vals])
    r_core_combined = np.concatenate([r_core_combined, np.linspace(0.0, 1.0, 100)])  # 保底：整個內部也取樣一點

    r_outer_full = np.linspace(1.00135, R_earth, 10000)
    step = 5  # 每 10 點取一個點，可調整
    r_outer_sub = r_outer_full[::step]

    r_vals_space = np.concatenate([r_core_combined,density_r_vals, r_outer_sub])

    
    r_vals_space = np.sort(np.unique(r_vals_space))
    print(f"模擬距離點總數：{len(r_vals_space)}")

    # Step 3: 並行模擬
    num_cores = max(1, mp.cpu_count() - 1)
    simulation_args = [(E, r_vals_space, density_r_vals, Ne_vals) for E in energies]

    results=[]
    failed=[]

    #with mp.Pool(num_cores) as pool:
    #    results = pool.map(simulate_single_energy, simulation_args)

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        future_to_energy = {
            executor.submit(simulate_single_energy, args): args[0]
            for args in simulation_args
        }
    '''
    for args in simulation_args:
        E = args[0]
        print(f"\n🔍 開始測試 E = {E:.2e}")
        try:
            result = simulate_single_energy(args)
            if result[1] is not None:
                print(f"✅ E = {E:.2e} 成功")
                results.append(result)
            else:
                print(f"⚠️ E = {E:.2e} 回傳空值")
                failed.append(E)
        except Exception as e:
            print(f"❌ E = {E:.2e} 崩潰，錯誤：{e}")
            failed.append(E)
        '''
    for future in as_completed(future_to_energy):
            E = future_to_energy[future]
            try:
                result = future.result(timeout=30)  # 每筆最多 30 秒
                if result[1] is not None:  # 確認不是 None
                    print(f"✅ E = {E:.1e} 完成")
                    results.append(result)
                else:
                    print(f"⚠️ E = {E:.1e} 回傳空值")
            except TimeoutError:
                print(f"⏰ E = {E:.1e} 超時跳過")
                failed.append(E)
            except Exception as e:
                print(f"❌ E = {E:.1e} 失敗：{e}")
                failed.append(E)
            
    print(f"\n✅ 成功 {len(results)} 筆，❌ 失敗 {len(failed)} 筆")

    # Step 6: 畫圖（如果有成功結果）
    if results:
        print("[DEBUG] ✅ 開始畫距離圖...")
        plot_neutrino_oscillations(results)
        #print("[DEBUG] ✅ 開始畫 E vs P_ee 圖...")
        #plot_energy_vs_pee(results)

# === 畫圖 ===
def plot_neutrino_oscillations(valid_results, model_name="BP2000"):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=False)

    colors = plt.cm.viridis(np.linspace(0, 1, len(valid_results)))

    model_title= "BP2000" if model_name == "BP2000" else "Custom n(r)=10^(14-4.3r) eV^3"

    # === [1] 全距離圖（太陽中心 → 地球） ===
    for (E_nu, r_vals, probs), color in zip(valid_results, colors):
        probs = np.array(probs)
        ax1.plot(r_vals, probs[:, 0], color='red', label=f"{E_nu:.2f} eV (νe → νe)", alpha=0.7)
        ax1.plot(r_vals, probs[:, 1], color='green', label=f"{E_nu:.2f} eV (νe → νμ)", alpha=0.7)
        ax1.plot(r_vals, probs[:, 2], color='blue', label=f"{E_nu:.2f} eV (νe → ντ)", alpha=0.7)

    ax1.set_xlabel('Distance from the Sun Center to Earth (in $R_\\odot$)', fontsize=10)
    ax1.set_ylabel('Oscillation Probability', fontsize=10)
    ax1.set_title('Neutrino Oscillation: Full Distance- {model_title}', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='upper right', fontsize=8)

    # === [2] 太陽內部圖（只看 r <= 1.0 R_sun）===
    for (E_nu, r_vals, probs), color in zip(valid_results, colors):
        probs = np.array(probs)
        r_vals = np.array(r_vals)
        for i in range(5):
            print(f"[DEBUG] r = {r_vals[i]:.3f}, P_ee = {probs[i][0]:.3f}")
        mask_inside_sun = r_vals <= 1.0
        r_sun = r_vals[mask_inside_sun]
        probs_sun = probs[mask_inside_sun]

        ax2.plot(r_sun, probs_sun[:, 0], color='red', label=f"{E_nu:.2f} MeV (νe → νe)", alpha=0.7)
        ax2.plot(r_sun, probs_sun[:, 1], color='green', label=f"{E_nu:.2f} MeV (νe → νμ)", alpha=0.7)
        ax2.plot(r_sun, probs_sun[:, 2], color='blue', label=f"{E_nu:.2f} MeV (νe → ντ)", alpha=0.7)

    ax2.set_xlabel('Distance Inside the Sun (in $R_\\odot$)', fontsize=10)
    ax2.set_ylabel('Oscillation Probability', fontsize=10)
    ax2.set_title('Neutrino Oscillation Inside the Sun- {model_title}', fontsize=12)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    filename = f"neutrino_oscillation_{model_name}.png"
    plt.savefig("neutrino_oscillation.png", dpi=300, bbox_inches='tight')
    path = os.path.abspath(filename)
    print(f"📁 圖檔儲存於：{path}")
    #plt.show()
    print("✅ 圖片已儲存為 neutrino_oscillation.png")



if __name__ == "__main__":
    # 您可以在這裡修改使用哪個模型
    # 選項: "BP2000" 或 "custom"
    density_model = "BP2000"  # 可以改成 "custom" 來使用自定義模型
    
    # 如果想從命令行參數取得模型名稱
    import sys
    if len(sys.argv) > 1:
        model_arg = sys.argv[1].lower()
        if model_arg in ["bp2000", "custom"]:
            density_model = "BP2000" if model_arg == "bp2000" else "custom"
            
    main(density_model)