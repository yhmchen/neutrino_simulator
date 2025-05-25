cmake_minimum_required(VERSION 3.10)
project(neutrino_simulator)

# 設定 C++ 標準
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 尋找 Python 解釋器
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# 下載和配置 pybind11
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.10.4
)
FetchContent_MakeAvailable(pybind11)

# 包含頭文件目錄
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${pybind11_INCLUDE_DIRS}
    )

# 收集所有源文件
set(SOURCES
    src/constants.cpp
    src/electron_density.cpp
    src/pmns_matrix.cpp
    src/hamiltonian.cpp
    src/simulator.cpp
)

# 創建靜態庫
add_library(neutrino_lib STATIC ${SOURCES})

# 創建 Python 模組
pybind11_add_module(neutrino_simulator MODULE src/binding.cpp)
target_link_libraries(neutrino_simulator PRIVATE neutrino_lib)

# 確保 PYTHON_SITE_PACKAGES 變數已正確設定
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path")
elseif(UNIX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
endif()

install(TARGETS neutrino_simulator
        #LIBRARY 
        DESTINATION "${PYTHON_SITE_PACKAGES}"
        #RUNTIME DESTINATION "${PYTHON_SITE_PACKAGES}"
        )


# 設置輸出目錄
set(CMAKE_BUILD_TYPE Release) # 或 Debug
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)