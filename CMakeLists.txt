# CMakeList.txt: alphaHiASM_cmake 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("alphaHiASM_cmake")

# 将源代码添加到此项目的可执行文件。
add_executable (alphaHiASM   "assembly.cpp" "overlap.cpp" "utility.cpp" "assembly.h" "overlap.h" "utility.h" "memory_info.h" "main.cpp")

# TODO: 如有需要，请添加测试并安装目标。
include_directories ("./seqan/include" "./boost")
target_link_libraries (alphaHiASM "pthread")