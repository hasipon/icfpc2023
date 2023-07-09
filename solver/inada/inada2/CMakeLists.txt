cmake_minimum_required (VERSION 3.8)

add_executable (inada2 "main.cpp" "gv.hpp")
set_property(TARGET inada2 PROPERTY CXX_STANDARD 17)
set_property(TARGET inada2 PROPERTY CXX_STANDARD_REQUIRED ON)

