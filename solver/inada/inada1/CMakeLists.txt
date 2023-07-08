cmake_minimum_required (VERSION 3.8)

add_executable (inada1 "inada1.cpp" "inada1.h" "gv.hpp")
set_property(TARGET inada1 PROPERTY CXX_STANDARD 17)
set_property(TARGET inada1 PROPERTY CXX_STANDARD_REQUIRED ON)

