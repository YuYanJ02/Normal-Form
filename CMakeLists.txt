cmake_minimum_required(VERSION 3.20)
project(crtbp_normal_form_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(nf
  src/main.cpp
  src/normal_form.cpp
)

target_include_directories(nf PRIVATE src)

if (MSVC)
  target_compile_options(nf PRIVATE /O2 /permissive-)
else()
  target_compile_options(nf PRIVATE -O3)
endif()

