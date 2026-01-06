cmake_minimum_required(VERSION 3.10)
project(ImageTriangulation)

# O código usa recursos modernos como structured bindings e filesystem
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Encontrar OpenCV
find_package(OpenCV REQUIRED)

# Encontrar CGAL
find_package(CGAL REQUIRED)

# Inclusão de diretórios (para achar os .h locais)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Definir o executável e as fontes
add_executable(t2 
    main.cpp 
    image_compression.cpp 
    triangulation.cpp
)

# Linkar as bibliotecas encontradas ao executável
target_link_libraries(t2 
    PRIVATE 
    ${OpenCV_LIBS} 
    CGAL::CGAL
)

# Mensagem de sucesso
message(STATUS "Build configurado para: ${CMAKE_PROJECT_NAME}")