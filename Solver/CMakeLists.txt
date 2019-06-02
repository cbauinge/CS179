include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/Domain/include)

set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/Solver.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/SolverEigen.cpp)

if (CUDA_FOUND)
    list(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/SolverGPU.cpp)
    message("Build uses CUDA.")
    message(${CUDA_INCLUDE_DIRS})    
    message(${CUDA_TOOLKIT_INCLUDE})
    message(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    cuda_add_library(CudaSolver ${CMAKE_CURRENT_SOURCE_DIR}/src/SolverGPU.cu)
endif()

add_library(Solver ${SOURCES})
target_link_libraries(Solver Domain CudaSolver 
    ${CUDA_cublas_LIBRARY}
    ${CUDA_cusolver_LIBRARY})

add_subdirectory(unit_test)