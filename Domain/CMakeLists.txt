include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Domain 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Domain.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Vec2D.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/BoundaryDataGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Boundary.cpp)

add_subdirectory(unit_test)