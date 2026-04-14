cmake_minimum_required(VERSION 3.20)

# -------------------------------------------------------------------------------------- Compiler options
set(PROJECT_NAME FWH)
project(${PROJECT_NAME})

# -------------------------------------------------------------------------------------- Compiler options
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wno-literal-suffix -Wno-invalid-offsetof")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -ffast-math -march=native -mtune=native")


# -------------------------------------------------------------------------------------- Setup paths
set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BIN_DIR ${ROOT_DIR}/bin)
set(INC_DIR ${ROOT_DIR}/include)
set(SRC_DIR ${ROOT_DIR}/src)
set(EXEC_DIR ${ROOT_DIR}/tools)

set(CMAKE_INSTALL_PREFIX ${BIN_DIR})
message(STATUS "\nInstall path: ${CMAKE_INSTALL_PREFIX}\n")

# -------------------------------------------------------------------------------------- Source files
file(GLOB SOURCE src/*.cpp)

file(GLOB LIB_HEADERS include/*.h
        include/*.hpp)

include_directories(${INC_DIR})

# -------------------------------------------------------------------------------- deps

find_package(VTK
        COMPONENTS
        CommonCore
        CommonDataModel
        IOCore
        IOLegacy
        IOXML
)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BIN_DIR})

# -------------------------------------------------------------------------------------- FWH
add_library(${PROJECT_NAME} ${SOURCE} ${LIB_HEADERS})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)





# -------------------------------------------------------------------------------------- tools
file(GLOB TOOLS ${EXEC_DIR}/*.cpp)
foreach( sourcefile ${TOOLS} )
    string( REPLACE ".cpp" "" file ${sourcefile} )
    get_filename_component(file ${file} NAME)
    if (${file} MATCHES ${PROJECT_NAME})
        set(file ${PROJECT_NAME}-bin)
    endif()
    add_executable( ${file} ${sourcefile} )
    target_link_libraries( ${file} ${PROJECT_NAME})
endforeach( sourcefile ${TOOLS} )



