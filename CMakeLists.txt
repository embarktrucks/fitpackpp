cmake_minimum_required(VERSION 2.8)

# add modules
#list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

find_package(Eigen3 3.2 REQUIRED)
#find_package(GTest GTEST_LIBRARY REQUIRED)
#find_package( Threads REQUIRED )

include_directories(
    ${EIGEN3_INCLUDE_DIR}
)

# FITPACK FORTRAN code
enable_language(Fortran)

if(${CMAKE_Fortran_COMPILER} MATCHES "ifort.*")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -fltconsistency")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -fltconsistency")
endif()

set(fitpack_SRC
    fitpack/bispev.f
    fitpack/curfit.f
    fitpack/curev.f
    fitpack/cualde.f
    fitpack/fpback.f
    fitpack/fpbisp.f
    fitpack/fpchec.f
    fitpack/fpcurf.f
    fitpack/fpbspl.f
    fitpack/fpdisc.f
    fitpack/fpgivs.f
    fitpack/fpader.f
    fitpack/fpknot.f
    fitpack/fporde.f
    fitpack/fprank.f
    fitpack/fprati.f
    fitpack/fprota.f
    fitpack/fpsurf.f
    fitpack/fppara.f
    fitpack/parder.f
    fitpack/parcur.f
    fitpack/splder.f
    fitpack/splev.f
    fitpack/surfit.f
)

add_library(fitpack STATIC 
            ${fitpack_SRC}
)

enable_language(CXX)

# C++ interface
include(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE ""
                         SYMBOLS curfit splev splder surfit bispev parder parcur curev cualde)
FortranCInterface_VERIFY(CXX)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(fitpackpp_SRC
    fitpackpp/BSplineCurve.cpp
    fitpackpp/BSplineSurface.cpp
    fitpackpp/BSplineParCurve.cpp
)

set(fitpackpp_HEADERS
    fitpackpp/BSplineCurve.h
    fitpackpp/BSplineSurface.h
    fitpackpp/BSplineParCurve.h
)

ADD_DEFINITIONS(-fPIC)
add_library(fitpackpp STATIC
            ${fitpackpp_SRC})

target_link_libraries(fitpackpp fitpack)

install(TARGETS fitpackpp fitpack
        RUNTIME DESTINATION /usr/bin
        LIBRARY DESTINATION /usr/lib
        ARCHIVE DESTINATION /usr/lib
        )

install(FILES ${fitpackpp_HEADERS}
        DESTINATION /usr/include/fitpackpp
        )

#add_executable(fitpack_tests
#        test/fitpack_tests.cpp
#)
#
#target_link_libraries(fitpack_tests
#    fitpackpp
#    gtest
#    ${CMAKE_THREAD_LIBS_INIT}
#)
