cmake_minimum_required(VERSION 3.16)
set(BIN "step_test")
PROJECT (step_test)

include_directories("./")
set(OCCT_INSTALL_DIR "F:/github/OCCT_install")

set(OCCT_DEBUG_LIB_DIR ${OCCT_INSTALL_DIR}/win64/vc14/libd)
set(OCCT_LIB_DIR ${OCCT_INSTALL_DIR}/win64/vc14/lib)

FILE(GLOB EXE_SRC "src/*.cpp" "src/*.c")
set(CMAKE_CXX_STANDARD 20)

add_executable(${BIN} ${EXE_SRC} )
target_include_directories(${BIN} PRIVATE "${OCCT_INSTALL_DIR}/inc")

set(OCCT_LIBS TKernel.lib TKMath.lib TKDESTEP.lib 
	TKXSBase.lib TKShHealing.lib TKBRep.lib TKTopAlgo.lib TKG3d.lib TKG2d.lib 
	TKGeomAlgo.lib TKGeomBase.lib TKMesh.lib)
set(OCCT_DEBUG_LIBS "")
set(OCCT_RELEASE_LIBS "")
foreach(item IN LISTS OCCT_LIBS)
  list(APPEND OCCT_DEBUG_LIBS "${OCCT_DEBUG_LIB_DIR}/${item}")
  list(APPEND OCCT_RELEASE_LIBS "${OCCT_LIB_DIR}/${item}")
endforeach()

target_link_libraries( ${BIN} PRIVATE  $<$<CONFIG:Debug>:${OCCT_DEBUG_LIBS}>)

target_link_libraries( ${BIN} PRIVATE $<$<CONFIG:Release>:${OCCT_RELEASE_LIBS}>)

