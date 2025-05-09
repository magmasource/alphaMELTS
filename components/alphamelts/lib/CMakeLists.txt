
cmake_minimum_required(VERSION 3.26)

project(libalphamelts C)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_LINKER clang++)

# temporary solution
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g3 -DBATCH_VERSION -DRHYOLITE_ADJUSTMENTS -DDO_NOT_PRODUCE_OUTPUT_FILES -DPHMELTS_ADJUSTMENTS")

#set(CMAKE_LIBRARY_PATH
#	/usr/local/lib
#	/usr/lib
#	/usr/local/opt/libxml2/lib
#)

#include_directories(
#	includes
#	/usr/local/include
#	/usr/local/opt/include/libxml2
#)

set(CMAKE_LIBRARY_PATH
	/opt/homebrew/lib
	/usr/lib
	/opt/homebrew/opt/libxml2/lib
)

include_directories(
	includes
	/opt/homebrew/include
	/opt/homebrew/opt/include/libxml2
)

add_library(alphamelts SHARED
	src/albite.c
	src/alloy-liquid.c
	src/alloy-solid.c
	src/clinoamphibole.c
	src/orthoamphibole.c
	src/biotite.c
	src/biotite-taj.c
	src/check_coexisting_liquids.c
	src/check_coexisting_solids.c
	src/clinopyroxene.c
	src/cummingtonite.c
	src/equality_constraints.c
	src/est_saturation_state.c
	src/est_sat-state_revised.c
	src/evaluate_saturation.c
	src/feldspar.c
	src/fluid.c
	src/fluid-phase.c
	src/garnet.c
	src/gibbs.c
	src/gradient_hessian.c
	src/hornblende.c
	src/ilmenite.c
	src/initial_guess.c
	src/kalsilite.c
	src/lawson_hanson.c
	src/leucite.c
	src/linear_search.c
	src/liquid.c
	src/liquid-co2.c
	src/liquid-co2-h2o.c
	src/liquidus.c
	src/liquid_v34.c
	src/majorite.c
	src/melilite.c
	src/melts_linalg.c
	src/melts_support.c
	src/melts_threads.c
	src/mthread.c
	src/nash.c
	src/nepheline.c
	src/olivine.c
	src/olivine-sx.c
	src/ortho-oxide.c
	src/orthopyroxene.c
	src/perovskite.c
	src/phmelts_support.c
	src/plagioclase.c
	src/read_write.c
	src/rhombohedral.c
	src/rhomsghiorso.c
	src/ringwoodite.c
	src/silmin.c
	src/silmin_support.c
	src/spinel.c
	src/sub-solidus-mu-o2.c
	src/sub-thermo-calc.c
	src/wadsleyite.c
	src/water.c
	src/wustite.c
	src/adiabat_1ph.c
	src/adiabat_support.c
	src/interface.c
	src/library.c
	src/subbasic_read_write.c
)

find_library(GSL_LIBRARY gsl)
if(GSL_LIBRARY)
    target_link_libraries(alphamelts PUBLIC ${GSL_LIBRARY})
endif()

find_library(XML2_LIBRARY xml2)
if(XML2_LIBRARY)
    target_link_libraries(alphamelts PUBLIC ${XML2_LIBRARY})
endif()
