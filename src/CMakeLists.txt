set (SOURCES
    u_derivative.cpp
    dissipative.cpp
    cell.cpp
    init.cpp
    reconst.cpp
    minmod.cpp
    music.cpp
    cornelius.cpp
    hydro_source_base.cpp
    hydro_source_strings.cpp
    hydro_source_ampt.cpp
    pretty_ostream.cpp
    grid_info.cpp
    util.cpp
    read_in_parameters.cpp
    advance.cpp
    eos.cpp
    eos_base.cpp
    eos_idealgas.cpp
    eos_EOSQ.cpp
    eos_s95p.cpp
    eos_WB.cpp
    eos_hotQCD.cpp
    eos_best.cpp
    eos_neos.cpp
    evolve.cpp
    emoji.cpp
    music_logo.cpp
    HydroinfoMUSIC.cpp
    )
if (GSL_FOUND)
    set(SOURCES ${SOURCES}
        freeze.cpp
        freeze_pseudo.cpp
        reso_decay.cpp
        )
endif()

add_library(music_lib SHARED ${SOURCES})
set(CompileFlags "")
if (APPLE)
    set(CompileFlags "${CompileFlags} -DAPPLE")
endif (APPLE)
if (GSL_FOUND)
    set(CompileFlags "${CompileFlags} -DGSL")
    target_link_libraries (music_lib ${GSL_LIBRARIES})
endif()
set_target_properties (music_lib PROPERTIES COMPILE_FLAGS "${CompileFlags}")
install(TARGETS music_lib DESTINATION ${CMAKE_HOME_DIRECTORY})

if (unittest)
    add_executable (unittest_reconst.e reconst_unittest.cpp)
    target_link_libraries (unittest_reconst.e music_lib)
    install(TARGETS unittest_reconst.e DESTINATION ${CMAKE_HOME_DIRECTORY})

    add_executable (unittest_grid.e grid_unittest.cpp)
    target_link_libraries (unittest_grid.e music_lib)
    install(TARGETS unittest_grid.e DESTINATION ${CMAKE_HOME_DIRECTORY})

    add_executable (unittest_minmod.e minmod_unittest.cpp)
    target_link_libraries (unittest_minmod.e music_lib)
    install(TARGETS unittest_minmod.e DESTINATION ${CMAKE_HOME_DIRECTORY})
else (unittest)
    add_executable (mpihydro main.cpp)
    set_target_properties (mpihydro PROPERTIES COMPILE_FLAGS "${CompileFlags}")
    target_link_libraries (mpihydro music_lib)
    install(TARGETS mpihydro DESTINATION ${CMAKE_HOME_DIRECTORY})
endif (unittest)
