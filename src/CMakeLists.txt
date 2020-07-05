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
    freeze.cpp
    grid_info.cpp
    grid.cpp
    util.cpp
    read_in_parameters.cpp
    freeze_pseudo.cpp
    reso_decay.cpp
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
    transport_coeffs.cpp
    )

add_library(${libname} SHARED ${SOURCES})
set(CompileFlags "")
if (APPLE)
    set(CompileFlags "${CompileFlags} -DAPPLE")
endif (APPLE)
set_target_properties (${libname} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
target_link_libraries (${libname} ${GSL_LIBRARIES})
install(TARGETS ${libname} DESTINATION ${CMAKE_HOME_DIRECTORY})

if (unittest)
    add_executable (unittest_grid.e grid.cpp)
    install(TARGETS unittest_grid.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_minmod.e minmod.cpp)
    install(TARGETS unittest_minmod.e DESTINATION ${CMAKE_HOME_DIRECTORY})
else (unittest)
    add_executable (${exename} main.cpp)
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
    target_link_libraries (${exename} ${libname})
    install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
endif (unittest)
