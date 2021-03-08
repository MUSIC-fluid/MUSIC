set (SOURCES
    u_derivative.cpp
    dissipative.cpp
    init.cpp
    reconst.cpp
    minmod.cpp
    music.cpp
    cornelius.cpp
    hydro_source_base.cpp
    hydro_source_strings.cpp
    hydro_source_ampt.cpp
    hydro_source_smash.cpp
    hydro_source_TATB.cpp
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
    eos_UH.cpp
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

add_library(${libname} SHARED ${SOURCES})
set(CompileFlags "")
if (APPLE)
    set(CompileFlags "${CompileFlags} -DAPPLE")
endif (APPLE)
if (GSL_FOUND)
    set(CompileFlags "${CompileFlags} -DGSL")
    target_link_libraries (${libname} ${GSL_LIBRARIES})
endif()
set_target_properties (${libname} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
install(TARGETS ${libname} DESTINATION ${CMAKE_HOME_DIRECTORY})

if (unittest)
    add_executable (unittest_reconst.e reconst_unittest.cpp)
    target_link_libraries (unittest_reconst.e ${libname})
    install(TARGETS unittest_reconst.e DESTINATION ${CMAKE_HOME_DIRECTORY})

    add_executable (unittest_grid.e grid_unittest.cpp)
    target_link_libraries (unittest_grid.e ${libname})
    install(TARGETS unittest_grid.e DESTINATION ${CMAKE_HOME_DIRECTORY})

    add_executable (unittest_cell.e cell_unittest.cpp)
    target_link_libraries (unittest_cell.e ${libname})
    install(TARGETS unittest_cell.e DESTINATION ${CMAKE_HOME_DIRECTORY})

    add_executable (unittest_minmod.e minmod_unittest.cpp)
    target_link_libraries (unittest_minmod.e ${libname})
    install(TARGETS unittest_minmod.e DESTINATION ${CMAKE_HOME_DIRECTORY})
else (unittest)
    if (link_with_lib)
        add_executable (${exename} main.cpp)
        target_link_libraries (${exename} ${libname})
    else (link_with_lib)
        add_executable (${exename} main.cpp ${SOURCES})
        if (GSL_FOUND)
            target_link_libraries (${exename} ${GSL_LIBRARIES})
        endif (GSL_FOUND)
    endif (link_with_lib)
    set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
    install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
endif (unittest)
