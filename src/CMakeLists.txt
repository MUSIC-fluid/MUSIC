set (SOURCES
    u_derivative.cpp
    dissipative.cpp
    cell.cpp
    init.cpp
    reconst.cpp
    minmod.cpp
    music.cpp
    cornelius.cpp
    hydro_source.cpp
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
    )

if (unittest)
    add_executable (unittest_grid.e grid.cpp)
    install(TARGETS unittest_grid.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_minmod.e minmod.cpp)
    install(TARGETS unittest_minmod.e DESTINATION ${CMAKE_HOME_DIRECTORY})
else (unittest)
    if (build_lib)
        add_library(music_lib SHARED ${SOURCES})
    endif()
    add_executable (mpihydro main.cpp ${SOURCES})
    if (APPLE)
        set_target_properties (mpihydro PROPERTIES COMPILE_FLAGS "-DAPPLE")
        if (build_lib)
            set_target_properties (music_lib PROPERTIES COMPILE_FLAGS "-DAPPLE")
        endif()
    endif (APPLE)
    target_link_libraries (mpihydro ${GSL_LIBRARIES})
    if (build_lib)
        target_link_libraries (music_lib ${GSL_LIBRARIES})
    endif()
    install(TARGETS mpihydro DESTINATION ${CMAKE_HOME_DIRECTORY})
endif (unittest)
