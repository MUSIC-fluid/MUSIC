add_executable (mpihydro
    u_derivative.cpp
    main.cpp
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
    evolve.cpp
    dissipative.cpp
    emoji.cpp
    )

target_link_libraries (mpihydro ${GSL_LIBRARIES})
if (APPLE)
    set_target_properties (mpihydro PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)

install(TARGETS mpihydro DESTINATION ${CMAKE_HOME_DIRECTORY})
