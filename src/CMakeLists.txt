add_executable (mpihydro
    main.cpp
    cell.cpp
    init.cpp
    grid_info.cpp
    hydro_source.cpp
    minmod.cpp
    reconst.cpp
    evolve.cpp
    eos.cpp
    advance.cpp
    u_derivative.cpp
    dissipative.cpp
    util.cpp
    read_in_parameters.cpp
    music.cpp
    cornelius.cpp
    pretty_ostream.cpp
    freeze.cpp
    freeze_pseudo.cpp
    reso_decay.cpp
    )

target_link_libraries (mpihydro ${GSL_LIBRARIES})
if (APPLE)
    set_target_properties (mpihydro PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)

install(TARGETS mpihydro DESTINATION ${CMAKE_HOME_DIRECTORY})
