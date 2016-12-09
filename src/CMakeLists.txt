add_executable (mpihydro
    main.cpp
    grid.cpp
    eos.cpp
    evolve.cpp
    init.cpp
    reconst.cpp
    freeze.cpp
    freeze_pseudo.cpp
    minmod.cpp
    glauber.cpp
    advance.cpp
    u_derivative.cpp
    dissipative.cpp
    random.cpp
    util.cpp
    grid_info.cpp
    read_in_parameters.cpp
    music.cpp
    )

target_link_libraries (mpihydro ${GSL_LIBRARIES})

install(TARGETS mpihydro DESTINATION ${CMAKE_HOME_DIRECTORY})
