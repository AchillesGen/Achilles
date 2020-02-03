add_library(utilities SHARED
    FourVector.cc
    Interpolation.cc
    Particle.cc
    ThreeVector.cc
    Utilities.cc
)
target_include_directories(utilities PRIVATE ${NUCHIC_INCLUDES} PUBLIC ${Python_INCLUDE_DIRS})

set(NUCHIC_LIBS ${NUCHIC_LIBS} utilities ${HDF5_LIBRARIES})

add_library(physics SHARED
    Cascade.cc
    Interactions.cc
    InteractionsFactory.cc
    Nucleus.cc
)
target_include_directories(physics PRIVATE ${NUCHIC_INCLUDES})
target_link_libraries(physics PRIVATE ${NUCHIC_LIBS})

set(NUCHIC_LIBS ${NUCHIC_LIBS} physics)

pybind11_add_module(vectors MODULE
    VectorModule.cc
)
pylink(vectors "${NUCHIC_INCLUDES}" "${NUCHIC_LIBS}")

pybind11_add_module(interpolation MODULE
    InterpolationModule.cc
)
pylink(interpolation "${NUCHIC_INCLUDES}" "${NUCHIC_LIBS}")

pybind11_add_module(particle MODULE
    ParticleModule.cc
)
pylink(particle "${NUCHIC_INCLUDES}" "${NUCHIC_LIBS}")

pybind11_add_module(nucleus MODULE
    NucleusModule.cc
)
pylink(nucleus "${NUCHIC_INCLUDES}" "${NUCHIC_LIBS}")

pybind11_add_module(interactions MODULE
    InteractionsModule.cc
)
pylink(interactions "${NUCHIC_INCLUDES}" "${NUCHIC_LIBS}")

pybind11_add_module(cascade MODULE
    CascadeModule.cc
)
pylink(cascade "${NUCHIC_INCLUDES}" "${NUCHIC_LIBS}")
