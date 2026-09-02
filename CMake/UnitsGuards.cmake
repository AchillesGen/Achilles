# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Negative compile tests for the strong units layer. Each snippet MUST FAIL to
# compile; if any starts compiling, a unit guard has regressed and we stop the
# build. Runs at configure time -- effectively free, and impossible to skip.
#
# Usage (from test/CMakeLists.txt):
#   include(${CMAKE_SOURCE_DIR}/CMake/UnitsGuards.cmake)
#   achilles_check_units_guards()

include(CheckCXXSourceCompiles)

function(_achilles_expect_fail tag code)
    # Save/restore so we don't leak flags into later checks.
    set(_save_flags "${CMAKE_REQUIRED_FLAGS}")
    set(_save_inc "${CMAKE_REQUIRED_INCLUDES}")
    set(CMAKE_REQUIRED_FLAGS "-std=c++17 -DSPDLOG_FMT_EXTERNAL")
    set(CMAKE_REQUIRED_INCLUDES "${ACHILLES_UNITS_INCLUDE_DIRS}")
    check_cxx_source_compiles("${code}" ACHILLES_GUARD_${tag})
    set(CMAKE_REQUIRED_FLAGS "${_save_flags}")
    set(CMAKE_REQUIRED_INCLUDES "${_save_inc}")
    if(ACHILLES_GUARD_${tag})
        message(FATAL_ERROR
            "Units guard '${tag}' REGRESSED: code that must not compile now compiles.")
    else()
        message(STATUS "Units guard '${tag}': correctly rejected")
    endif()
endfunction()

function(achilles_check_units_guards)
    # Point at the headers under test.
    if(NOT DEFINED ACHILLES_UNITS_INCLUDE_DIRS)
        # FourVector.hh pulls in fmt/spdlog, so those headers are needed too.
        set(ACHILLES_UNITS_INCLUDE_DIRS
            "${CMAKE_SOURCE_DIR}/include"
            "${fmt_SOURCE_DIR}/include"
            "${spdlog_SOURCE_DIR}/include")
    endif()

    # A positive control: if this one fails, the harness itself is broken (bad
    # include path, wrong compiler) and every "correctly rejected" below is a
    # false pass.
    set(_save_flags "${CMAKE_REQUIRED_FLAGS}")
    set(_save_inc "${CMAKE_REQUIRED_INCLUDES}")
    set(CMAKE_REQUIRED_FLAGS "-std=c++17 -DSPDLOG_FMT_EXTERNAL")
    set(CMAKE_REQUIRED_INCLUDES "${ACHILLES_UNITS_INCLUDE_DIRS}")
    check_cxx_source_compiles([=[
        #include "Achilles/Units.hh"
        using namespace achilles::units; using namespace achilles::units::literals;
        int main(){ Energy e = 2.2_GeV; return e.in(GeV) > 0 ? 0 : 1; }
    ]=] ACHILLES_GUARD_positive_control)
    set(CMAKE_REQUIRED_FLAGS "${_save_flags}")
    set(CMAKE_REQUIRED_INCLUDES "${_save_inc}")
    if(NOT ACHILLES_GUARD_positive_control)
        message(FATAL_ERROR
            "Units guard harness is broken: valid units code does not compile. "
            "The negative guards below cannot be trusted.")
    endif()

    # Second positive control, for the guards that include FourVector.hh.
    set(_save_flags "${CMAKE_REQUIRED_FLAGS}")
    set(_save_inc "${CMAKE_REQUIRED_INCLUDES}")
    set(CMAKE_REQUIRED_FLAGS "-std=c++17 -DSPDLOG_FMT_EXTERNAL")
    set(CMAKE_REQUIRED_INCLUDES "${ACHILLES_UNITS_INCLUDE_DIRS}")
    check_cxx_source_compiles([=[
        #include "Achilles/FourVector.hh"
        using namespace achilles; using namespace achilles::units::literals;
        int main(){ FourVector p{1200.0_MeV,300.0_MeV,0.0_MeV,700.0_MeV};
                    FourPosition x{0.0_fm,1.5_fm,0.0_fm,2.0_fm};
                    return p.Boost(p.BoostVector()).E().in(units::MeV) +
                           x.ToArray(units::fm)[3] > 0 ? 0 : 1; }
    ]=] ACHILLES_GUARD_positive_control_vectors)
    set(CMAKE_REQUIRED_FLAGS "${_save_flags}")
    set(CMAKE_REQUIRED_INCLUDES "${_save_inc}")
    if(NOT ACHILLES_GUARD_positive_control_vectors)
        message(FATAL_ERROR
            "Units guard harness is broken: valid FourVector code does not compile. "
            "The vector guards below cannot be trusted.")
    endif()

    _achilles_expect_fail(add_energy_length [=[
        #include "Achilles/Units.hh"
        using namespace achilles::units; using namespace achilles::units::literals;
        int main(){ auto bad = 2.2_GeV + 1.0_fm; (void)bad; }
    ]=])

    _achilles_expect_fail(assign_bare_double [=[
        #include "Achilles/Units.hh"
        using namespace achilles::units;
        int main(){ Energy e = 3.0; (void)e; }
    ]=])

    _achilles_expect_fail(energy_in_fm [=[
        #include "Achilles/Units.hh"
        using namespace achilles::units; using namespace achilles::units::literals;
        int main(){ double x = (2.2_GeV).in(fm); (void)x; }
    ]=])

    _achilles_expect_fail(mix_length_energy [=[
        #include "Achilles/Units.hh"
        using namespace achilles::units; using namespace achilles::units::literals;
        int main(){ Energy e = 1.0_fm; (void)e; }
    ]=])

    _achilles_expect_fail(energy_to_double [=[
        #include "Achilles/Units.hh"
        using namespace achilles::units; using namespace achilles::units::literals;
        int main(){ double x = 2.2_GeV; (void)x; }
    ]=])

    _achilles_expect_fail(boost_by_momentum [=[
        #include "Achilles/FourVector.hh"
        using namespace achilles; using namespace achilles::units::literals;
        int main(){ FourVector p{1200.0_MeV,300.0_MeV,0.0_MeV,700.0_MeV};
                    auto b = p.Boost(p.Vec3()); (void)b; }
    ]=])

    _achilles_expect_fail(add_position_momentum [=[
        #include "Achilles/FourVector.hh"
        using namespace achilles; using namespace achilles::units::literals;
        int main(){ FourVector p{1200.0_MeV,300.0_MeV,0.0_MeV,700.0_MeV};
                    FourPosition x{0.0_fm,1.5_fm,0.0_fm,2.0_fm};
                    auto bad = p + x; (void)bad; }
    ]=])

    _achilles_expect_fail(position_native_as_fm [=[
        #include "Achilles/FourVector.hh"
        using namespace achilles; using namespace achilles::units::literals;
        int main(){ FourPosition x{0.0_fm,1.5_fm,0.0_fm,2.0_fm};
                    double bad = x.ToArray(units::MeV)[3]; (void)bad; }
    ]=])

    _achilles_expect_fail(sqrt_odd_dimension [=[
        #include "Achilles/Units.hh"
        using namespace achilles::units; using namespace achilles::units::literals;
        int main(){ auto bad = sqrt(2.2_GeV); (void)bad; }
    ]=])
endfunction()
