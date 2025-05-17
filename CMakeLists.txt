cmake_minimum_required(VERSION 3.28.1)
project(BEM_Turbine)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_definitions(XFOIL_PROJECT_ROOT="${CMAKE_SOURCE_DIR}")

# -------------------- LIBS -------------------- #
include(FetchContent)
set(CMAKE_MINIMUM_REQUIRED_VERSION 3.5)
# ---------------------- GSL-LITE -------------------- #
find_package(gsl-lite QUIET)
if (NOT gsl-lite_FOUND)
    message(STATUS "gsl-lite not found – fetching from GitHub")
    FetchContent_Declare(
            gsl-lite
            GIT_REPOSITORY https://github.com/gsl-lite/gsl-lite.git
            GIT_TAG        v0.41.0
    )
    FetchContent_MakeAvailable(gsl-lite)
endif()

# ---------------------- nholthaus/units ---------------------- #
find_package(units CONFIG QUIET)

if (NOT units_FOUND)
    message(STATUS "Holthaus units not found – fetching from GitHub")
    include(FetchContent)

    FetchContent_Declare(
            units_src
            GIT_REPOSITORY https://github.com/nholthaus/units.git
            GIT_TAG        v2.3.1
    )

    FetchContent_MakeAvailable(units_src)

    if (TARGET units AND NOT TARGET units::units)
        add_library(units::units ALIAS units)
    endif()
endif()

# ---------------------- MATPLOTPP -------------------- #
find_package(Matplot++ CONFIG QUIET)

if (NOT Matplot++_FOUND)
    message(STATUS "Matplot++ not found – fetching from GitHub")
    include(FetchContent)

    FetchContent_Declare(
            matplotplusplus
            GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus.git
            GIT_TAG        v1.2.0
    )
    set(MATPLOTPP_BUILD_TESTS      OFF CACHE INTERNAL "")
    set(MATPLOTPP_BUILD_EXAMPLES   OFF CACHE INTERNAL "")
    set(MATPLOTPP_BUILD_BENCHMARK  OFF CACHE INTERNAL "")
    FetchContent_MakeAvailable(matplotplusplus)

    # The source build creates a target called `matplot`; add an ALIAS so the
    # imported-name is always available.
    if (TARGET matplot AND NOT TARGET Matplot++::matplot)
        add_library(Matplot++::matplot ALIAS matplot)
    endif()
endif()

# --------------------  EVE (Expressive Vector Engine) -------------------- #
find_package(eve CONFIG QUIET)

if (NOT eve_FOUND)
    message(STATUS "EVE not found – fetching from GitHub")
    include(FetchContent)

    FetchContent_Declare(
            eve
            GIT_REPOSITORY https://github.com/jfalcou/eve.git
            GIT_TAG        v2023.02.15
    )
    set(EVE_BUILD_TEST OFF          CACHE INTERNAL "")
    set(EVE_BUILD_BENCHMARKS OFF    CACHE INTERNAL "")
    set(EVE_BUILD_DOCUMENTATION OFF CACHE INTERNAL "")
    FetchContent_MakeAvailable(eve)

    if (TARGET eve AND NOT TARGET eve::eve)
        add_library(eve::eve ALIAS eve)
    endif()
endif()


# ---------------------- NLOHMANN JSON ---------------------- #
find_package(nlohmann_json CONFIG QUIET)

if (NOT nlohmann_json_FOUND)
    message(STATUS "nlohmann_json not found – fetching from GitHub")
    include(FetchContent)

    FetchContent_Declare(
            nlohmann_json
            GIT_REPOSITORY https://github.com/nlohmann/json.git
            GIT_TAG        v3.11.2
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif()



# -------------------- ADD -------------------- #
add_executable(BEM_Turbine src/main.cpp
        include/headers.h
        include/rotor.h
        src/blade_section.cpp
        include/blade_section.h
        src/rotor.cpp
        include/xfoil_hook.h
        src/xfoil_hook.cpp
        include/helpers.h
        include/bem_solver.h
        src/bem_solver.cpp
        include/bem_solver.h
        include/blade_section.h
        include/headers.h
        include/helpers.h
        include/rotor.h
        include/xfoil_hook.h
        src/bem_solver.cpp
        src/blade_section.cpp
        src/main.cpp
        src/rotor.cpp
        src/xfoil_hook.cpp
        include/viterna_extrapolator.h
        include/fusion_hook.h
        include/bem_corrections.h
)


# -------------------- LINK -------------------- #
target_link_libraries(BEM_Turbine PRIVATE Matplot++::matplot)
target_link_libraries(BEM_Turbine PRIVATE eve::eve)
target_link_libraries(BEM_Turbine PRIVATE units::units)
target_link_libraries(BEM_Turbine PRIVATE nlohmann_json::nlohmann_json)

# -------------------- COMPILE-OPTS -------------------- #




