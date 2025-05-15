# mingw-toolchain.cmake

set(CMAKE_SYSTEM_NAME Windows)

# Specify the compiler
set(CMAKE_C_COMPILER /opt/homebrew/opt/mingw-w64/bin/x86_64-w64-mingw32-gcc)
set(CMAKE_CXX_COMPILER /opt/homebrew/opt/mingw-w64/bin/x86_64-w64-mingw32-g++)

# Avoid Apple flags
set(CMAKE_OSX_SYSROOT "")
set(CMAKE_OSX_ARCHITECTURES "")

# Optional: remove default flags that might break
set(CMAKE_C_FLAGS "")
set(CMAKE_CXX_FLAGS "")