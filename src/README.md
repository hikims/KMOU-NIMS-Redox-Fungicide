# Shared solver files

This directory contains the shared C++ solver files used by all reproducibility drivers in `drivers/`.

Files:

- `functions.hpp`
- `functions.cpp`
- `initial.hpp`
- `initial.cpp`
- `parameters.h`

The driver files include these headers as:

```cpp
#include "functions.hpp"
#include "initial.hpp"
#include "parameters.h"