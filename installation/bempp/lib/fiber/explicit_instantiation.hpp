// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef fiber_explicit_instantiation_hpp
#define fiber_explicit_instantiation_hpp

#include "../common/common.hpp"

#include "bempp/common/config_data_types.hpp"
#include <complex>

// Invoke an arbitrary macro for all valid basis types
#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_SP_REAL(MACRO) MACRO(float)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_SP_REAL(MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_SP_COMPLEX(MACRO)                \
  MACRO(std::complex<float>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_SP_COMPLEX(MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_DP_REAL(MACRO) MACRO(double)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_DP_REAL(MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_DP_COMPLEX(MACRO)                \
  MACRO(std::complex<double>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_DP_COMPLEX(MACRO)
#endif

#define FIBER_ITERATE_OVER_BASIS_TYPES(MACRO)                                  \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_SP_REAL(MACRO);                        \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_SP_COMPLEX(MACRO);                     \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_DP_REAL(MACRO);                        \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_DP_COMPLEX(MACRO)

// Invoke an arbitrary macro for all valid basis types
#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_SP_REAL(MACRO) MACRO(float)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_SP_REAL(MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_SP_COMPLEX(MACRO)               \
  MACRO(std::complex<float>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_SP_COMPLEX(MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_DP_REAL(MACRO) MACRO(double)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_DP_REAL(MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_DP_COMPLEX(MACRO)               \
  MACRO(std::complex<double>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_DP_COMPLEX(MACRO)
#endif

#define FIBER_ITERATE_OVER_KERNEL_TYPES(MACRO)                                 \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_SP_REAL(MACRO);                       \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_SP_COMPLEX(MACRO);                    \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_DP_REAL(MACRO);                       \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_KERNEL_DP_COMPLEX(MACRO)

// Invoke arbitrary macro for all valid ValueTypes

#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_SP_REAL(MACRO) MACRO(float)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_SP_REAL(MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_SP_COMPLEX(MACRO)                \
  MACRO(std::complex<float>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_SP_COMPLEX(MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_DP_REAL(MACRO) MACRO(double)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_DP_REAL(MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_DP_COMPLEX(MACRO)                \
  MACRO(std::complex<double>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_DP_COMPLEX(CLASSNAME)
#endif

#define FIBER_ITERATE_OVER_VALUE_TYPES(MACRO)                                  \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_SP_REAL(MACRO);                        \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_SP_COMPLEX(MACRO);                     \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_DP_REAL(MACRO);                        \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_VALUE_DP_COMPLEX(MACRO)

// Instantiation of classes templated on basis type

#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_SP_REAL(CLASSNAME)          \
  template class CLASSNAME<float>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_SP_REAL(CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_SP_COMPLEX(CLASSNAME)       \
  template class CLASSNAME<std::complex<float>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_SP_COMPLEX(CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_DP_REAL(CLASSNAME)          \
  template class CLASSNAME<double>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_DP_REAL(CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_DP_COMPLEX(CLASSNAME)       \
  template class CLASSNAME<std::complex<double>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_DP_COMPLEX(CLASSNAME)
#endif

#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(CLASSNAME)                  \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_SP_REAL(CLASSNAME);               \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_SP_COMPLEX(CLASSNAME);            \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_DP_REAL(CLASSNAME);               \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_DP_COMPLEX(CLASSNAME)

// Instantiation of classes templated on kernel type

#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_SP_REAL(CLASSNAME)         \
  template class CLASSNAME<float>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_SP_REAL(CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_SP_COMPLEX(CLASSNAME)      \
  template class CLASSNAME<std::complex<float>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_SP_COMPLEX(CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_DP_REAL(CLASSNAME)         \
  template class CLASSNAME<double>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_DP_REAL(CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_DP_COMPLEX(CLASSNAME)      \
  template class CLASSNAME<std::complex<double>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_DP_COMPLEX(CLASSNAME)
#endif

#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL(CLASSNAME)                 \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_SP_REAL(CLASSNAME);              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_SP_COMPLEX(CLASSNAME);           \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_DP_REAL(CLASSNAME);              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL_DP_COMPLEX(CLASSNAME)

// Instantiation of classes templated on result type

#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_SP_REAL(CLASSNAME)         \
  template class CLASSNAME<float>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_SP_REAL(CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_SP_COMPLEX(CLASSNAME)      \
  template class CLASSNAME<std::complex<float>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_SP_COMPLEX(CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_DP_REAL(CLASSNAME)         \
  template class CLASSNAME<double>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_DP_REAL(CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_DP_COMPLEX(CLASSNAME)      \
  template class CLASSNAME<std::complex<double>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_DP_COMPLEX(CLASSNAME)
#endif

#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(CLASSNAME)                 \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_SP_REAL(CLASSNAME);              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_SP_COMPLEX(CLASSNAME);           \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_DP_REAL(CLASSNAME);              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_DP_COMPLEX(CLASSNAME)

#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_REAL_ONLY(CLASSNAME)       \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_SP_REAL(CLASSNAME);              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_DP_REAL(CLASSNAME)

// Invoke an arbitrary macro for all valid (basis, result) type combinations
#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_REAL_REAL(MACRO)   \
  MACRO(float, float)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_REAL_REAL(MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_REAL_COMPLEX(      \
    MACRO)                                                                     \
  MACRO(float, std::complex<float>);
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_REAL_COMPLEX(MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_COMPLEX_COMPLEX(   \
    MACRO)                                                                     \
  MACRO(std::complex<float>, std::complex<float>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_COMPLEX_COMPLEX(   \
    MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_REAL_REAL(MACRO)   \
  MACRO(double, double)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_REAL_REAL(MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_REAL_COMPLEX(      \
    MACRO)                                                                     \
  MACRO(double, std::complex<double>);
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_REAL_COMPLEX(MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_COMPLEX_COMPLEX(   \
    MACRO)                                                                     \
  MACRO(std::complex<double>, std::complex<double>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_COMPLEX_COMPLEX(   \
    MACRO)
#endif

#define FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(MACRO)                       \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_REAL_REAL(MACRO);        \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_REAL_COMPLEX(MACRO);     \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_SP_COMPLEX_COMPLEX(MACRO);  \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_REAL_REAL(MACRO);        \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_REAL_COMPLEX(MACRO);     \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_AND_RESULT_DP_COMPLEX_COMPLEX(MACRO);

// Instantiation of classes templated on basis type and result type

#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_REAL(    \
    CLASSNAME)                                                                 \
  template class CLASSNAME<float, float>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_REAL(    \
    CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_COMPLEX( \
    CLASSNAME)                                                                 \
  template class CLASSNAME<float, std::complex<float>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_COMPLEX_COMPLEX( \
    CLASSNAME)                                                                    \
  template class CLASSNAME<std::complex<float>, std::complex<float>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_COMPLEX_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_REAL(    \
    CLASSNAME)                                                                 \
  template class CLASSNAME<double, double>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_REAL(    \
    CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_COMPLEX( \
    CLASSNAME)                                                                 \
  template class CLASSNAME<double, std::complex<double>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_COMPLEX_COMPLEX( \
    CLASSNAME)                                                                    \
  template class CLASSNAME<std::complex<double>, std::complex<double>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_COMPLEX_COMPLEX( \
    CLASSNAME)
#endif

#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASSNAME)       \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_REAL(          \
      CLASSNAME);                                                              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_COMPLEX(       \
      CLASSNAME);                                                              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_COMPLEX_COMPLEX(    \
      CLASSNAME);                                                              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_REAL(          \
      CLASSNAME);                                                              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_COMPLEX(       \
      CLASSNAME);                                                              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_COMPLEX_COMPLEX(    \
      CLASSNAME)

#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_REAL_ONLY(       \
    CLASSNAME)                                                                 \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_REAL(          \
      CLASSNAME);                                                              \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_REAL(CLASSNAME)

// Invoke an arbitrary macro for all valid (basis, kernel, result) type
// combinations
#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_REAL( \
    MACRO)                                                                         \
  MACRO(float, float, float)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_REAL( \
    MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_COMPLEX( \
    MACRO)                                                                            \
  MACRO(float, float, std::complex<float>);
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_COMPLEX( \
    MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_COMPLEX_COMPLEX( \
    MACRO)                                                                               \
  MACRO(float, std::complex<float>, std::complex<float>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_COMPLEX_COMPLEX( \
    MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_REAL_COMPLEX( \
    MACRO)                                                                               \
  MACRO(std::complex<float>, float, std::complex<float>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_REAL_COMPLEX( \
    MACRO)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS) &&     \
    defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_COMPLEX_COMPLEX( \
    MACRO)                                                                                  \
  MACRO(std::complex<float>, std::complex<float>, std::complex<float>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_COMPLEX_COMPLEX( \
    MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_REAL( \
    MACRO)                                                                         \
  MACRO(double, double, double)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_REAL( \
    MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_COMPLEX( \
    MACRO)                                                                            \
  MACRO(double, double, std::complex<double>);
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_COMPLEX( \
    MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_COMPLEX_COMPLEX( \
    MACRO)                                                                               \
  MACRO(double, std::complex<double>, std::complex<double>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_COMPLEX_COMPLEX( \
    MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_REAL_COMPLEX( \
    MACRO)                                                                               \
  MACRO(std::complex<double>, double, std::complex<double>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_REAL_COMPLEX( \
    MACRO)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS) &&     \
    defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_COMPLEX_COMPLEX( \
    MACRO)                                                                                  \
  MACRO(std::complex<double>, std::complex<double>, std::complex<double>)
#else
#define FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_COMPLEX_COMPLEX( \
    MACRO)
#endif

#define FIBER_ITERATE_OVER_BASIS_KERNEL_AND_RESULT_TYPES(MACRO)                       \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_REAL(          \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_COMPLEX(       \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_COMPLEX_COMPLEX(    \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_REAL_COMPLEX(    \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_COMPLEX_COMPLEX( \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_REAL(          \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_COMPLEX(       \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_COMPLEX_COMPLEX(    \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_REAL_COMPLEX(    \
      MACRO);                                                                         \
  FIBER_INVOKE_MACRO_DEPENDENT_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_COMPLEX_COMPLEX( \
      MACRO)

// Instantiation of classes templated on basis type, kernel type and result type

#if defined(ENABLE_SINGLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_REAL( \
    CLASSNAME)                                                                          \
  template class CLASSNAME<float, float, float>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_REAL( \
    CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_COMPLEX( \
    CLASSNAME)                                                                             \
  template class CLASSNAME<float, float, std::complex<float>>;
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_COMPLEX_COMPLEX( \
    CLASSNAME)                                                                                \
  template class CLASSNAME<float, std::complex<float>, std::complex<float>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_COMPLEX_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_REAL_COMPLEX( \
    CLASSNAME)                                                                                \
  template class CLASSNAME<std::complex<float>, float, std::complex<float>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_REAL_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS) &&     \
    defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_COMPLEX_COMPLEX( \
    CLASSNAME)                                                                                   \
  template class CLASSNAME<std::complex<float>, std::complex<float>,                             \
                           std::complex<float>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_COMPLEX_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_REAL( \
    CLASSNAME)                                                                          \
  template class CLASSNAME<double, double, double>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_REAL( \
    CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_COMPLEX( \
    CLASSNAME)                                                                             \
  template class CLASSNAME<double, double, std::complex<double>>;
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_COMPLEX_COMPLEX( \
    CLASSNAME)                                                                                \
  template class CLASSNAME<double, std::complex<double>, std::complex<double>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_COMPLEX_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_REAL_COMPLEX( \
    CLASSNAME)                                                                                \
  template class CLASSNAME<std::complex<double>, double, std::complex<double>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_REAL_COMPLEX( \
    CLASSNAME)
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS) &&     \
    defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_COMPLEX_COMPLEX( \
    CLASSNAME)                                                                                   \
  template class CLASSNAME<std::complex<double>, std::complex<double>,                           \
                           std::complex<double>>
#else
#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_COMPLEX_COMPLEX( \
    CLASSNAME)
#endif

#define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(                      \
    CLASSNAME)                                                                             \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_REAL(          \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_REAL_COMPLEX(       \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_REAL_COMPLEX_COMPLEX(    \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_REAL_COMPLEX(    \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_SP_COMPLEX_COMPLEX_COMPLEX( \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_REAL(          \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_REAL_COMPLEX(       \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_REAL_COMPLEX_COMPLEX(    \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_REAL_COMPLEX(    \
      CLASSNAME);                                                                          \
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT_DP_COMPLEX_COMPLEX_COMPLEX( \
      CLASSNAME)

#endif // fiber_explicit_instantiation_hpp
