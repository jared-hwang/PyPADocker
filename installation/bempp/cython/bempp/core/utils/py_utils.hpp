#ifndef BEMPP_PYTHON_UTILS_HPP
#define BEMPP_PYTHON_UTILS_HPP

#include "bempp/common/common.hpp"
#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>
#include <dune/common/exceptions.hh>
#include <Python.h>
#include <algorithm>
#include <sstream>
#include <complex>
#include "bempp/common/eigen_support.hpp"
#include "bempp/common/shared_ptr.hpp"

namespace Bempp {
    inline static void catch_exception() {
      try {
        // latest Python exn passes through, ignore the current one
        if (not PyErr_Occurred()) throw;
      } catch (const Dune::IOError& exn) {
#if DUNE_VERSION <= 241
        PyErr_SetString(PyExc_IOError, exn.what().c_str());
#else 
        PyErr_SetString(PyExc_IOError, exn.what());
#endif
      } catch (const Dune::Exception& exn) {
#if DUNE_VERSION <= 241
        PyErr_SetString(PyExc_RuntimeError, exn.what().c_str());
#else
        PyErr_SetString(PyExc_RuntimeError, exn.what());
#endif
      } catch (const std::bad_alloc& exn) {
        PyErr_SetString(PyExc_MemoryError, exn.what());
      } catch (const std::bad_cast& exn) {
        PyErr_SetString(PyExc_TypeError, exn.what());
      } catch (const std::domain_error& exn) {
        PyErr_SetString(PyExc_ValueError, exn.what());
      } catch (const std::invalid_argument& exn) {
        PyErr_SetString(PyExc_ValueError, exn.what());
      } catch (const std::ios_base::failure& exn) {
        PyErr_SetString(PyExc_IOError, exn.what());
      } catch (const std::out_of_range& exn) {
        PyErr_SetString(PyExc_IndexError, exn.what());
      } catch (const std::overflow_error& exn) {
        PyErr_SetString(PyExc_OverflowError, exn.what());
      } catch (const std::range_error& exn) {
        PyErr_SetString(PyExc_ArithmeticError, exn.what());
      } catch (const std::underflow_error& exn) {
        PyErr_SetString(PyExc_ArithmeticError, exn.what());
      } catch (const std::exception& exn) {
        PyErr_SetString(PyExc_RuntimeError, exn.what());
      }
      catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
      }
    }

    template <typename T>
    Vector<T> copy_buf_to_vec(T* buf, int n){

        Vector<T> res(n);
        for (int i = 0; i < n; ++i)
            res(i) = buf[i];

        return res;

    }

    template <typename T>
    RowVector<T> copy_buf_to_row_vec(T* buf, int n){

        RowVector<T> res(n);
        for (int i = 0; i < n; ++i)
            res(i) = buf[i];

        return res;

    }

    template <typename T>
    Matrix<T> copy_buf_to_mat(T* buf, int m, int n){

        Matrix<T> res(m,n);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i)
                res(i,j) = buf[j*m+i];

        return res;

    }

    template <typename T>
    shared_ptr<T> py_const_pointer_cast(const shared_ptr<const T>& p){
        return const_pointer_cast<T>(p);

    }

    template <typename T>
    shared_ptr<const T> py_reverse_const_pointer_cast(const shared_ptr<T>& p){
        return p;

    }

}
#endif
