# Description

This is a mix of different stuff used in my atlas experiments. This whole work is very *experimental* and must be understood as a proof of concept only!

- `_atlas4py.cpp` is the file providing the atlas binding
- `computation.cpp` is the file providing the dawn binding for unstructured grid

These two files can be compiled with the `CMakeLists.txt` (note there are absolute paths inside) or if used in the jupyter notebook, with `cppimport` (there are absolute paths inside the cpp files)

`generated_Diffusion.hpp` is the file generated in the prototypes of dawn (dawn version: a67b66afb9182e0b52d2c4c8497af69e166eed79). For easy of use, it is checked into this repository.

`Atlas4py-Test.ipynb` is a python notebook showing how to use the atlas bindings in the current experimental version. 

# Conclusions

- We should create a proper python package from `_atlas4py.cpp`
- We should either merge the package into atlas, or into a standalone repository in one of the organizations
- We should think about which features we are going to support in this package
- The bindings for dawn are trivial and can be written properly whenever needed.
