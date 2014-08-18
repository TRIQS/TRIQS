Parameters
===============

.. warning::

   Library of beta quality.

   Documentation : just intro/tutorial, reference doc in progress.


Introduction
--------------

The purpose of the class parameters is to handle program input parameters.
It provides a convenient way to pass several parameters of different types to a program.

Parameters can be stored in a parameters object and accessed in a manner which is reminiscent of a Python dict:

.. code-block:: c
   
   parameters P;
   P["x"] = 3.14;
   y = P["x"];

In Python the syntax is similar:

.. code-block:: python

   import parameters as params
   
   P = params.Parameters()
   P['x'] = 3.14
   y = P['x']
   
The class implements the following features:

* storage of parameters of any type, including arrays
* C-like cast operations, such as integral type to double conversion as well as lexical cast to arithmetic types
* Python interface
* type checking
* hdf5 I/O operations
* write contents to std::ostream
* the class is boost-serializable

parameter_defaults
-----------------------------

When parameters are used within a program, a number of them are usually optional, in
the sense that they need not be provided by the user. This either means that the program 
can run without them, or meaningful default values can be provided.

Other parameters, on the other hand, may be required to be provided by the user.
These typically are parameters for which one cannot provide meaningful default values.

These cases can be handled in a consistent manner using the parameter_defaults class. 

An object of this class allows to specify the name, type and documentation for all
required and optional parameters. Default values are provided for optional parameters:

.. code-block:: c
   
   parameter_defaults pdef;
   pdef.required("Beta", double(), "Inverse temperature")
       .optional("N_time", int(100), "Number of time points")
       .optional("N_freq", int(200), "Number of frequencies")
    ; 

The default_parameters object serves two main purposes: Firstly, the input parameters can
be checked for completeness and correctness, i.e. they can be checked if all required parameters
are provided and have correct type.

Secondly, the optional parameters provided through the input parameters can be checked whether they
have correct type. Those which are not provided can be added to the input parameter set with
their respective default values.
Both steps are performed through the update method of the parameter class:

.. code-block:: c
   
   parameters P;
   P["Beta"] = 3.14;
   P.update(pdef);

Given that the parameter_defaults are complete, the resulting parameter object contains all
parameters that the program will eventually access.

Finally, the information stored in the defaults can be used to generate the help of the program,
by printing a list of its parameters:

.. code-block:: c
   
   std::cout<<pdef<<std::endl;

See below for more detailed working code examples.

For detailed information on how to use the parameter and parameter_defaults classes within a program, refer to the **TRIQS solver coding guidelines**.

Headers
--------------
The parameter class and the parameter_defaults class are declared in the header

  triqs/utility/parameters.hpp

Doxygen documentation
-------------------------

The :doxy:`full C++ parameter documentation<triqs::utility::parameters>` and
the :doxy:`parameter_defaults documentation<triqs::utility::parameter_defaults>` are available here.

Parameter example  
-------------

.. triqs_example:: ./parameters_0.cpp
Parameter default value example
-------------

.. triqs_example:: ./parameters_1.cpp
Python example
-------------

.. code-block:: python
   
   import parameters as params

   P = params.Parameters()
   P['num'] = 2.3 
   P['str'] = "hello world"
   print P
   
   dct = { 'a' : 10, 'b' : "a string"}
   P.update(dct) # update from a dictionary
   print P
   
   import numpy as np
   a = np.array([[1,2],[3,4]])
   P['A'] = a # add an array
   print P
   
   P.update2( c = 4, d= "u") # update from an initializer list
   print P

  
