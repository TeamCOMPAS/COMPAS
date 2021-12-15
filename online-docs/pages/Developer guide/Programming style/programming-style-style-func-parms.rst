Function parameters
===================

In most cases, function parameters should be input only – meaning that the values of function parameters should not be changed by the function. 
Anything that needs to be changed and returned to the caller should be returned as a functional return. There are a few exceptions to this in COMPAS – 
all were done for performance reasons, and are documented in the code.

To avoid unexpected side-effects, developers should expect (in most cases) that any variables they pass to a function will remain unchanged – all 
changes should be returned as a functional return.
