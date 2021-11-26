Abstraction
===========

For any entity, product, or service, the goal of abstraction is to handle the complexity of the implementation by hiding details that 
don’t need to be known in order to use, or consume, the entity, product, or service. In the OOP paradigm, hiding details in this way 
enables the consumer to implement more complex logic on top of the provided abstraction without needing to understand the hidden 
implementation details and complexity. (There is no suggestion that consumers shouldn’t understand the implementation details, but they
shouldn’t need to in order to consume the entity, product, or service).

Abstraction in ``C++`` is achieved via the use of objects – an object is an instance of a class, and typically corresponds to a 
real-world object or entity (in COMPAS, usually a star or binary star). An object maintains the state of an object (via class member 
variables), and provides all necessary means of changing the state of the object (by exposing public class member functions (methods)). 
A class may expose public functions to allow consumers to determine the value of class member variables (“getters”), and to set the value
of class member variables (“setters”).
