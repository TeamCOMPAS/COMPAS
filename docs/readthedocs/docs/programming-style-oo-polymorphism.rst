Polymorphism
============

Polymorphism means having many forms. In OOP, polymorphism occurs when there is a hierarchy of classes and they are related by inheritance.

Following the discussion earlier regarding inheritance, in the OOP paradigm, and ``C++`` specifically, derived classes can override methods 
defined by ancestor classes, allowing a derived class to implement functions specific to its circumstances. This means that a call to a class 
member function will cause a different function to be executed depending on the type of object that invokes the function. Descendent classes 
of a class that has overridden a base class member function inherit the overridden function (but can override it themselves).

COMPAS makes heavy use of inheritance and polymorphism, especially for the implementation of the different stellar types.