Inheritance
===========

Inheritance allows classes to be arranged in a hierarchy that represents is-a-type-of relationships. All nonprivate class member variables and
functions of the parent (base) class are available to the child (derived) class (and, therefore, child classes of the child class). This allows
easy re-use of the same procedures and data definitions, in addition to describing real-world relationships in an intuitive way. ``C++`` allows
multiple inheritance – a class may inherit from multiple parent classes.

Derived classes can define additional class member variables (using the private, protected, and public access restrictions), which will be 
available to any descendent classes (subject to inheritance rules), but will only be available to ancestor classes via the normal access methods
(getters and setters).
