Encapsulation
=============

Encapsulation binds together the data and functions that manipulate the data in an attempt to keep both safe from outside interference and
accidental misuse. An encapsulation paradigm that does not allow calling code to access internal object data and permits access through 
functions only is a strong form of abstraction. ``C++`` allows developers to enforce access restrictions explicitly by defining class member
variables and functions as private, protected, or public. These keywords are used throughout COMPAS to enforce encapsulation.

There are very few circumstances in which a consumer should change the value of a class member variable directly (via the use of a setter 
function) – almost always consumers should present new situational information to an object (via a public member function), and allow the 
object to respond to the new information. For example, in COMPAS, there should be almost no reason for a consumer of a star object to 
directly change (say) the radius of the star – the consumer should inform the star object of new circumstances or events, and allow the star 
object to respond to those events (perhaps changing the value of the radius of the star). Changing a single class member variable directly 
introduces the possibility that related class member variables (e.g. other attributes of stars) will not be changed accordingly. Moreover, 
developers changing the code in the future should, in almost all cases, expect that the state of an object is maintained consistently by the
object, and that there should be no unexpected side-effects caused by calling non class-member functions.

In short, changing the state of an object outside the object is potentially unsafe and should be avoided where possible.
