Constant declarations (constants.h)
===================================

``constants.h`` is the COMPAS ``C++`` constants header file.

``constants.h`` is where constants used throughout the COMPAS code base are declared.  The constants declared here are basic data
types (such as ``double``, ``int``, ``bool``, etc.) as well as more complex data types (such as vectors, ``COMPASUnorderedMaps``, etc.).
Any new constants being added to COMPAS should be added here.

Refer to ``EnumHash.h`` for the definition of the type alias ``COMPASUnorderedMap``.

``constants.h`` also contains the global declaration for the global object identifier ``globalObjectId``, and a few ``typedef`` declarations
(at the head of the file) that are commonly used type definitions that we make glbally available by declaring them in ``constants.h`` (which
is included in all other source files).
