Object identifiers
==================

All objects (instantiations of a class) are assigned unique object identifiers of type ``OBJECT_ID`` (unsigned long int - see 
``constants.h`` for the typedef). The purpose of the unique object id is to aid in object tracking and debugging.

**Note that the object id is not the same as, nor it superseded by, the RANDOM_SEED value assigned to each single or binary star**.

The RANDOM_SEED is used to seed the random number generator, and can be used to uniquely identify a single or binary star. The object
id is more granular the the RANDOM_SEED. For example, each binary star is comprised of multiple objects: the ``BinaryStar`` object, 
which contains two ``BaseBinaryStar`` objects (the object undergoing evolution, and a saved copy); each ``BaseBinaryStar`` object 
contains two ``BinaryConstituentStar`` objects (one for each of the constituent stars), and each ``BinaryConstituentStar`` object 
inherits from the ``Star`` class, which contains two ``BaseStar`` objects (the underlying star and a saved copy). Whereas the RANDOM_SEED
uniquely identifies (for example) a binary star, and so identifies the collection of objects that comprise the binary star, the object ids
uniquely identify the constituent objects of the binary star.  Object tracking at this lower level cannot be achieved using the RANDOM_SEED,
hence the need for object ids when debugging.

As well as unique object ids, all objects are assigned an object type (of type ``OBJECT_TYPE`` – see ``constants.h`` for the enum class
declaring ``OBJECT_TYPE``), and a stellar type where applicable (of type ``STELLAR_TYPE`` – see ``constants.h`` for the enum class declaring
``STELLAR_TYPE``).

Objects should expose the following functions::

    OBJECT_ID ObjectId() const { return m ObjectId; }
    OBJECT_TYPE ObjectType() const { return m ObjectType; }
    STELLAR_TYPE StellarType() const { return m StellarType; }

If any of the functions are not applicable to the object, then they must return "\*::NONE" (all objects should implement ``ObjectId()``
correctly).

Any object that uses the Errors service (i.e. the ``ERRORS`` and ``WARNINGS`` macros) must expose these functions: the functions are 
called by the ``ERRORS`` and ``WARNINGS`` macros (:doc:`./Services/services-error-handling`).