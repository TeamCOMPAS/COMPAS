Braces
======

The placement of braces in ``C++`` code (actually, any code that uses braces to enclose scope) is a contentious issue, with many developers having 
long-held, often dogmatic preferences. COMPAS (so far) uses the K&R style (”the one true brace style”) - the style used in the original Unix kernel
and Kernighan and Ritchie’s book :doc:`The C Programming Language <../../references>`.

The K&R style puts the opening brace on the same line as the control statement:

::

    while (x == y) {
        call_something();
        var1 = var2
        call_somethingelse();
    }

Note also the space between the keyword while and the opening parenthesis, surrounding the ``==`` operator, and between the closing parenthesis 
and the opening brace. Spaces in those places help with code readability. Surrounding all arithmetic operators with spaces is preferred.
