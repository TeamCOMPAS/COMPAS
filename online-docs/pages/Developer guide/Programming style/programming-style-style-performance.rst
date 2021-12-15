Performance & optimisation
==========================

In general, COMPAS developers should code for performance – within reason. Bear in mind that many functions will be called many, many thousands of 
times (in some cases, millions) in one execution of the program.

- Avoid calculating values inside loops that could be calculated once outside the loop.
- Try to use constants where possible.
- Use multiplication in preference to functions such as pow() and sqrt() (note that pow() is very expensive computationally; sqrt() is expensive, but much less expensive than pow()).
- Don’t optimise to the point that readability and maintainability is compromised. Bear in mind that most compilers are good at optimising, and are very forgiving of less-than-optimally-written code (though they are not miracle workers...).
