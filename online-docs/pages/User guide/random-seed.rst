Random seed
===========

The ``--random-seed`` option allows users to specify the initial value to be used to seed the pseudo-random number generator. Once set, the
random seed values increments from its initial value for each star, or binary star, evolved. How the random seed increments depends upon the
context.

The ``--random-seed`` option can be specified on either, or both, the command line and a :doc:`grid file <./grid-files>` line. If the option 
is not specified on one or the other, the default value is used (see :doc:`./Program options/program-options-list-defaults`).

In general, if the ``--random-seed`` option is specified, the pseudo-random number generator will be seeded using the specified value for 
the first star, or binary star, evolved, then for each subsequent star or binary star, the seed value will be incremented by one and the 
pseudo-random number generator re-seeded. Seeding the pseudo-random number generator with a known seed for each star, or binary star, 
evolved ensures that the evolution of specific stars, or binary stars, can be reproduced.

Consider a single execution of COMPAS effected with the command::

    ./COMPAS --random-seed 15 --number-of-systems 100 --metallicity 0.015

This would evolve 100 binary stars, each with `metallicity = 0.015`, and other initial attributes set to their defaults. The first of the 
100 binary stars will be evolved using the random seed 15, the second 16, the third 17, and so on - each binary star will evolve using
a unique random seed.

In the example shown above (see Section :doc:`./Program options/program-options-mixing-ranges-sets`), all 104 binary stars would evolve with
unique random seed values, ranging from 0 (the default, since the option was not specified on either the command line or in the grid file), to 103.

In both these examples, the random seed was incremented in the context of the command line. In the first example, the random seed was 
explicitly specified on the command line, and in the second example the random seed defaulted to the command line default.

Consider now a single execution of COMPAS, using the grid file ``mygrid.txt``::

    ./COMPAS --random-seed 12 --grid mygrid.txt

where the contents of the grid file ``mygrid.txt`` are::

    --allow-rlof-at-birth true --metallicity 0.1
    --semi-major-axis 23.4 --random-seed 107
    --random-seed 63 --metallicity 0.12 --eccentricity s[0.1,0.2,0.3,0.4]
    --initial-mass-1 12.3

This would evolve 7 binary stars with random seed values 12, 107, 63, 64, 65, 66, and 18.

The first binary star evolved is the first line of the grid file. This line does not specify the ``--random-seed`` option, so the random seed
defaults to the command line value. The command line did specify a value for the random seed (12), so that value is used. Since the first
line of the grid file is the first binary star evolved, the random seed is not incremented, and the value of 12 is used.

The second binary star evolved is the second line of the grid file. This line does specify the ``--random-seed`` option. Since this is the 
first binary star evolved in the context of the random seed specified on the grid file line, the random seed is not incremented, and the value 
of 107 is used.

The third binary star evolved is the third line of the grid file. This line does specify the ``--random-seed`` option. Since this is the first
binary star evolved in the context of the random seed specified on the grid file line, the random seed is not incremented, and the value of 63
is used.

The fourth, fifth, and sixth binary stars evolved are also from the third line of the grid file - a set of four values for eccentricity was 
specified. Since these are subsequent to the first binary star evolved in the context of the random seed specified on the grid file line, the
random seed is incremented, and the values of 64, 65, and 66 are used.

The seventh binary star evolved is the fourth line of the grid file. This line does not specify the ``--random-seed`` option, so the random 
seed defaults to the command line value. The command line did specify a value for the random seed (12), so that value is used, but since this 
binary star is subsequent to the first binary star evolved in the context of the random seed specified on the command line, the random seed is 
incremented. This is the sixth subsequent binary star evolved in the context of the command line (all stars, or binary stars, evolved in a 
single execution of COMPAS are evolved in the context of the command line), so the random seed is incremented from 12 to 18 (by 1 for each 
binary star evolved), and the value used for this binary star is 18.

Note that in this example, all binary stars were evolved using a unique random seed. This is because the values specified for the random seeds 
via the ``--random-seed`` option were ’well-behaved’. Unfortunately there is no reasonable way to protect the user against specifying duplicate 
random seeds – especially since the random seed increments for each star or binary star. If the user chooses to specify multiple grid file lines 
with the same random seed, or initial random seeds that would collide with other random seed values and cause duplicates as they increment 
through ranges and sets, then there will be duplicate random seeds in the output files. Users should take care when specifying random seeds in 
grid files via the ``--random-seed`` option.
