Constants source file
=====================

``constants.h`` is the COMPAS ``C++`` constants source file.

As well as plain constant values, many distribution and prescription identifiers are declared in ``constants.h``. These are mostly 
declared as enum classes, with each enum class having a corresponding map of labels. The benefit is that the values of a particular
(e.g.) prescription are limited to the values declared in the enum class, rather than any integer value, so the compiler will complain 
if an incorrect value is inadvertently used to reference that prescription.

For example, the Common_Envelope Accretion Prescriptions are declared in ``constants.h`` thus::

    enum class CE_ACCRETION_PRESCRIPTION: int { ZERO, CONSTANT, UNIFORM, MACLEOD };

    const COMPASUnorderedMap<CE_ACCRETION_PRESCRIPTION, std::string> CE_ACCRETION_PRESCRIPTION_LABEL = {
        { CE_ACCRETION_PRESCRIPTION::ZERO, ”ZERO” },
        { CE_ACCRETION_PRESCRIPTION::CONSTANT, ”CONSTANT” },
        { CE_ACCRETION_PRESCRIPTION::UNIFORM, ”UNIFORM” },
        { CE_ACCRETION_PRESCRIPTION::MACLEOD, ”MACLEOD” },
    };

Refer to ``constants.h`` for the definition of ``COMPASUnorderedMap``.

Note that the values allowed for variables of type ``CE_ACCRETION_PRESCRIPTION`` are limited to ZERO, CONSTANT, UNIFORM, and MACLEOD – 
anything else will cause a compilation error.

The unordered map ``CE_ACCRETION_PRESCRIPTION_LABEL`` declares a string label for each ``CE_ACCRETION_PRESCRIPTION``, and is indexed by 
``CE_ACCRETION_PRESCRIPTION``. The strings declared in ``CE_ACCRETION_PRESCRIPTION_LABEL`` are used by the Options service to match user
input to the required ``CE_ACCRETION_PRESCRIPTION``. These strings can also be used if an English description of the value of a variable
is required: instead of just printing an integer value that maps to a ``CE_ACCRETION_PRESCRIPTION``, the string label associated with the
prescription can be printed.


Stellar types are also declared in ``constants.h`` via an enum class and associated label map. This allows stellar types to be referenced
using symbolic names rather than an ordinal number. The stellar types enum class is ``STELLAR_TYPE``, and is declared as::

    enum class STELLAR_TYPE: int {
        MS_LTE_07,
        MS_GT_07,
        HERTZSPRUNG_GAP,
        FIRST_GIANT_BRANCH,
        CORE_HELIUM_BURNING,
        EARLY_ASYMPTOTIC_GIANT_BRANCH,
        THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH,
        NAKED_HELIUM_STAR_MS,
        NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
        NAKED_HELIUM_STAR_GIANT_BRANCH,
        HELIUM_WHITE_DWARF,
        CARBON_OXYGEN_WHITE_DWARF,
        OXYGEN_NEON_WHITE_DWARF,
        NEUTRON_STAR,
        BLACK_HOLE,
        MASSLESS_REMNANT,
        CHEMICALLY_HOMOGENEOUS,
        STAR,
        BINARY_STAR,
        NONE
    };

Ordinal numbers can still be used to reference the stellar types, and because of the order of definition in the enum class the ordinal numbers
match those given in :cite:`Hurley2000`.

The label map ``STELLAR_TYPE_LABEL`` can be used to print text descriptions of the stellar types, and is declared as::

    const std::unordered map<STELLAR_TYPE, std::string> STELLAR_TYPE_LABEL = {
        { STELLAR TYPE::MS_LTE_07,                                 ”Main Sequence <= 0.7” },
        { STELLAR_TYPE::MS_GT_07,                                  ”Main Sequence > 0.7” },
        { STELLAR_TYPE::HERTZSPRUNG_GAP,                           ”Hertzsprung Gap” },
        { STELLAR_TYPE::FIRST_GIANT_BRANCH,                        ”First Giant Branch” },
        { STELLAR_TYPE::CORE_HELIUM_BURNING,                       ”Core Helium Burning” },
        { STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,             ”Early Asymptotic Giant Branch” },
        { STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH, ”Thermally Pulsing Asymptotic Giant Branch” },
        { STELLAR_TYPE::NAKED_HELIUM_STAR_MS,                      ”Naked Helium Star MS” },
        { STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,         ”Naked Helium Star Hertzsprung Gap” },
        { STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,            ”Naked Helium Star Giant Branch” },
        { STELLAR_TYPE::HELIUM_WHITE_DWARF,                        ”Helium White Dwarf” },
        { STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,                 ”Carbon-Oxygen White Dwarf” },
        { STELLAR_TYPE::OXYGEN NEON WHITE DWARF,                   ”Oxygen-Neon White Dwarf” },
        { STELLAR_TYPE::NEUTRON_STAR,                              ”Neutron Star” },
        { STELLAR_TYPE::BLACK_HOLE,                                ”Black Hole” },
        { STELLAR_TYPE::MASSLESS_REMNANT,                          ”Massless Remnant” },
        { STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS,                    ”Chemically Homogeneous” },
        { STELLAR_TYPE::STAR,                                      ”Star” },
        { STELLAR_TYPE::BINARY_STAR,                               ”Binary Star” },
        { STELLAR_TYPE::NONE,                                      ”Not a Star!” }
    };

To print the ordinal number of a stellar type as an integer (sometimes referred to as the "Hurley type"), use ``static_cast``.  For example::

    std::cout << "CHeB ordinal number = " << static_cast<int>(STELLAR_TYPE::CORE_HELIUM_BURNING) << "\n";

To print the text label of a stellar type::

    std::cout << "CHeB label = " << STELLAR_TYPE_LABEL.at(STELLAR_TYPE::CORE_HELIUM_BURNING) << "\n";
