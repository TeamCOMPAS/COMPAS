def chirp_mass(m1, m2):
    """ Calculate the chirp mass of a binary """
    return (m1 * m2) ** (3 / 5) / (m1 + m2) ** (1 / 5)


def eta(m1, m2):
    """ Calculate the symmetric mass ratio of a binary """
    return m1 * m2 / (m1 + m2) ** 2
